#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN

try:
    from scipy.spatial.distance import pdist
    from scipy.spatial import cKDTree
    from scipy.stats import norm
except Exception as e:
    raise ImportError("This script requires scipy. Please install it (pip install scipy).") from e


# ======================================================================================
# Params
# ======================================================================================

@dataclass(frozen=True)
class BacteriaObjectsParams:
    # Transcripts / object detection (GLOBAL)
    input_csv: Path
    out_dir: Path

    target_gene: str
    gene_col: str

    tx_x_col: str
    tx_y_col: str
    tx_z_col: Optional[str]  # None => 2D (X,Y only)

    tx_cell_id_col: str  # transcript cell_id string (empty/NaN if not inside a cell)

    eps_min: int
    eps_max: int
    min_samples: int
    n_perms: int
    alpha: float
    eps_final: Optional[int]  # None => auto by global max F1 across eps
    seed: Optional[int]  # None => random seed

    # NEW: resume/reuse DBSCAN outputs if present
    reuse_existing_objects: bool

    # DESeq2 (optional; GLOBAL)
    summary_table_path: Optional[Path]

    celltype_col: str

    cell_x_col: str
    cell_y_col: str
    cell_z_col: Optional[str]  # None => 2D (X,Y only) for cells


    dist_threshold: float
    min_cells_per_group: int

    # explicit gene list for DESeq2
    genes_names_path: Optional[Path]


# ======================================================================================
# Core utilities
# ======================================================================================

def sum_pairwise_distances(coords: np.ndarray) -> float:
    """Sum of all pairwise Euclidean distances within coords (N x 3)."""
    if coords.shape[0] < 2:
        return 0.0
    return float(pdist(coords, metric="euclidean").sum())


def compute_group_pvalues(
    coords_target: np.ndarray,
    labels: np.ndarray,
    coords_null_pool: np.ndarray,
    n_perms: int,
    rng: np.random.Generator,
) -> Dict[int, float]:
    """
    Per cluster id (>=1), compute left-tail p-value:
        z = (D_obs - mean(D_perm)) / std(D_perm)
        p = norm.cdf(z)
    where D = sum pairwise distances in the cluster.
    Null perms sample same-size sets from coords_null_pool.
    """
    cluster_ids = sorted([cid for cid in np.unique(labels) if cid >= 1])
    pvals: Dict[int, float] = {}

    if len(cluster_ids) == 0:
        return pvals

    if coords_null_pool.shape[0] < 2:
        for cid in cluster_ids:
            pvals[cid] = float("nan")
        return pvals

    n_pool = coords_null_pool.shape[0]

    for cid in cluster_ids:
        idx = np.where(labels == cid)[0]
        m = idx.size

        if m < 2:
            pvals[cid] = 1.0
            continue

        d_obs = sum_pairwise_distances(coords_target[idx, :])

        replace = m > n_pool
        d_perm = np.empty(n_perms, dtype=float)
        for p in range(n_perms):
            pick = rng.choice(n_pool, size=m, replace=replace)
            d_perm[p] = sum_pairwise_distances(coords_null_pool[pick, :])

        mu = float(d_perm.mean())
        sd = float(d_perm.std(ddof=0))

        if sd == 0.0 or np.isnan(sd):
            pvals[cid] = float("nan")
        else:
            z = (d_obs - mu) / sd
            pvals[cid] = float(norm.cdf(z))  # left-tail

    return pvals


def most_frequent_gene_pool(df_all: pd.DataFrame, gene_col: str, xyz_cols: List[str]) -> Tuple[str, np.ndarray]:
    """Return (most frequent gene name, coords pool for that gene)."""
    counts = df_all[gene_col].value_counts(dropna=False)
    most_gene = str(counts.index[0])
    pool = df_all.loc[df_all[gene_col] == most_gene, xyz_cols].to_numpy(dtype=float)
    return most_gene, pool


def sanitize_filename(s: str) -> str:
    s = s.strip()
    s = re.sub(r"[^\w\-_\. ]+", "_", s)
    s = s.replace(" ", "_")
    return s


def normalize_cell_id_series(s: pd.Series) -> pd.Series:
    """Normalize cell_id strings: strip; empty -> NA; keep NA."""
    out = s.astype("string")
    out = out.str.strip()
    out = out.replace({"": pd.NA})
    return out


def load_gene_list(path: Path) -> List[str]:
    genes: List[str] = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if not g:
                continue
            genes.append(g)
    # keep order, drop duplicates
    seen = set()
    uniq: List[str] = []
    for g in genes:
        if g not in seen:
            uniq.append(g)
            seen.add(g)
    return uniq


def select_gene_columns_from_list(cells_df: pd.DataFrame, genes: List[str]) -> List[str]:
    present = [g for g in genes if g in cells_df.columns]
    missing = [g for g in genes if g not in cells_df.columns]
    if len(present) == 0:
        raise ValueError("None of the genes from genes_names.txt were found as columns in the cells summary CSV.")
    if missing:
        preview = ", ".join(missing[:20])
        raise ValueError(
            "Some genes from genes_names.txt are missing from the cells summary CSV columns "
            f"(showing up to 20): {preview}"
        )
    return present


# ======================================================================================
# DESeq2 step (GLOBAL) - GROUPS: near_sig vs far_sig
# ======================================================================================

def build_cell_groups_global(
    transcripts_df: pd.DataFrame,
    cells_df: pd.DataFrame,
    params: BacteriaObjectsParams,
) -> pd.DataFrame:
    """
    Adds 'bacteria_group' to cells_df:
      - near_sig: min distance to ANY significant target transcript <= dist_threshold
                  OR contains >=1 significant target transcript
      - far_sig : min distance to significant target transcripts > dist_threshold

    Significant target transcript:
      gene == target_gene AND group_id>=1 AND p_value_group<=alpha
    """
    out = cells_df.copy()
    out["bacteria_group"] = np.nan


    # Cell IDs are always taken from the dataframe index (Var1)
    out.index = out.index.astype("string").str.strip()

    # Normalize transcripts cell_id
    tx = transcripts_df.copy()
    tx[params.tx_cell_id_col] = normalize_cell_id_series(tx[params.tx_cell_id_col])

    # Significant target transcripts (global)
    sig_tx = tx.loc[
        (tx[params.gene_col] == params.target_gene)
        & (tx["group_id"].notna())
        & (tx["group_id"] >= 1)
        & (tx["p_value_group"].notna())
        & (tx["p_value_group"] <= params.alpha)
    ].copy()

    if sig_tx.empty:
        print("[DESeq2] No significant target transcripts found (p<=alpha & clustered). Skipping grouping.")
        return out

    # Ensure transcripts/cells use the same dimensionality for KDTree queries
    tx_dim = 2 if params.tx_z_col is None else 3
    cell_dim = 2 if params.cell_z_col is None else 3
    if tx_dim != cell_dim:
        raise ValueError(
            f"Dimension mismatch: transcripts are {tx_dim}D but cells are {cell_dim}D. "
            "Either provide both Z columns or set both to 2D."
        )

       
    # KDTree from significant transcript coords (2D or 3D depending on tx_z_col)
    tx_cols = [params.tx_x_col, params.tx_y_col]
    if params.tx_z_col is not None:
        tx_cols.append(params.tx_z_col)

    tx_coords = sig_tx[tx_cols].to_numpy(dtype=float)
    tree = cKDTree(tx_coords)


    # Only label "real" cells: index notna
    real_cells = out.index.notna()
    if not real_cells.any():
        print("[DESeq2] No real cells found (index is NA for all). Skipping grouping.")
        return out

    cell_cols = [params.cell_x_col, params.cell_y_col]
    if params.cell_z_col is not None:
        cell_cols.append(params.cell_z_col)

    cell_coords = out.loc[real_cells, cell_cols].to_numpy(dtype=float)
    dists, _ = tree.query(cell_coords, k=1)

    near_by_dist = dists <= float(params.dist_threshold)
    far_by_dist = dists > float(params.dist_threshold)

    real_idx = out.loc[real_cells].index
    out.loc[real_idx[near_by_dist], "bacteria_group"] = "near_sig"
    out.loc[real_idx[far_by_dist], "bacteria_group"] = "far_sig"

    # Ensure: if cell contains a significant transcript => near_sig
    sig_tx_in_cells = sig_tx.loc[sig_tx[params.tx_cell_id_col].notna()].copy()
    if not sig_tx_in_cells.empty:
        contains_ids = set(sig_tx_in_cells[params.tx_cell_id_col].astype(str).unique().tolist())
        in_contains = out.index.astype(str).isin(contains_ids)
        out.loc[in_contains, "bacteria_group"] = "near_sig"

    n_near = int((out["bacteria_group"] == "near_sig").sum())
    n_far = int((out["bacteria_group"] == "far_sig").sum())
    print(f"[DESeq2] Grouping done: near_sig={n_near:,} | far_sig={n_far:,} | dist_threshold={params.dist_threshold}")

    return out


def run_deseq2_for_celltype(counts_df: pd.DataFrame, meta_df: pd.DataFrame) -> pd.DataFrame:
    """
    Run pyDESeq2 comparing:
      near_sig  vs  far_sig
    with far_sig as reference.
    """
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except Exception as e:
        raise ImportError("pyDESeq2 is required for this step. Install with: pip install pydeseq2") from e

    # IMPORTANT: +1 pseudocount
    counts = counts_df.round().astype(int) + 1

    meta_df = meta_df.copy()

    # FORCE SAME INDEX TYPE + SAME ORDER
    counts.index = counts.index.astype(str)
    meta_df.index = meta_df.index.astype(str)
    meta_df = meta_df.loc[counts.index]

    assert counts.index.equals(meta_df.index), "counts/meta index mismatch right before DESeq2"

    print(counts.index.dtype, meta_df.index.dtype)
    print(counts.index[:5].tolist())
    print(meta_df.index[:5].tolist())

    dds = DeseqDataSet(
        counts=counts,
        metadata=meta_df,
        design_factors=["bacteria_group"],
        ref_level=["bacteria_group", "far_sig"],
        n_cpus=1,
    )
    dds.deseq2()

    # FIX: handle newer/older pydeseq2 APIs
    try:
        stat_res = DeseqStats(dds, contrast=("bacteria_group", "near_sig", "far_sig"))
    except TypeError:
        # older pydeseq2 API
        stat_res = DeseqStats(dds)

    stat_res.summary()

    res = stat_res.results_df.copy()
    res = res.reset_index().rename(columns={"index": "gene"})
    return res


# ======================================================================================
# Main runner
# ======================================================================================

def run_bacteria_objects(params: BacteriaObjectsParams) -> None:
    rng = np.random.default_rng(params.seed)

    # Resolve output paths (single out_dir)
    out_dir = params.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    out_transcripts_csv = out_dir / "transcripts_with_objects.csv"
    out_eps_summary_csv = out_dir / "eps_sweep_summary.csv"
    out_deseq_dir = out_dir / "deseq2_results"

    print("[1/4] Loading transcripts...")
    tx = pd.read_csv(params.input_csv)

    required_tx = {params.gene_col, params.tx_x_col, params.tx_y_col}
    if params.tx_z_col is not None:
        required_tx.add(params.tx_z_col)

    missing_tx = [c for c in required_tx if c not in tx.columns]
    if missing_tx:
        raise ValueError(f"Missing required columns in transcripts CSV: {missing_tx}")

    tx = tx.copy()
    tx[params.tx_cell_id_col] = normalize_cell_id_series(tx[params.tx_cell_id_col])

    if "transcript_id" in tx.columns:
        tx = tx.sort_values(["transcript_id"], kind="mergesort").reset_index(drop=True)

    xy_cols = [params.tx_x_col, params.tx_y_col]

    if params.tx_z_col is not None:
        xyz_cols = xy_cols + [params.tx_z_col]
        ndim = 3
    else:
        xyz_cols = xy_cols
        ndim = 2

    print(f"[INFO] Running analysis in {ndim}D using coordinates: {', '.join(xyz_cols)}")


    # If asked to reuse existing DBSCAN outputs, load them and skip steps [2/4]-[3/4]
    if params.reuse_existing_objects:
        if not out_transcripts_csv.exists() or not out_eps_summary_csv.exists():
            raise FileNotFoundError(
                "reuse_existing_objects=True but required outputs are missing in out_dir: "
                f"{out_transcripts_csv} and/or {out_eps_summary_csv}"
            )

        print("[2/4] Reusing existing DBSCAN outputs (skipping DBSCAN + permutations)...")
        tx_out = pd.read_csv(out_transcripts_csv)
        summary_df = pd.read_csv(out_eps_summary_csv)

        # Basic sanity checks
        for col in ["group_id", "p_value_group"]:
            if col not in tx_out.columns:
                raise ValueError(f"Expected column '{col}' in {out_transcripts_csv} but it is missing.")

        # Determine eps_selected from summary if present; otherwise keep None
        eps_selected = None
        if "eps_selected" in summary_df.columns and summary_df["eps_selected"].notna().any():
            eps_selected = int(summary_df.loc[summary_df["eps_selected"].notna(), "eps_selected"].iloc[0])
        elif "eps_selected_global" in summary_df.columns and summary_df["eps_selected_global"].notna().any():
            eps_selected = int(summary_df.loc[summary_df["eps_selected_global"].notna(), "eps_selected_global"].iloc[0])

        # For DESeq2 downstream, we still need params.target_gene and params.alpha (already in params).
        tx_target = tx_out.loc[tx_out[params.gene_col] == params.target_gene].copy()
        n_target = int(len(tx_target))

        print(f"[2/4] Loaded existing outputs: n={len(tx_out):,} | target='{params.target_gene}' n_target={n_target:,}")
        if eps_selected is not None:
            print(f"[2/4] eps_selected loaded from summary: {eps_selected}")
        else:
            print("[2/4] eps_selected not found in summary; continuing without it (DESeq2 does not require eps).")

        # Skip directly to DESeq2 step (which is [4/4] in logs)
        print("[3/4] Skipped (reuse_existing_objects=True).")
        # Set tx_out variable for DESeq2 usage below
    else:
        tx_target = tx.loc[tx[params.gene_col] == params.target_gene].copy()
        n_target = int(len(tx_target))
        most_gene, null_pool = most_frequent_gene_pool(tx, params.gene_col, xyz_cols)

        print(f"[1/4] Loaded transcripts: n={len(tx):,} | target='{params.target_gene}' n_target={n_target:,}")
        print(f"[1/4] Output directory: {out_dir}")

        # ---------- Sweep eps globally ----------
        print("[2/4] Sweeping eps (DBSCAN + permutations)...")
        eps_stats: Dict[int, Dict[str, float]] = {}
        eps_aux: Dict[int, Dict[str, float]] = {}  # NEW: store informative counts per eps

        if n_target == 0:
            for eps in range(params.eps_min, params.eps_max + 1):
                eps_stats[eps] = {"precision": 0.0, "recall": 0.0, "f1": 0.0}
                eps_aux[eps] = {
                    "n_target_transcripts": 0.0,
                    "n_target_in_clusters": 0.0,
                    "n_clusters_total": 0.0,
                    "n_clusters_significant": 0.0,
                }
        else:
            coords_target = tx_target[xyz_cols].to_numpy(dtype=float)

            for eps in range(params.eps_min, params.eps_max + 1):
                labels0 = DBSCAN(eps=float(eps), min_samples=params.min_samples).fit_predict(coords_target)
                labels = labels0.copy()
                labels[labels0 >= 0] = labels0[labels0 >= 0] + 1  # 1..k, -1 noise

                n_target_in_clusters = int(np.sum(labels != -1))
                recall = float(n_target_in_clusters / n_target)

                cluster_ids = sorted([cid for cid in np.unique(labels) if cid >= 1])
                n_clusters_total = int(len(cluster_ids))

                n_clusters_significant = 0
                if n_clusters_total > 0:
                    pvals = compute_group_pvalues(coords_target, labels, null_pool, params.n_perms, rng)
                    pvals_list = np.array([pvals[cid] for cid in cluster_ids], dtype=float)

                    if pvals_list.size == 0 or np.all(np.isnan(pvals_list)):
                        precision = 0.0
                        n_clusters_significant = 0
                    else:
                        n_clusters_significant = int(np.sum(pvals_list <= params.alpha))
                        precision = float(n_clusters_significant / n_clusters_total)
                else:
                    precision = 0.0
                    n_clusters_significant = 0

                f1 = 0.0 if (precision + recall) == 0 else float(2 * (precision * recall) / (precision + recall))
                eps_stats[eps] = {"precision": precision, "recall": recall, "f1": f1}
                eps_aux[eps] = {
                    "n_target_transcripts": float(n_target),
                    "n_target_in_clusters": float(n_target_in_clusters),
                    "n_clusters_total": float(n_clusters_total),
                    "n_clusters_significant": float(n_clusters_significant),
                }

                # light progress
                if eps in (params.eps_min, params.eps_max) or (eps - params.eps_min) % 25 == 0:
                    print(f"  eps={eps:>3d} | precision={precision:.3f} | recall={recall:.3f} | f1={f1:.3f}")

        # ---------- Choose eps_selected (global max F1) ----------
        if params.eps_final is None:
            best_f1 = max(v["f1"] for v in eps_stats.values()) if eps_stats else 0.0
            eps_selected = min(eps for eps, v in eps_stats.items() if v["f1"] == best_f1) if eps_stats else params.eps_min
            eps_mode = "auto (max F1)"
        else:
            eps_selected = int(params.eps_final)
            eps_mode = "user"

        print(f"[2/4] eps selection done: eps_selected={eps_selected} ({eps_mode})")

        # ---------- Final assignment at eps_selected ----------
        print("[3/4] Assigning final clusters and p-values to transcripts...")
        tx_out = tx.copy()
        tx_out["group_id"] = np.nan
        tx_out["p_value_group"] = np.nan

        if n_target > 0:
            coords_target = tx_target[xyz_cols].to_numpy(dtype=float)

            labels0 = DBSCAN(eps=float(eps_selected), min_samples=params.min_samples).fit_predict(coords_target)
            labels = labels0.copy()
            labels[labels0 >= 0] = labels0[labels0 >= 0] + 1

            cluster_ids = [cid for cid in np.unique(labels) if cid >= 1]
            pvals_final: Dict[int, float] = {}
            if len(cluster_ids) > 0:
                pvals_final = compute_group_pvalues(coords_target, labels, null_pool, params.n_perms, rng)

            target_idx = tx_out.index[tx_out[params.gene_col] == params.target_gene].to_numpy()
            tx_out.loc[target_idx, "group_id"] = labels.astype(float)

            pval_per_row = np.full(labels.shape[0], np.nan, dtype=float)
            for i, cid in enumerate(labels):
                if cid >= 1:
                    pval_per_row[i] = pvals_final.get(int(cid), np.nan)
            tx_out.loc[target_idx, "p_value_group"] = pval_per_row

            n_sig_tx = int(
                (
                    (tx_out[params.gene_col] == params.target_gene)
                    & (tx_out["group_id"].notna()) & (tx_out["group_id"] >= 1)
                    & (tx_out["p_value_group"].notna()) & (tx_out["p_value_group"] <= params.alpha)
                ).sum()
            )
            print(f"[3/4] Final assignment done: significant target transcripts (p<=alpha) = {n_sig_tx:,}")
        else:
            print("[3/4] No target transcripts found; output will contain empty group_id/p_value_group columns.")

        # ---------- Write eps sweep summary (with requested columns) ----------
        rows = []
        for eps in range(params.eps_min, params.eps_max + 1):
            rows.append(
                {
                    "eps": eps,
                    "n_target_transcripts": float(eps_aux[eps]["n_target_transcripts"]) if "eps_aux" in locals() else float(n_target),
                    "n_target_in_clusters": float(eps_aux[eps]["n_target_in_clusters"]) if "eps_aux" in locals() else 0.0,
                    "recall": eps_stats[eps]["recall"],
                    "n_clusters_total": float(eps_aux[eps]["n_clusters_total"]) if "eps_aux" in locals() else 0.0,
                    "n_clusters_significant": float(eps_aux[eps]["n_clusters_significant"]) if "eps_aux" in locals() else 0.0,
                    "precision": eps_stats[eps]["precision"],
                    "f1": eps_stats[eps]["f1"],
                    "eps_selected": float(eps_selected),
                }
            )
        summary_df = pd.DataFrame(rows)
        summary_df.to_csv(out_eps_summary_csv, index=False)
        print(f"[3/4] Wrote eps sweep summary: {out_eps_summary_csv}")

        # ---------- Write transcripts output ----------
        tx_out.to_csv(out_transcripts_csv, index=False)
        print(f"[3/4] Wrote transcripts output: {out_transcripts_csv}")

    # ---------- Optional DESeq2 ----------
    if params.summary_table_path is None:
        print("[4/4] DESeq2 step skipped (summary table not provided).")
        print("DONE")
        return

    if params.genes_names_path is None:
        raise ValueError("genes_names_path must be provided when running DESeq2 (use --genes-names-path).")

    print("[4/4] Running DESeq2 per cell type (near_sig vs far_sig)...")

    # Load cells
    # require a 'Var1' column that we explicitly set as index)
    labels_full = pd.read_csv(params.summary_table_path, dtype={"Var1": str})
    if "tissue" in labels_full.columns:
        labels_full["tissue"] = labels_full["tissue"].astype(str)

    # Validate presence of Var1 column and set as index explicitly (as requested)
    if "Var1" not in labels_full.columns:
        raise ValueError(
            f"summary_table_path must contain a 'Var1' column with unique cell IDs. Missing 'Var1' in: {params.summary_table_path}"
        )
    # set Var1 as the index explicitly, ensure string dtype
    labels_full["Var1"] = labels_full["Var1"].astype(str)
    labels_full = labels_full.set_index("Var1", drop=True)

    # Ensure tissue column is treated as string (user requested)
    if "tissue" in labels_full.columns:
        labels_full["tissue"] = labels_full["tissue"].astype(str)

    cells = labels_full.copy()


    required_cells = {params.celltype_col, params.cell_x_col, params.cell_y_col}
    if params.cell_z_col is not None:
        required_cells.add(params.cell_z_col)

    missing_cells = [c for c in required_cells if c not in cells.columns]
    if missing_cells:
        raise ValueError(f"Missing required columns in summary table CSV: {missing_cells}")

    # Build groups
    cells_with_groups = build_cell_groups_global(
        transcripts_df=tx_out,
        cells_df=cells,
        params=params,
    )

    cells_filtered = cells_with_groups.loc[cells_with_groups["bacteria_group"].isin(["near_sig", "far_sig"])].copy()
    if cells_filtered.empty:
        print("[4/4] No cells assigned to near_sig / far_sig. Skipping DESeq2.")
        print("DONE")
        return

    # Load gene list and select gene columns
    genes = load_gene_list(params.genes_names_path)
    print(f"[4/4] Loaded gene list: {len(genes):,} genes from {params.genes_names_path}")

    gene_cols = select_gene_columns_from_list(cells_filtered, genes)
    print(f"[4/4] Verified gene columns present in summary table CSV: {len(gene_cols):,}")

    # Output dir for DESeq2
    out_deseq_dir.mkdir(parents=True, exist_ok=True)
    print(f"[4/4] DESeq2 outputs will be written to: {out_deseq_dir}")

    n_ran = 0
    n_skipped = 0
    for ct, df_ct in cells_filtered.groupby(params.celltype_col, sort=False):
        n_near = int((df_ct["bacteria_group"] == "near_sig").sum())
        n_far = int((df_ct["bacteria_group"] == "far_sig").sum())

        if n_near < params.min_cells_per_group or n_far < params.min_cells_per_group:
            n_skipped += 1
            continue

        print(f"  [DESeq2] cell_type='{ct}' | near_sig={n_near:,} | far_sig={n_far:,} -> running...")

        counts_df = df_ct[gene_cols].copy()
        meta_df = df_ct[["bacteria_group"]].copy()

        print(f"\n[DEBUG DESeq2] cell_type = {ct}")
        print("[DEBUG DESeq2] bacteria_group counts:")
        print(meta_df["bacteria_group"].value_counts(dropna=False))

        print("[DEBUG DESeq2] counts_df shape:", counts_df.shape)
        print("[DEBUG DESeq2] meta_df shape:", meta_df.shape)
        print("[DEBUG DESeq2] index identical:", counts_df.index.equals(meta_df.index))

        res_df = run_deseq2_for_celltype(counts_df=counts_df, meta_df=meta_df)

        out_name = f"deseq2_{sanitize_filename(str(ct))}.csv"
        out_path = out_deseq_dir / out_name
        res_df.to_csv(out_path, index=False)
        n_ran += 1

    print(f"[4/4] DESeq2 done: wrote {n_ran} result files to: {out_deseq_dir} (skipped {n_skipped} cell types due to min_cells_per_group)")
    print("DONE")


