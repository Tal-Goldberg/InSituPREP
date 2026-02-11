from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple, Mapping

import numpy as np
import pandas as pd
import warnings
import os
import contextlib
from scipy.spatial import cKDTree

from statsmodels.stats.multitest import multipletests
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


@dataclass
class DeseqRunParams:
    tissue: str
    dist_threshold: float = 1
    n_perm: int = 1000
    pv_thresh: float = 0.05
    fdr_thresh: float = 0.1
    lfc_thresholds: Tuple[float, ...] = (0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)
    add_pseudocount: int = 1


def read_genes_list(genes_names_path: Path) -> List[str]:
    genes = []
    with open(genes_names_path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.append(g)
    return genes


def load_marker_genes_for_tissue(
    marker_genes_by_tissue: Optional[Mapping[str, str]],
    tissue: str,
) -> Optional[pd.DataFrame]:
    """
    Load marker genes for a specific tissue from a mapping: tissue -> csv path.
    If the tissue is missing in the mapping (or path is empty/None/"None") -> return None.
    If file path exists but file is missing -> warn and return None.
    """
    if not marker_genes_by_tissue:
        return None

    path_str = marker_genes_by_tissue.get(str(tissue), None)

    if path_str is None:
        return None
    if str(path_str).strip() == "" or str(path_str).strip().lower() == "none":
        return None

    marker_path = Path(path_str)
    if not marker_path.exists():
        warnings.warn(f"Marker genes file not found for tissue {tissue}: {marker_path}. Continuing without markers.")
        return None

    mg = pd.read_csv(marker_path)
    if ("CellType" not in mg.columns) or ("Marker" not in mg.columns):
        raise ValueError(f"Marker genes file {marker_path} must contain columns: 'CellType', 'Marker'")
    return mg


def get_proximity(
    distance_matrix: pd.DataFrame,
    primary_cells: pd.Index,
    neighbor_cells: pd.Index,
    threshold: float,
) -> pd.DataFrame:
    """Return boolean proximity matrix (primary x neighbor) where True means distance <= threshold."""
    primary_cells = primary_cells.astype(str)
    neighbor_cells = neighbor_cells.astype(str)

    if not set(primary_cells).issubset(distance_matrix.index):
        missing = list(set(primary_cells) - set(distance_matrix.index))
        raise ValueError(f"Primary cell IDs missing in distance matrix index (example): {missing[:5]}")
    if not set(neighbor_cells).issubset(distance_matrix.index):
        missing = list(set(neighbor_cells) - set(distance_matrix.index))
        raise ValueError(f"Neighbor cell IDs missing in distance matrix index (example): {missing[:5]}")

    proximity = distance_matrix.loc[primary_cells, neighbor_cells]
    return proximity <= threshold

def compute_proximal_from_coordinates(
    labels_tissue: pd.DataFrame,
    primary_cells: pd.Index,
    neighbor_cells: pd.Index,
    dist_threshold: float,
    x_col: str = "X_space",
    y_col: str = "Y_space",
    z_col: str = "Z_space",
) -> pd.Series:
    """
    Return boolean Series indexed by primary_cells:
    True if min distance(primary -> any neighbor) <= dist_threshold.
    Uses 3D if z_col exists, else 2D.
    """
    primary_cells = primary_cells.astype(str)
    neighbor_cells = neighbor_cells.astype(str)

    df_ct1 = labels_tissue.loc[labels_tissue.index.intersection(primary_cells)].copy()
    df_ct2 = labels_tissue.loc[labels_tissue.index.intersection(neighbor_cells)].copy()

    if df_ct1.empty or df_ct2.empty:
        # no neighbors or no primary -> all False
        return pd.Series(False, index=primary_cells)

    # coordinates
    if z_col in labels_tissue.columns:
        cols = [x_col, y_col, z_col]
    else:
        cols = [x_col, y_col]

    # coerce numeric + drop cells with missing coords
    df_ct1 = df_ct1[cols].apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")
    df_ct2 = df_ct2[cols].apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")

    if df_ct1.empty or df_ct2.empty:
        return pd.Series(False, index=primary_cells)

    tree = cKDTree(df_ct2.values)
    dists, _ = tree.query(df_ct1.values, k=1)

    prox = pd.Series(False, index=primary_cells)
    prox.loc[df_ct1.index.astype(str)] = (dists <= float(dist_threshold))
    return prox


def validate_counts_matrix(
    counts_int: pd.DataFrame,
    meta_df: pd.DataFrame,
    distance_matrix: Optional[pd.DataFrame] = None,
) -> None:
    """
    Ensure counts_int is samples x genes as PyDESeq2 expects:
      - rows (index) = samples (cells)
      - columns = genes
    """
    if not counts_int.index.equals(meta_df.index):
        raise ValueError(
            "counts_int index must match meta_df index exactly (samples alignment). "
            f"counts_int.index[:3]={list(counts_int.index[:3])}, meta_df.index[:3]={list(meta_df.index[:3])}"
        )
    if distance_matrix is not None:
        col_as_cells_overlap = len(set(counts_int.columns.astype(str)) & set(distance_matrix.index.astype(str)))
        if col_as_cells_overlap > 0:
            overlap = list(set(counts_int.columns.astype(str)) & set(distance_matrix.index.astype(str)))[:5]
            raise ValueError(
                "counts_int columns overlap with distance_matrix cell IDs. "
                "This suggests counts_int might be transposed (genes should be columns). "
                f"Example overlapping IDs: {overlap}"
            )

    if counts_int.shape[0] < 2:
        raise ValueError("counts_int must have at least 2 samples (rows).")
    if counts_int.shape[1] < 1:
        raise ValueError("counts_int must have at least 1 gene (column).")


def _bh_adjust(pvals: pd.Series) -> pd.Series:
    adj = multipletests(pvals.values, method="fdr_bh")[1]
    return pd.Series(adj, index=pvals.index)


def filter_genes_for_threshold(
    res_df: pd.DataFrame,
    marker_genes_df: Optional[pd.DataFrame],
    neighbor_cell_type: str,
    lfc_thresh: Optional[float],
    pv_thresh: float,
    fdr_thresh: float,
) -> pd.DataFrame:
    """
    Filtering order (requested):
    1) Filter by absolute log2 fold-change threshold (if provided).
    2) Compute BH-adjusted p-values (analytic and permutation) on the remaining genes.
    3) Apply final thresholds:
       pvalue <= pv_thresh AND pvalue_perm <= pv_thresh AND
       p_value_adj_new <= fdr_thresh AND p_value_adj_new_perm <= fdr_thresh
    4) Optionally remove neighbor marker genes (if marker_genes_df is provided).
    """
    df = res_df.copy()

    # Coerce numeric columns
    df["pvalue"] = pd.to_numeric(df.get("pvalue"), errors="coerce")
    df["pvalue_perm"] = pd.to_numeric(df.get("pvalue_perm"), errors="coerce")
    df["log2FoldChange"] = pd.to_numeric(df.get("log2FoldChange"), errors="coerce")

    # Initialize output columns
    df["p_value_adj_new"] = np.nan
    df["p_value_adj_new_perm"] = np.nan

    # 1) LFC filter FIRST (abs(log2FC) >= threshold)
    if lfc_thresh is not None:
        df = df.loc[df["log2FoldChange"].notna()].copy()
        df = df.loc[df["log2FoldChange"].abs() >= float(lfc_thresh)].copy()

    if df.empty:
        return df

    # 2) Compute BH-adjusted p-values on the remaining genes
    finite_p = df["pvalue"].notna() & np.isfinite(df["pvalue"].values)
    if finite_p.any():
        df.loc[finite_p, "p_value_adj_new"] = _bh_adjust(df.loc[finite_p, "pvalue"]).values

    finite_pp = df["pvalue_perm"].notna() & np.isfinite(df["pvalue_perm"].values)
    if finite_pp.any():
        df.loc[finite_pp, "p_value_adj_new_perm"] = _bh_adjust(df.loc[finite_pp, "pvalue_perm"]).values

    # 3) Apply final thresholds (raw p AND BH-adjusted)
    df = df.loc[
        df["pvalue"].notna() &
        (df["pvalue"] <= pv_thresh) &
        df["pvalue_perm"].notna() &
        (df["pvalue_perm"] <= pv_thresh) &
        df["p_value_adj_new"].notna() &
        (df["p_value_adj_new"] <= fdr_thresh) &
        df["p_value_adj_new_perm"].notna() &
        (df["p_value_adj_new_perm"] <= fdr_thresh)
    ].copy()

    if df.empty:
        return df

    # 4) Optional: remove neighbor marker genes
    if marker_genes_df is not None:
        markers = (
            marker_genes_df.loc[marker_genes_df["CellType"] == neighbor_cell_type, "Marker"]
            .astype(str)
        )
        df = df.loc[~df["gene"].astype(str).isin(set(markers))].copy()

    return df


def run_deseq_pairs(
    summary_table_path: Path,
    distance_matrix_path: Optional[Path],  
    genes_names_path: Path,
    out_dir: Path,
    params: DeseqRunParams,
    marker_genes_by_tissue: Optional[Mapping[str, str]] = None,
    primary_cell_types: Optional[List[str]] = None,
    neighbor_cell_types: Optional[List[str]] = None,
    quiet_permutations: bool = True,
    perm_print_every: int = 50,
) -> Dict[str, Any]:
    """
    Run DESeq-like analysis over cell-type pairs for a given tissue.
    If this tissue has no marker genes path in marker_genes_by_tissue, analysis runs normally
    and marker filtering is skipped.
    """
    if params.tissue is None or str(params.tissue).strip() == "":
        raise ValueError(
            "Parameter 'tissue' must be provided. "
            "Running on all tissues is not supported in this analysis."
        )

    out_dir.mkdir(parents=True, exist_ok=True)
    genes_names = read_genes_list(genes_names_path)

    marker_genes_df = load_marker_genes_for_tissue(marker_genes_by_tissue, tissue=str(params.tissue))

    labels_full = pd.read_csv(summary_table_path, header=0, dtype={"Var1": str})
    if "tissue" in labels_full.columns:
        labels_full["tissue"] = labels_full["tissue"].astype(str)

    if "Var1" not in labels_full.columns:
        raise ValueError(
            f"summary_table_path must contain a 'Var1' column with unique cell IDs. "
            f"Missing 'Var1' in: {summary_table_path}"
        )

    # set Var1 as the index explicitly
    labels_full["Var1"] = labels_full["Var1"].astype(str)
    labels_full = labels_full.set_index("Var1", drop=True)

    # optional: ensure uniqueness
    if not labels_full.index.is_unique:
        dup = labels_full.index[labels_full.index.duplicated()].astype(str)[:5].tolist()
        raise ValueError(f"'Var1' must contain unique cell IDs. Example duplicates: {dup}")

    labels_full.index = labels_full.index.astype(str)

    
    labels_full["tissue"] = labels_full["tissue"].astype(str)
    labels = labels_full[labels_full["tissue"] == str(params.tissue)].copy()
    labels.index = labels.index.astype(str)
    
    distance_matrix = None
    if distance_matrix_path is not None:
        if not Path(distance_matrix_path).exists():
            raise ValueError(f"distance_matrix_path does not exist: {distance_matrix_path}")
        
        distance_matrix = pd.read_csv(distance_matrix_path, header=0, index_col=0, dtype={0: str})
        distance_matrix.index = distance_matrix.index.astype(str)
        distance_matrix.columns = distance_matrix.columns.astype(str)
        print("[INFO] Using distance matrix for proximity.")
    else:
        print("[INFO] No distance matrix provided -> will use coordinate-based proximity.")
    

    cell_types_all = labels["cell_type"].dropna().unique().tolist()
    primary_list = primary_cell_types if primary_cell_types is not None else cell_types_all
    neighbor_list = neighbor_cell_types if neighbor_cell_types is not None else cell_types_all

    genes_missing = [g for g in genes_names if g not in labels.columns]
    if genes_missing:
        warnings.warn(
            f"{len(genes_missing)} genes from genes_names were not found in labels columns. Using intersection only."
        )

    common_genes = [g for g in genes_names if g in labels.columns]
    if len(common_genes) == 0:
        raise ValueError("No overlap between genes_names and labels gene columns.")

    results_index = []

    for primary_ct in primary_list:
        primary_cells = labels.index[labels["cell_type"] == primary_ct].astype(str)

        raw_data = labels.loc[primary_cells, common_genes].copy()
        raw_data = raw_data.apply(pd.to_numeric, errors="coerce")
        raw_data = raw_data.dropna(axis=0, how="any")
        if raw_data.shape[0] < 2:
            warnings.warn(f"[SKIP] Not enough primary cells for {primary_ct} (n={raw_data.shape[0]}).")
            continue

        for neighbor_ct in neighbor_list:
            if primary_ct == neighbor_ct:
                continue

            # print which pair is running
            print(f"\n[PAIR] Tissue {params.tissue} | Primary: {primary_ct} | Neighbor: {neighbor_ct}")

            neighbor_cells = labels.index[labels["cell_type"] == neighbor_ct].astype(str)

            if distance_matrix is not None:
                proximity = get_proximity(distance_matrix, raw_data.index, neighbor_cells, params.dist_threshold)
                if proximity.shape[1] == 0:
                    warnings.warn(f"[SKIP] No neighbor cells for {neighbor_ct}.")
                    continue
                prox_any = proximity.any(axis=1).reindex(raw_data.index).fillna(False).astype(bool)

            else:
                # Fallback: compute proximal from coordinates
                needed = {"X_space", "Y_space"}
                if not needed.issubset(labels.columns):
                    raise ValueError(
                        "distance_matrix_path not provided, so labels must contain "
                        "X_space and Y_space columns.")

                labels_tissue = labels  # per tissue
                prox_any = compute_proximal_from_coordinates(
                    labels_tissue=labels_tissue,
                    primary_cells=raw_data.index,
                    neighbor_cells=neighbor_cells,
                    dist_threshold=params.dist_threshold,
                ).reindex(raw_data.index).fillna(False).astype(bool)

            # prox/nonprox
            if prox_any.astype(int).nunique() < 2:
                warnings.warn(f"[SKIP] Only one proximity state for {primary_ct} vs {neighbor_ct}.")
                continue


            counts = raw_data.loc[prox_any.index].copy()
            counts = counts.dropna(axis=1, how="any")
            if counts.shape[1] == 0:
                warnings.warn(f"[SKIP] No genes left after NaN filtering for {primary_ct} vs {neighbor_ct}.")
                continue

            counts_int = counts.round().astype(int) + params.add_pseudocount
            prox_any_aligned = prox_any.reindex(counts_int.index).fillna(False).astype(bool)
            meta_df = pd.DataFrame({"proximal": prox_any_aligned.astype(int)}, index=counts_int.index)
            meta_df["proximal"] = meta_df["proximal"].map({0: "nonprox", 1: "prox"}).astype("category")

            validate_counts_matrix(counts_int=counts_int, meta_df=meta_df, distance_matrix=distance_matrix)

            try:
                dds = DeseqDataSet(counts=counts_int, metadata=meta_df, design="~ proximal")
                dds.deseq2()
                ds = DeseqStats(dds, contrast=("proximal", "prox", "nonprox"))
                ds.summary()
                res_df = ds.results_df.copy().reset_index().rename(columns={"index": "gene"})
            except Exception as e:
                warnings.warn(f"[SKIP] PyDESeq2 failed for {primary_ct} vs {neighbor_ct}: {e}")
                continue

            # ---- permutations ----
            p_orig = pd.to_numeric(res_df["pvalue"], errors="coerce").fillna(1.0).astype(float).values
            genes_order = res_df["gene"].astype(str).values
            sum_perm = np.zeros(len(res_df), dtype=np.int64)
            n_perm_done = 0

            for perm_i in range(params.n_perm):
                # ✅ print permutation progress every N
                if perm_print_every > 0 and (perm_i + 1) % perm_print_every == 0:
                    print(
                        f"[PERM] Tissue {params.tissue} | {primary_ct} vs {neighbor_ct} | "
                        f"{perm_i + 1}/{params.n_perm}"
                    )

                shuffled = meta_df["proximal"].sample(frac=1.0, replace=False).reset_index(drop=True)
                meta_perm = meta_df.copy()
                meta_perm["proximal"] = shuffled.values
                meta_perm["proximal"] = meta_perm["proximal"].astype("category")

                try:
                    if quiet_permutations:
                        with open(os.devnull, "w") as f, \
                                contextlib.redirect_stdout(f), \
                                contextlib.redirect_stderr(f), \
                                warnings.catch_warnings():
                            warnings.simplefilter("ignore")
                            dds_perm = DeseqDataSet(counts=counts_int, metadata=meta_perm, design="~ proximal")
                            dds_perm.deseq2()
                            ds_perm = DeseqStats(dds_perm, contrast=("proximal", "prox", "nonprox"))
                            ds_perm.summary()
                    else:
                        dds_perm = DeseqDataSet(counts=counts_int, metadata=meta_perm, design="~ proximal")
                        dds_perm.deseq2()
                        ds_perm = DeseqStats(dds_perm, contrast=("proximal", "prox", "nonprox"))
                        ds_perm.summary()

                    res_perm = ds_perm.results_df.copy().reset_index().rename(columns={"index": "gene"})
                    res_perm["gene"] = res_perm["gene"].astype(str)
                    res_perm = res_perm.set_index("gene").reindex(genes_order)

                    p_perm = pd.to_numeric(res_perm["pvalue"], errors="coerce").fillna(1.0).astype(float).values
                    sum_perm += (p_perm <= p_orig).astype(np.int64)
                    n_perm_done += 1

                except Exception as e:
                    warnings.warn(
                        f"[PERM-SKIP] perm {perm_i + 1}/{params.n_perm} failed for {primary_ct} vs {neighbor_ct}: {e}"
                    )

            if n_perm_done == 0:
                warnings.warn(f"[PERM-SKIP] No successful permutations for {primary_ct} vs {neighbor_ct}.")
                res_df["pvalue_perm"] = np.nan
            else:
                res_df["pvalue_perm"] = (sum_perm + 1) / n_perm_done

            # ---- choose threshold + filter ----
            pv_thresh = params.pv_thresh
            fdr_thresh = params.fdr_thresh

            baseline_mask = (
                (pd.to_numeric(res_df["pvalue"], errors="coerce") <= pv_thresh) &
                (pd.to_numeric(res_df["pvalue_perm"], errors="coerce") <= pv_thresh)
            )
            baseline_n = int(baseline_mask.sum())

            counts_by_thr = {}
            for thr in params.lfc_thresholds:
                df_thr = filter_genes_for_threshold(
                    res_df=res_df,
                    marker_genes_df=marker_genes_df,
                    neighbor_cell_type=neighbor_ct,
                    lfc_thresh=thr,
                    pv_thresh=pv_thresh,
                    fdr_thresh=fdr_thresh,
                )
                counts_by_thr[thr] = int(df_thr.shape[0])

            diffs = {thr: baseline_n - counts_by_thr[thr] for thr in params.lfc_thresholds}
            chosen_thr = min(params.lfc_thresholds, key=lambda t: diffs[t])

            final_df = filter_genes_for_threshold(
                res_df=res_df,
                marker_genes_df=marker_genes_df,
                neighbor_cell_type=neighbor_ct,
                lfc_thresh=chosen_thr,
                pv_thresh=pv_thresh,
                fdr_thresh=fdr_thresh,
            )

            final_df["chosen_lfc_threshold"] = chosen_thr
            final_df["primary_cell_type"] = primary_ct
            final_df["neighbor_cell_type"] = neighbor_ct
            final_df["tissue"] = str(params.tissue)

            out_path = out_dir / (
                f"final_genes_{params.tissue}_{primary_ct}_vs_{neighbor_ct}_LFC_{str(chosen_thr).replace('.','')}.csv"
            )
            final_df.to_csv(out_path, index=False)

            n_final_genes = int(final_df.shape[0])

            print(
                f"[DONE] Tissue {params.tissue} | {primary_ct} vs {neighbor_ct} | "
                f"Final genes: {n_final_genes} | Permutations: {n_perm_done}/{params.n_perm}"
            )

            results_index.append({
                "tissue": params.tissue,
                "primary": primary_ct,
                "neighbor": neighbor_ct,
                "baseline_n": baseline_n,
                "chosen_lfc": chosen_thr,
                "n_final_genes": n_final_genes,
                "out_csv": str(out_path),
                "n_perm_done": n_perm_done,
                "marker_genes_used": bool(marker_genes_df is not None),
            })

    summary_df = pd.DataFrame(results_index)
    summary_path = out_dir / f"deseq_summary_tissue_{params.tissue}.csv"
    summary_df.to_csv(summary_path, index=False)

    return {"summary_csv": str(summary_path), "n_pairs_done": int(summary_df.shape[0])}
