# insituprep/analyses/RNA_velocity/rna_velocity_stage2.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import norm
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

from scipy.spatial import cKDTree

try:
    from statsmodels.stats.multitest import multipletests
except Exception:  # pragma: no cover
    multipletests = None


ProximityMode = Literal["continuous", "binary"]
YFeature = Literal["magnitude_delta", "magnitude_s0", "phase_deg_delta", "phase_deg_s0"]


@dataclass
class VelocityStage2Params:
    t: float = 3.0  # t to analyze (should match Stage1 )
    n_pcs: int = 3  # PCA dims used internally to project s(t)

    # columns in pc-origin CSV
    pc1_col: str = "PC1"
    pc2_col: str = "PC2"
    pc3_col: str = "PC3"  # optional
    cell_type_col: str = "cell_type"

    # tissue handling
    parse_tissue_from_index: bool = True
    tissue_col: str = "tissue"

    # proximity definition
    distance_threshold: float = 1.0  # µm; used for binary "close" definition
    maximum_distance_threshold: float = np.inf #maximum distance to take into account for calculating distances between cells - usually fits the size of field of view 
    proximity_mode: ProximityMode = "continuous"  # which prox_df to analyze

    # stats
    min_cells: int = 10
    n_permutations: int = 10_000
    perm_seed: Optional[int] = None
    fdr_alpha: float = 0.05

    # gene-level (optional)
    gene_min_nonzero: int = 3
    gene_fdr_alpha: float = 0.05
    run_gene_permutations: bool = True
    gene_n_permutations: int = 10_000
    gene_perm_seed: Optional[int] = None

    # UX
    verbose: bool = True


# -----------------------------
# Logging / IO helpers
# -----------------------------
def _log(params: VelocityStage2Params, msg: str) -> None:
    if params.verbose:
        print(msg, flush=True)


def ensure_dir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p

def load_summary_table(summary_table_path: Path) -> pd.DataFrame:
    labels_full = pd.read_csv(summary_table_path, dtype={"Var1": str})  # enforce Var1 as string
    if "tissue" in labels_full.columns:
        labels_full["tissue"] = labels_full["tissue"].astype(str)  

    if "Var1" not in labels_full.columns:
        raise ValueError(
            "summary_table_path must contain a 'Var1' column with unique cell IDs. "
            f"Missing 'Var1' in: {summary_table_path}"
        )

    labels_full["Var1"] = labels_full["Var1"].astype(str)
    labels_full = labels_full.set_index("Var1", drop=True)

    return labels_full

def load_cells_df(pc_origin_path: Path, params: VelocityStage2Params) -> pd.DataFrame:
    df = pd.read_csv(pc_origin_path, index_col=0, dtype={0: str})
    df.index = df.index.astype(str)

    # basic checks
    for c in [params.pc1_col, params.pc2_col]:
        if c not in df.columns:
            raise ValueError(f"pc_origin CSV must contain '{c}' column.")

    # tissue (optional)
    if params.parse_tissue_from_index and params.tissue_col not in df.columns:
        df[params.tissue_col] = df.index.to_series().astype(str).str.split(".").str[0].astype(str)

    return df


def load_pickle_dict(gene_dict_pkl: Path) -> dict:
    import pickle

    with open(gene_dict_pkl, "rb") as f:
        return pickle.load(f)


def _safe_unique_cell_types(x: pd.Series) -> List[str]:
    s = x.replace(["nan", "NA", "None", ""], np.nan).dropna().astype(str)
    return sorted(s.unique().tolist())


def bh_fdr(pvals: pd.Series, alpha: float = 0.05) -> pd.Series:
    """
    Benjamini–Hochberg adjusted p-values (q-values).
    Returns a Series aligned to pvals.index.
    """
    out = pd.Series(np.nan, index=pvals.index, dtype=float)

    pv = pvals.dropna().astype(float)
    if pv.empty:
        return out

    if multipletests is None:
        raise ImportError(
            "statsmodels is required for BH-FDR. Install with: pip install statsmodels"
        )

    _, q, _, _ = multipletests(pv.values, alpha=alpha, method="fdr_bh")
    out.loc[pv.index] = q
    return out


# -----------------------------
# Proximity computation 
# -----------------------------
def compute_proximity_continuous_and_binary(
    summary_table_path: Path,
    distance_dir: Path,
    distance_threshold: float,
    tissue_col: str = "tissue",
    cell_type_col: str = "cell_type",
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Computes:
      1) prox_df_continuous: per cell, min distance to each cell type (NaN if none)
      2) prox_df_binary: per cell, 1 if min distance <= threshold else 0

    This follows exactly the logic you provided (with safety casts and directory checks).
    """
    if verbose:
        print(f"[PROX] Loading summary table: {summary_table_path}", flush=True)
        print(f"[PROX] Loading distance matrices from: {distance_dir}", flush=True)
        print(f"[PROX] distance_threshold = {distance_threshold}", flush=True)

    full = load_summary_table(summary_table_path)
    full.index = full.index.astype(str)

    if tissue_col not in full.columns:
        raise ValueError(f"summary table must include '{tissue_col}' column.")
    if cell_type_col not in full.columns:
        raise ValueError(f"summary table must include '{cell_type_col}' column.")

    full[tissue_col] = full[tissue_col].astype(str)

    cell_types = pd.Series(full[cell_type_col].unique()).dropna().astype(str).tolist()
    tissues = sorted(pd.Series(full[tissue_col].unique()).dropna().astype(str).tolist())

    # continuous
    prox_cont = pd.DataFrame(np.nan, index=full.index.astype(str), columns=cell_types)

    # binary
    prox_bin = pd.DataFrame(0, index=full.index.astype(str), columns=cell_types, dtype=int)

    for tissue in tissues:
        dm_path = distance_dir / f"distance_matrix_{tissue}.csv"
        if not dm_path.exists():
            if verbose:
                print(f"[PROX] WARNING: missing {dm_path} (skip tissue={tissue})", flush=True)
            continue

        if verbose:
            print(f"[PROX] tissue={tissue} | reading {dm_path.name}", flush=True)

        dm = pd.read_csv(dm_path, index_col=0, dtype={0: str})
        #dm.index = dm.index.astype(str)
        #dm.columns = dm.columns.astype(str)
        dm.index = dm.index.map(str)
        dm.columns = dm.columns.map(str)

        tissue_cells = full.loc[full[tissue_col].astype(str) == str(tissue)].index.astype(str).tolist()
        tissue_cells = list(set(tissue_cells).intersection(dm.index))
        if len(tissue_cells) == 0:
            continue

        for ctype in cell_types:
            cells_X = (
                full.loc[
                    (full[tissue_col].astype(str) == str(tissue))
                    & (full[cell_type_col].astype(str) == str(ctype))
                ]
                .index.astype(str)
                .tolist()
            )
            cells_X = list(set(cells_X).intersection(dm.columns))
            if len(cells_X) == 0:
                # continuous stays NaN; binary stays 0
                continue

            sub = dm.loc[tissue_cells, cells_X]

            # exclude self distance
            common = set(sub.index).intersection(sub.columns)
            if common:
                for cid in common:
                    sub.loc[cid, cid] = np.nan

            min_dist = sub.min(axis=1)  # per-cell min distance to that type

            # continuous assignment
            prox_cont.loc[min_dist.index, str(ctype)] = min_dist.values

            # binary assignment
            is_close = (min_dist <= float(distance_threshold)).fillna(False).astype(int)
            prox_bin.loc[is_close.index, str(ctype)] = is_close.values

    # add metadata column at left
    prox_cont = prox_cont.join(full[cell_type_col].rename(cell_type_col), how="left")
    prox_bin = prox_bin.join(full[cell_type_col].rename(cell_type_col), how="left")

    cols = [cell_type_col] + [c for c in prox_cont.columns if c != cell_type_col]
    prox_cont = prox_cont[cols]
    prox_bin = prox_bin[cols]

    return prox_cont, prox_bin



def create_min_dist_vector_from_coordinates(
    tissue_list,
    cell_type_1: str,
    cell_type_2: str,
    labels: pd.DataFrame,
    maximum_distance: float = np.inf,
    tissue_col: str = "tissue",
    cell_type_col: str = "cell_type",
    x_col: str = "X_space",
    y_col: str = "Y_space",
    z_col: str = "Z_space",
):
    """
    Compute the minimum Euclidean distance from each cell in cell_type_1 to the nearest cell in cell_type_2,
    for each tissue in the list, using X_space and Y_space coordinates.
    Distances above `maximum_distance` are replaced with NaN.
    """
    min_dist_vector = pd.Series(dtype=float)

    for tissue in tissue_list:
        print(f"Processing tissue: {tissue}")

        # Filter cells for this tissue
        df_tissue = labels.loc[labels[tissue_col].astype(str) == str(tissue)]

        # Subset by cell type
        df_ct1 = df_tissue[df_tissue[cell_type_col] == cell_type_1].copy()
        df_ct2 = df_tissue[df_tissue[cell_type_col] == cell_type_2].copy()

        if df_ct1.empty or df_ct2.empty:
            print(f"Skipping tissue {tissue}: one or both cell types are missing.")
            continue

        # Get coordinates
        if z_col in df_tissue.columns:
            coords_ct1 = df_ct1[[x_col, y_col, z_col]].values
            coords_ct2 = df_ct2[[x_col, y_col, z_col]].values
            print("Using 3D coordinates (X, Y, Z).")
        else:
            print("Warning: 'Z_space' column not found, falling back to 2D (X, Y).")
            coords_ct1 = df_ct1[[x_col, y_col]].values
            coords_ct2 = df_ct2[[x_col, y_col]].values

        # Build KDTree for fast nearest neighbor search
        tree_ct2 = cKDTree(coords_ct2)

        # For each cell in ct1, find the closest cell in ct2
        distances, _ = tree_ct2.query(coords_ct1, k=1)

        # Apply threshold: set distances > maximum_distance to NaN
        distances = np.where(distances > maximum_distance, np.nan, distances)

        # Assign distances using ct1 cell IDs, prefixed with tissue number
        min_dist_vector = pd.concat([min_dist_vector, pd.Series(distances, index=df_ct1.index)])

        #print(
        #    f"The initial number of {cell_type_1} cells is {(len(df_ct1))} and the number of {cell_type_1} cells with neighbors is {len(distances)} and within the maximum distance is {np.sum(~np.isnan(distances))}.")
        #print(
        #    f"{cell_type_1} → {cell_type_2} | {len(distances)} cells | Min: {np.nanmin(distances):.2f}, Max: {np.nanmax(distances):.2f}")

    return min_dist_vector


def compute_min_distance_from_x_y(
    summary_table_path: Path,
    tissue_col: str = "tissue",
    cell_type_col: str = "cell_type",
    maximum_distance_threshold: float = np.inf,   # ~FOV size
    distance_threshold: float = 1.0,  
    verbose: bool = True,
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    labels = load_summary_table(summary_table_path)
    labels.index = labels.index.astype(str)
    labels[tissue_col] = labels[tissue_col].astype(str)
    all_cell_types = labels[cell_type_col].dropna().astype(str).unique()
    tissue_list = labels[tissue_col].dropna().astype(str).unique()
    

    # Output DF: rows=cells, cols=neighbor cell types (ct2)
    min_dist_df = pd.DataFrame(
        np.nan,
        index=labels.index.astype(str),
        columns=all_cell_types
    )

    # Optional: keep cell_type as metadata column (like you did elsewhere)
    #min_dist_df.insert(0, "cell_type", labels["cell_type"].astype(object))

    # Fill DF using the SAME nested loops logic you already use
    for ct1 in all_cell_types:
        for ct2 in all_cell_types:
            if ct1 == ct2:
                continue

            min_dist = create_min_dist_vector_from_coordinates(
                tissue_list=tissue_list,
                cell_type_1=ct1,
                cell_type_2=ct2,
                labels=labels,
                maximum_distance=maximum_distance_threshold
            )

            # Write into the single DF: only ct1 rows get values for this ct2 column
            min_dist_df.loc[min_dist.index.astype(str), ct2] = min_dist.values


    prox_cont = min_dist_df.copy()
    prox_bin = (min_dist_df <= float(distance_threshold)).astype(int)

    prox_cont = prox_cont.join(labels[cell_type_col].rename(cell_type_col), how="left")
    prox_bin  = prox_bin.join(labels[cell_type_col].rename(cell_type_col), how="left")
    cols = [cell_type_col] + [c for c in prox_cont.columns if c != cell_type_col]
    prox_cont = prox_cont[cols]
    prox_bin  = prox_bin[cols]


    return prox_cont, prox_bin

# -----------------------------
# Math helpers (fast correlations)
# -----------------------------
def _pearsonr_fast(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """
    Pearson r + asymptotic p-value (two-sided) via t approximation.
    We use it only as a rough p-value (main significance comes from permutations + FDR).
    """
    x = x.astype(float)
    y = y.astype(float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    n = x.size
    if n < 3:
        return np.nan, np.nan

    x0 = x - x.mean()
    y0 = y - y.mean()
    denom = np.sqrt((x0**2).sum()) * np.sqrt((y0**2).sum())
    if denom == 0:
        return np.nan, np.nan
    r = float((x0 * y0).sum() / denom)

    # t-stat
    if abs(r) >= 1.0:
        return r, 0.0
    t = r * np.sqrt((n - 2) / (1 - r * r))
    # two-sided p using normal approx on t for simplicity (OK for large n)
    p = 2 * (1 - norm.cdf(abs(t)))
    return r, float(p)


def _pointbiserial_fast(x01: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """
    Point-biserial correlation for binary x (0/1) and continuous y.
    """
    x01 = x01.astype(float)
    y = y.astype(float)
    m = np.isfinite(x01) & np.isfinite(y)
    x01 = x01[m]
    y = y[m]
    n = x01.size
    if n < 3:
        return np.nan, np.nan

    g1 = y[x01 > 0.5]
    g0 = y[x01 <= 0.5]
    if g1.size == 0 or g0.size == 0:
        return np.nan, np.nan

    m1 = g1.mean()
    m0 = g0.mean()
    s = y.std(ddof=1)
    if s == 0:
        return np.nan, np.nan

    p = g1.size / n
    q = 1.0 - p
    r = (m1 - m0) / s * np.sqrt(p * q)

    # rough p-value via normal approx
    t = r * np.sqrt((n - 2) / (1 - r * r)) if abs(r) < 1 else np.inf
    pval = 2 * (1 - norm.cdf(abs(t)))
    return float(r), float(pval)


def _one_tailed_perm_p_from_z(observed: float, perm_mean: float, perm_std: float) -> float:
    if not np.isfinite(observed) or not np.isfinite(perm_mean) or not np.isfinite(perm_std) or perm_std == 0:
        return np.nan
    z = (observed - perm_mean) / perm_std
    if observed < 0:
        return float(norm.cdf(z))
    if observed > 0:
        return float(1 - norm.cdf(z))
    return np.nan


# -----------------------------
# Stage2 core computations
# -----------------------------
def compute_projected_end_state_pcs(
    st_norm_filtered: pd.DataFrame,
    n_components: int = 3,
) -> pd.DataFrame:
    """
    Compute PCs for s(t) end state (from Stage2 outputs inside the pkl dict),
    using standardized st and returning PC1..PCn for all cells.
    """
    # match your original scaling: St * mean(sum per cell)
    norm_st = st_norm_filtered * st_norm_filtered.sum(axis=1).mean()
    pcs = PCA(n_components=int(n_components)).fit_transform(StandardScaler().fit_transform(norm_st))
    cols = [f"PC{i+1}_end" for i in range(pcs.shape[1])]
    return pd.DataFrame(pcs, index=st_norm_filtered.index.astype(str), columns=cols)


def build_cells_info(
    cells_df: pd.DataFrame,
    prox_df: pd.DataFrame,
    st_end_pcs: pd.DataFrame,
    params: VelocityStage2Params,
) -> pd.DataFrame:
    """
    Merge origin PCs + projected end PCs + proximity columns.
    Add delta vectors and magnitudes and phase angles.
    """
    df = cells_df.copy()
    df.index = df.index.astype(str)

    # attach end PCs
    df = df.join(st_end_pcs, how="left")

    # ensure end PC columns exist for PC1/PC2
    if "PC1_end" not in df.columns or "PC2_end" not in df.columns:
        raise ValueError("Missing PC1_end/PC2_end in end-state PCs table.")

    # delta vector (PC space)
    df["delta_x"] = df["PC1_end"] - df[params.pc1_col]
    df["delta_y"] = df["PC2_end"] - df[params.pc2_col]

    # magnitudes
    df["magnitude_s0"] = np.sqrt(df[params.pc1_col] ** 2 + df[params.pc2_col] ** 2)
    df["magnitude_delta"] = np.sqrt(df["delta_x"] ** 2 + df["delta_y"] ** 2)

    # phases (degrees)
    df["phase_deg_delta"] = (np.degrees(np.arctan2(df["delta_y"], df["delta_x"])) % 360).astype(float)
    df["phase_deg_s0"] = (np.degrees(np.arctan2(df[params.pc2_col], df[params.pc1_col])) % 360).astype(float)

    # merge proximity (exclude duplicate cell_type column handling)
    if params.cell_type_col not in prox_df.columns:
        raise ValueError(f"prox_df must include '{params.cell_type_col}' column.")
    prox_cols = [c for c in prox_df.columns if c != params.cell_type_col]
    df = df.join(prox_df[prox_cols], how="left")

    return df


def _get_neighbor_columns(prox_df: pd.DataFrame, cell_type_col: str) -> List[str]:
    return [c for c in prox_df.columns if c != cell_type_col]


def compute_pairwise_correlations(
    cells_info: pd.DataFrame,
    neighbor_cols: List[str],
    params: VelocityStage2Params,
    y_feature: YFeature,
) -> pd.DataFrame:
    """
    For each primary cell type and each neighbor cell type:
      correlate proximity(primary->neighbor) with selected y_feature among primary cells only.
    """
    if params.cell_type_col not in cells_info.columns:
        raise ValueError(f"cells_info missing '{params.cell_type_col}'.")

    primary_types = _safe_unique_cell_types(cells_info[params.cell_type_col])

    rows = []
    for primary_ct in primary_types:
        df_primary = cells_info.loc[cells_info[params.cell_type_col].astype(str) == str(primary_ct)]

        for neighbor_ct in neighbor_cols:
            if neighbor_ct == params.cell_type_col:
                continue

            # use only cells with defined proximity value for this neighbor
            df_use = df_primary.loc[df_primary[neighbor_ct].dropna().index]
            n_cells = int(df_use.shape[0])
            if n_cells < params.min_cells:
                continue

            x = df_use[neighbor_ct].to_numpy()
            y = df_use[y_feature].to_numpy()

            try:
                if params.proximity_mode == "continuous":
                    r, p = _pearsonr_fast(x, y)
                else:
                    r, p = _pointbiserial_fast(x, y)
            except Exception:
                r, p = np.nan, np.nan

            rows.append(
                {
                    "primary_cell_type": str(primary_ct),
                    "neighbor_cell_type": str(neighbor_ct),
                    "correlation": r,
                    "p_value": p,
                    "n_cells": n_cells,
                    "mean_proximity": float(np.nanmean(x)) if np.isfinite(x).any() else np.nan,
                    "mean_y": float(np.nanmean(y)) if np.isfinite(y).any() else np.nan,
                    "y_feature": y_feature,
                }
            )

    res = pd.DataFrame(rows)
    if res.empty:
        return res

    # BH-FDR on raw p-values
    res["p_value_fdr"] = bh_fdr(res["p_value"], alpha=params.fdr_alpha)
    return res.sort_values(["y_feature", "p_value"], ascending=[True, True])


def _permute_correlations(
    df_primary: pd.DataFrame,
    neighbor_ct: str,
    y_feature: YFeature,
    params: VelocityStage2Params,
    rng: np.random.Generator,
) -> np.ndarray:
    """
    Permute x (proximity) among primary cells, recompute correlation, return vector of permuted r.
    """
    x = df_primary[neighbor_ct].to_numpy().astype(float)
    y = df_primary[y_feature].to_numpy().astype(float)

    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < params.min_cells:
        return np.array([], dtype=float)

    out = np.full(params.n_permutations, np.nan, dtype=float)
    for i in range(params.n_permutations):
        x_perm = x.copy()
        rng.shuffle(x_perm)
        if params.proximity_mode == "continuous":
            r, _ = _pearsonr_fast(x_perm, y)
        else:
            r, _ = _pointbiserial_fast(x_perm, y)
        out[i] = r
    return out


def add_pairwise_permutation_pvalues(
    cells_info: pd.DataFrame,
    pairwise_df: pd.DataFrame,
    params: VelocityStage2Params,
) -> pd.DataFrame:
    """
    Add one-tailed permutation p-values + BH-FDR for each (primary, neighbor, y_feature).
    """
    if pairwise_df.empty:
        return pairwise_df

    rng = np.random.default_rng(params.perm_seed)
    pairwise_df = pairwise_df.copy()
    pairwise_df["p_value_perm"] = np.nan

    # group by primary to avoid repeated filtering
    grouped = pairwise_df.groupby(["y_feature", "primary_cell_type"], sort=False)
    for (y_feature, primary_ct), g in grouped:
        df_primary = cells_info.loc[cells_info[params.cell_type_col].astype(str) == str(primary_ct)]
        for idx, row in g.iterrows():
            neighbor_ct = row["neighbor_cell_type"]

            df_use = df_primary.loc[df_primary[neighbor_ct].dropna().index]
            if df_use.shape[0] < params.min_cells:
                continue

            perm_r = _permute_correlations(df_use, neighbor_ct, y_feature, params, rng)
            if perm_r.size == 0 or not np.isfinite(perm_r).any():
                continue

            obs = float(row["correlation"])
            perm_mean = float(np.nanmean(perm_r))
            perm_std = float(np.nanstd(perm_r, ddof=1))
            p_perm = _one_tailed_perm_p_from_z(obs, perm_mean, perm_std)
            pairwise_df.at[idx, "p_value_perm"] = p_perm

    pairwise_df["p_value_perm_fdr"] = bh_fdr(pairwise_df["p_value_perm"], alpha=params.fdr_alpha)
    return pairwise_df


def run_pairwise_by_tissue(
    cells_info: pd.DataFrame,
    neighbor_cols: List[str],
    params: VelocityStage2Params,
) -> pd.DataFrame:
    """
    Tissue-stratified pairwise correlations. Same outputs + tissue column.
    """
    if params.tissue_col not in cells_info.columns:
        raise ValueError(f"cells_info missing '{params.tissue_col}'.")

    out_rows = []
    tissues = sorted(cells_info[params.tissue_col].dropna().astype(str).unique().tolist())
    for tissue in tissues:
        _log(params, f"[TISSUE] pairwise tissue={tissue}")
        df_t = cells_info.loc[cells_info[params.tissue_col].astype(str) == str(tissue)]
        for y_feature in ["magnitude_delta", "magnitude_s0", "phase_deg_delta", "phase_deg_s0"]:
            pw = compute_pairwise_correlations(df_t, neighbor_cols, params, y_feature=y_feature)  # raw + fdr
            pw = add_pairwise_permutation_pvalues(df_t, pw, params)
            if not pw.empty:
                pw["tissue"] = str(tissue)
                out_rows.append(pw)

    if not out_rows:
        return pd.DataFrame()

    return pd.concat(out_rows, ignore_index=True)


# -----------------------------
# Gene-level (optional)
# -----------------------------
def correlate_gene_velocity_with_proximity(
    v_abs: pd.DataFrame,
    x: pd.Series,
    params: VelocityStage2Params,
) -> pd.DataFrame:
    """
    For each gene: corr(|v_gene|, proximity_x) across the selected cells.
    Uses Pearson (continuous x) or point-biserial (binary x).
    """
    x_arr = x.to_numpy().astype(float)
    rows = []
    for gene in v_abs.columns:
        y_arr = v_abs[gene].to_numpy().astype(float)

        # skip mostly-zero genes
        if np.sum(y_arr != 0) < params.gene_min_nonzero:
            continue

        if params.proximity_mode == "continuous":
            r, p = _pearsonr_fast(x_arr, y_arr)
        else:
            r, p = _pointbiserial_fast(x_arr, y_arr)

        rows.append({"gene": gene, "correlation": r, "p_value": p})

    df = pd.DataFrame(rows)
    if df.empty:
        return df
    df["p_value_fdr"] = bh_fdr(df["p_value"], alpha=params.gene_fdr_alpha)
    return df.sort_values("p_value", ascending=True)


def add_gene_permutation_pvalues(
    gene_df: pd.DataFrame,
    v_abs: pd.DataFrame,
    x: pd.Series,
    params: VelocityStage2Params,
) -> pd.DataFrame:
    """
    Optional: permutation p-values per gene (VERY slow for many genes).
    One-tailed, using z-score vs permuted correlations.
    """
    if gene_df.empty:
        return gene_df

    rng = np.random.default_rng(params.gene_perm_seed)
    x_arr = x.to_numpy().astype(float)
    gene_df = gene_df.copy()
    gene_df["p_value_perm"] = np.nan

    for i, row in gene_df.iterrows():
        gene = row["gene"]
        y = v_abs[gene].to_numpy().astype(float)

        m = np.isfinite(x_arr) & np.isfinite(y)
        x0 = x_arr[m]
        y0 = y[m]
        if x0.size < 3:
            continue

        perm_rs = np.full(params.gene_n_permutations, np.nan, dtype=float)
        for k in range(params.gene_n_permutations):
            x_perm = x0.copy()
            rng.shuffle(x_perm)
            if params.proximity_mode == "continuous":
                r, _ = _pearsonr_fast(x_perm, y0)
            else:
                r, _ = _pointbiserial_fast(x_perm, y0)
            perm_rs[k] = r

        obs = float(row["correlation"])
        perm_mean = float(np.nanmean(perm_rs))
        perm_std = float(np.nanstd(perm_rs, ddof=1))
        gene_df.at[i, "p_value_perm"] = _one_tailed_perm_p_from_z(obs, perm_mean, perm_std)

    gene_df["p_value_perm_fdr"] = bh_fdr(gene_df["p_value_perm"], alpha=params.gene_fdr_alpha)
    return gene_df


def run_gene_level_for_significant_pairs(
    cells_info: pd.DataFrame,
    gene_dict: dict,
    pairwise_feature_df: pd.DataFrame,
    params: VelocityStage2Params,
    select_pval_col: str = "p_value_perm_fdr",
    alpha: float = 0.05,
) -> Dict[str, pd.DataFrame]:
    """
    For each significant primary→neighbor pair (based on pairwise feature table):
      correlate abs velocity per gene with proximity to that neighbor within primary cells.

    Returns dict: "<primary>__<neighbor>" -> gene_df
    """
    if pairwise_feature_df.empty:
        return {}

    if "filtered_v" not in gene_dict:
        raise ValueError("gene_dict missing key 'filtered_v'.")
    v = gene_dict["filtered_v"]
    if not isinstance(v, pd.DataFrame):
        raise ValueError("gene_dict['filtered_v'] must be a pandas DataFrame.")
    v.index = v.index.astype(str)

    # significant pairs
    if select_pval_col not in pairwise_feature_df.columns:
        raise ValueError(f"pairwise_feature_df missing '{select_pval_col}' column.")

    sig = pairwise_feature_df.loc[pairwise_feature_df[select_pval_col] <= float(alpha)].copy()
    if sig.empty:
        return {}

    out: Dict[str, pd.DataFrame] = {}
    for _, row in sig.iterrows():
        primary_ct = str(row["primary_cell_type"])
        neighbor_ct = str(row["neighbor_cell_type"])
        key = f"{primary_ct}__{neighbor_ct}"

        df_primary = cells_info.loc[cells_info[params.cell_type_col].astype(str) == primary_ct]
        df_primary = df_primary.loc[df_primary[neighbor_ct].dropna().index]
        if df_primary.shape[0] < params.min_cells:
            continue

        x = df_primary[neighbor_ct]
        idx = df_primary.index.astype(str)

        v_sub = v.loc[v.index.intersection(idx)]
        # align x to v_sub
        x = x.loc[v_sub.index]
        v_abs = v_sub.abs()

        _log(params, f"[GENE] {primary_ct}→{neighbor_ct} | cells={v_abs.shape[0]} | genes={v_abs.shape[1]}")
        gene_df = correlate_gene_velocity_with_proximity(v_abs, x, params)

        if params.run_gene_permutations and not gene_df.empty:
            gene_df = add_gene_permutation_pvalues(gene_df, v_abs, x, params)

        out[key] = gene_df

    return out


# -----------------------------
# Orchestrator + saving
# -----------------------------
def save_stage2_outputs(
    out_dir: Path,
    params: VelocityStage2Params,
    cells_info: pd.DataFrame,
    prox_cont: pd.DataFrame,
    prox_bin: pd.DataFrame,
    pairwise_all: Dict[str, pd.DataFrame],
    pairwise_tissue_df: Optional[pd.DataFrame],
    gene_level: Dict[str, pd.DataFrame],
) -> None:
    out_dir = ensure_dir(out_dir)
    inputs_dir = ensure_dir(out_dir / "stage2_inputs")
    outputs_dir = ensure_dir(out_dir / "stage2_outputs")

    # save proximity tables (requested)
    prox_cont_path = inputs_dir / "proximity_to_cell_types_continuous.csv"
    prox_bin_path = inputs_dir / f"proximity_to_cell_types_binary_thr{params.distance_threshold:g}.csv"
    prox_cont.to_csv(prox_cont_path)
    prox_bin.to_csv(prox_bin_path)
    _log(params, f"[SAVE] {prox_cont_path}")
    _log(params, f"[SAVE] {prox_bin_path}")

    # cells_info
    cells_info_path = outputs_dir / f"cells_info_{params.proximity_mode}.csv"
    cells_info.to_csv(cells_info_path)
    _log(params, f"[SAVE] {cells_info_path}")

    # pairwise feature tables
    for y_feature, df in pairwise_all.items():
        if df is None or df.empty:
            continue
        p = outputs_dir / f"pairwise_{y_feature}_{params.proximity_mode}.csv"
        df.to_csv(p, index=False)
        _log(params, f"[SAVE] {p}")

    # tissue table (optional)
    if pairwise_tissue_df is not None and not pairwise_tissue_df.empty:
        p = outputs_dir / f"pairwise_by_tissue_{params.proximity_mode}.csv"
        pairwise_tissue_df.to_csv(p, index=False)
        _log(params, f"[SAVE] {p}")

    # gene-level dict -> per-pair CSVs
    if gene_level:
        gene_dir = ensure_dir(outputs_dir / "gene_level")
        for key, gdf in gene_level.items():
            if gdf is None or gdf.empty:
                continue
            p = gene_dir / f"gene_velocity_vs_proximity__{key}__{params.proximity_mode}.csv"
            gdf.to_csv(p, index=False)
            _log(params, f"[SAVE] {p}")


def run_velocity_stage2_pipeline(
    pc_origin_path: Path,
    summary_table_path: Path,
    distance_dir: Optional[Path],
    gene_dict_pkl: Path,
    out_dir: Path,
    params: VelocityStage2Params,
    compute_pair_permutations: bool = True,
    compute_tissue: bool = True,
    compute_gene_level: bool = True,
    gene_from_feature: Literal["magnitude_delta", "magnitude_s0", "phase_deg_delta", "phase_deg_s0"] = "magnitude_delta",
    gene_pval_col: str = "p_value_perm_fdr",
    gene_alpha: float = 0.05,
) -> None:
    """
    Stage 2 pipeline (updated):
      1) Load origin PCs (pc_origin_path)
      2) Compute prox_df_continuous and prox_df_binary from distance matrices (distance_dir + summary_table)
         and save them as CSV under out_dir/stage2_inputs
      3) Pick one of them (params.proximity_mode) and run:
         - pairwise proximity↔(magnitude_delta, magnitude_s0, phase_deg_delta, phase_deg_s0)
           including permutations + BH-FDR
         - optional tissue-stratified analysis
         - optional gene-level velocity↔proximity for significant primary→neighbor pairs

    NOTE: "primary" = the cell type whose cells are being correlated;
          "neighbor" = the proximity column / neighboring cell type.
    """
    _log(params, f"[LOAD] pc_origin_path={pc_origin_path}")
    cells_df = load_cells_df(pc_origin_path, params)

    _log(params, f"[LOAD] gene_dict_pkl={gene_dict_pkl}")
    gene_dict = load_pickle_dict(gene_dict_pkl)

    if distance_dir is None:
        _log(params, "[PROX] distance_dir not provided -> computing proximity from X/Y(/Z) coordinates")
        prox_cont, prox_bin = compute_min_distance_from_x_y(
            summary_table_path=summary_table_path,
            tissue_col=params.tissue_col,
            cell_type_col=params.cell_type_col,
            distance_threshold=params.distance_threshold,
            maximum_distance_threshold=params.maximum_distance_threshold,
            verbose=params.verbose,
        )
        _log(params, f"[PROX] DONE building prox_cont/prox_bin from coordinates | shape_cont={prox_cont.shape} | shape_bin={prox_bin.shape}")
    else:
        _log(params, "[PROX] distance_dir provided -> computing proximity from distance matrices")
        prox_cont, prox_bin = compute_proximity_continuous_and_binary(
            summary_table_path=summary_table_path,
            distance_dir=distance_dir,
            distance_threshold=params.distance_threshold,
            tissue_col=params.tissue_col,
            cell_type_col=params.cell_type_col,
            verbose=params.verbose,
        )
        _log(params, f"[PROX] DONE building prox_cont/prox_bin from distance matrices | shape_cont={prox_cont.shape} | shape_bin={prox_bin.shape}")


    # choose prox df to analyze
    prox_df = prox_cont if params.proximity_mode == "continuous" else prox_bin
    neighbor_cols = _get_neighbor_columns(prox_df, params.cell_type_col)

    # ---- compute end-state PCs from gene_dict at t ----
    # your pkl is "all_gene_level_dfs_per_t" OR sometimes stored directly;
    # support both patterns robustly.
    if "all_gene_level_dfs_per_t" in gene_dict:
        per_t = gene_dict["all_gene_level_dfs_per_t"]
    else:
        per_t = gene_dict

    if params.t not in per_t:
        # sometimes keys are strings
        if str(params.t) in per_t:
            t_key = str(params.t)
        else:
            raise ValueError(f"t={params.t} not found in gene_dict keys. Available: {list(per_t.keys())[:10]} ...")
    else:
        t_key = params.t

    t_block = per_t[t_key]
    if "st_norm_filtered" not in t_block:
        raise ValueError("gene_dict for selected t missing 'st_norm_filtered'.")
    if "filtered_v" not in t_block:
        raise ValueError("gene_dict for selected t missing 'filtered_v'.")

    st_norm_filtered = t_block["st_norm_filtered"]
    st_norm_filtered.index = st_norm_filtered.index.astype(str)

    _log(params, f"[PCA] computing end-state PCs from st_norm_filtered (cells={st_norm_filtered.shape[0]}, genes={st_norm_filtered.shape[1]})")
    st_end_pcs = compute_projected_end_state_pcs(st_norm_filtered, n_components=params.n_pcs)

    # ---- build cells_info ----
    _log(params, f"[MERGE] building cells_info ({params.proximity_mode})")
    cells_info = build_cells_info(cells_df=cells_df, prox_df=prox_df, st_end_pcs=st_end_pcs, params=params)

    # ---- pairwise analysis (4 y-features) ----
    pairwise_all: Dict[str, pd.DataFrame] = {}
    for y_feature in ["magnitude_delta", "magnitude_s0", "phase_deg_delta", "phase_deg_s0"]:
        _log(params, f"[PAIRWISE] y_feature={y_feature} | mode={params.proximity_mode}")
        pw = compute_pairwise_correlations(cells_info, neighbor_cols, params, y_feature=y_feature)
        if compute_pair_permutations and not pw.empty:
            _log(params, f"[PERM] pairwise permutations for y_feature={y_feature} (n_perm={params.n_permutations})")
            pw = add_pairwise_permutation_pvalues(cells_info, pw, params)
        pairwise_all[y_feature] = pw

    # ---- tissue stratified ----
    pairwise_tissue_df: Optional[pd.DataFrame] = None
    if compute_tissue:
        _log(params, f"[TISSUE] running tissue-stratified analysis")
        pairwise_tissue_df = run_pairwise_by_tissue(cells_info, neighbor_cols, params)

    # ---- gene-level ----
    gene_level: Dict[str, pd.DataFrame] = {}
    if compute_gene_level:
        _log(params, f"[GENE] selecting significant pairs from '{gene_from_feature}' using '{gene_pval_col}' <= {gene_alpha}")
        feature_df = pairwise_all.get(gene_from_feature, pd.DataFrame())
        if feature_df is not None and not feature_df.empty:
            gene_level = run_gene_level_for_significant_pairs(
                cells_info=cells_info,
                gene_dict=t_block,  # IMPORTANT: use the selected t-block (contains filtered_v)
                pairwise_feature_df=feature_df,
                params=params,
                select_pval_col=gene_pval_col,
                alpha=gene_alpha,
            )
        else:
            _log(params, "[GENE] skipped (no pairwise feature table / empty results).")

    # ---- save ----
    save_stage2_outputs(
        out_dir=out_dir,
        params=params,
        cells_info=cells_info,
        prox_cont=prox_cont,
        prox_bin=prox_bin,
        pairwise_all=pairwise_all,
        pairwise_tissue_df=pairwise_tissue_df,
        gene_level=gene_level,
    )

    _log(params, f"[DONE] Stage 2 completed. Outputs under: {out_dir}")


# Convenience wrapper (kept for CLI call style)
def run_velocity_stage2(
    pc_origin_path: Path,
    summary_table_path: Path,
    distance_dir: Optional[Path],
    gene_dict_pkl: Path,
    out_dir: Path,
    t: float,
    proximity_mode: ProximityMode = "continuous",
    distance_threshold: float = 1.0,
    maximum_distance_threshold: float = np.inf,
    min_cells: int = 10,
    n_perm: int = 10_000,
    rng_seed: Optional[int] = None,
    fdr_alpha: float = 0.05,
    compute_pair_permutations: bool = True,
    compute_tissue: bool = True,
    compute_gene_level: bool = True,
    gene_from_feature: Literal["magnitude_delta", "magnitude_s0", "phase_deg_delta", "phase_deg_s0"] = "magnitude_delta",
    gene_pval_col: str = "p_value_perm_fdr",
    gene_alpha: float = 0.05,
    run_gene_permutations: bool = True,
    gene_n_perm: int = 10_000,
    gene_rng_seed: Optional[int] = None,
    verbose: bool = True,
) -> None:
    params = VelocityStage2Params(
        t=float(t),
        proximity_mode=proximity_mode,
        distance_threshold=float(distance_threshold),
        maximum_distance_threshold=float(maximum_distance_threshold),
        min_cells=int(min_cells),
        n_permutations=int(n_perm),
        perm_seed=rng_seed,
        fdr_alpha=float(fdr_alpha),
        run_gene_permutations=bool(run_gene_permutations),
        gene_n_permutations=int(gene_n_perm),
        gene_perm_seed=gene_rng_seed,
        verbose=bool(verbose),
    )

    run_velocity_stage2_pipeline(
        pc_origin_path=pc_origin_path,
        summary_table_path=summary_table_path,
        distance_dir=distance_dir,
        gene_dict_pkl=gene_dict_pkl,
        out_dir=out_dir,
        params=params,
        compute_pair_permutations=compute_pair_permutations,
        compute_tissue=compute_tissue,
        compute_gene_level=compute_gene_level,
        gene_from_feature=gene_from_feature,
        gene_pval_col=gene_pval_col,
        gene_alpha=gene_alpha,
    )
