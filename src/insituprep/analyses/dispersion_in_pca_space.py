"""
dispersion_in_pca_space.py

Implements the "Dispersion in PCA space" analysis.

GOAL (conceptual):
For each tissue and each ordered pair of cell types (cell_type_1 -> cell_type_2),
we embed the primary cell type (cell_type_1) into a 2D PCA space built from its gene-expression profiles.
We then label each primary cell as Proximal / Non-proximal to the neighbor cell type (cell_type_2)
based on the physical distance matrix and a user-defined threshold.

We quantify how "separated" proximal vs non-proximal primary cells are in PCA space using:
  - Single-centroid statistic: sum of distances of proximal points to the overall centroid
  - Two-centroid statistic: distance between proximal centroid and non-proximal centroid
Then we run permutations to build null distributions, convert each statistic to a z-score,
take max(z_single, z_two), and compute a Gaussian-approx p-value for the max-z statistic.

This module also supports:
  1) Shuffled-proximity sanity-check across all pairs
  2) Marker-removal sensitivity test on significant pairs
  3) t-SNE sensitivity test on significant pairs
  4) t-SNE sensitivity with shuffled proximity labels on all pairs
  5) A multi-panel PCA plot for significant pairs (saved to output folder)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple, Mapping

import json
import warnings

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.metrics.pairwise import euclidean_distances
from sklearn.manifold import TSNE
from scipy.spatial import cKDTree

from scipy.stats import norm
from statsmodels.stats.multitest import multipletests

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

# =========================
# Parameters / configuration
# =========================

@dataclass
class DispersionPCARunParams:
    dist_threshold: float = 1.0
    n_perm: int = 10000
    fdr_alpha: float = 0.05

    # filters
    tissue_ids: Optional[List[str]] = None
    primary_cell_types: Optional[List[str]] = None
    neighbor_cell_types: Optional[List[str]] = None

    # TSNE options (sensitivity test)
    tsne_perplexity: int = 30
    tsne_random_state: Optional[int] = 0

    # Plot options (significant pairs plot)
    plot_n_cols: int = 4
    plot_max_panels: Optional[int] = None  # None => plot all
    figsize_per_panel: Tuple[float, float] = (4.2, 5.2)

    # Behavior
    skip_if_no_markers: bool = True  # for marker-removal sensitivity stage


# =========================
# IO helpers
# =========================

def read_genes_list(genes_names_path: Path) -> List[str]:
    genes: List[str] = []
    with open(genes_names_path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.append(g)
    return genes


def load_marker_map(marker_genes_by_tissue_json: Optional[Path]) -> Optional[Dict[str, str]]:
    """
    Expects JSON mapping tissue -> marker CSV path, like:
      {"100": "/path/to/markers.csv", "313": "/path/to/markers_313.csv", ...}
    """
    if marker_genes_by_tissue_json is None:
        return None

    with open(marker_genes_by_tissue_json, "r", encoding="utf-8") as f:
        d = json.load(f)

    out: Dict[str, str] = {}
    for k, v in d.items():
        out[str(k)] = "" if v is None else str(v)
    return out


def load_marker_genes_for_tissue(
    marker_genes_by_tissue: Optional[Mapping[str, str]],
    tissue: str,
) -> Optional[pd.DataFrame]:
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


def load_labels(summary_table_path: Path) -> pd.DataFrame:
    """
    Load the summary table and enforce:
      - 'Var1' exists and is used as the index (cell IDs)
      - 'tissue' is treated as a string identifier
    """
    labels_full = pd.read_csv(summary_table_path, header=0, dtype={"Var1": str})

    if "Var1" not in labels_full.columns:
        raise ValueError(
            f"summary_table_path must contain a 'Var1' column with unique cell IDs. "
            f"Missing 'Var1' in: {summary_table_path}"
        )

    # set Var1 as the index explicitly
    labels_full["Var1"] = labels_full["Var1"].astype(str)
    labels_full = labels_full.set_index("Var1", drop=True)
    labels_full.index = labels_full.index.astype(str)

    if "tissue" not in labels_full.columns:
        raise ValueError(
            f"summary_table_path must contain a 'tissue' column. "
            f"Missing 'tissue' in: {summary_table_path}"
        )

    # treat tissue as string IDs
    labels_full["tissue"] = labels_full["tissue"].astype(str)

    return labels_full


def load_distance_matrix(distance_dir: Path, tissue: str) -> pd.DataFrame:
    path = distance_dir / f"distance_matrix_{str(tissue)}.csv"
    dm = pd.read_csv(path, header=0, index_col=0, dtype={0: str})
    #dm.index = dm.index.astype(str)
    #dm.columns = dm.columns.astype(str)
    dm.index = dm.index.map(str)
    dm.columns = dm.columns.map(str)
    return dm


# =========================
# Filtering helpers
# =========================

def _normalize_cell_type_list(x: Optional[List[str]]) -> Optional[List[str]]:
    if x is None:
        return None
    out = [str(v).strip() for v in x if str(v).strip() != ""]
    return out if len(out) > 0 else None


def _filter_cell_types(available: List[str], allowed: Optional[List[str]]) -> List[str]:
    if allowed is None:
        return available
    allowed_set = set(allowed)
    return [ct for ct in available if ct in allowed_set]


# =========================
# Core proximity helper
# =========================

def get_proximity(
    distance_matrix: pd.DataFrame,
    primary_cells: pd.Index,
    neighbor_cells: pd.Index,
    threshold: float,
) -> pd.DataFrame:
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


def get_proximity_any_from_coordinates(
    labels_tissue: pd.DataFrame,
    ct1_index: pd.Index,          # ct1 cells used in embedding (emb_df.index)
    ct2_index: pd.Index,          # all ct2 cells in this tissue
    dist_threshold: float,
    x_col: str = "X_space",
    y_col: str = "Y_space",
    z_col: str = "Z_space",
) -> pd.Series:
    """
    Returns a boolean Series indexed by ct1_index:
      True  => min distance(ct1 cell -> any ct2 cell) <= dist_threshold
      False => otherwise (including missing coordinates)

    Uses 3D coordinates if z_col exists in labels_tissue; otherwise falls back to 2D.
    """

    # Decide whether to use 2D or 3D coordinates
    use_3d = z_col in labels_tissue.columns
    cols = [x_col, y_col, z_col] if use_3d else [x_col, y_col]

    # Build KDTree on ct2 coordinates (drop ct2 cells with missing coords)
    ct2_coords_df = (
        labels_tissue.loc[ct2_index, cols]
        .apply(pd.to_numeric, errors="coerce")
        .dropna(axis=0, how="any")
    )
    if ct2_coords_df.empty:
        # No valid ct2 coordinates => all ct1 are non-proximal
        return pd.Series(False, index=ct1_index.astype(str))

    tree = cKDTree(ct2_coords_df.values)

    # Prepare ct1 coordinates; keep only rows with complete coordinates
    ct1_coords_df = labels_tissue.loc[ct1_index, cols].apply(pd.to_numeric, errors="coerce")
    valid_mask = ~ct1_coords_df.isna().any(axis=1)
    valid_ct1 = ct1_coords_df.index[valid_mask]

    # Default: all False
    prox_any = pd.Series(False, index=ct1_index.astype(str))

    if len(valid_ct1) == 0:
        return prox_any

    # Query nearest ct2 neighbor for each ct1 cell
    dists, _ = tree.query(ct1_coords_df.loc[valid_ct1].values, k=1)

    # Threshold to create boolean proximity labels
    prox_any.loc[valid_ct1.astype(str)] = (dists <= dist_threshold)

    return prox_any


# =========================
# Statistics: PCA/tSNE + max-z Gaussian p-value
# =========================

def _compute_embedding_2d(
    raw_data: pd.DataFrame,
    method: str,
    tsne_perplexity: int,
    tsne_random_state: Optional[int],
) -> pd.DataFrame:
    X = StandardScaler().fit_transform(raw_data.values)

    if method == "pca":
        emb = PCA(n_components=2).fit_transform(X)
        return pd.DataFrame(emb, columns=["Dim1", "Dim2"], index=raw_data.index.astype(str))

    if method == "tsne":
        n_samples = X.shape[0]
        if n_samples < 3:
            raise ValueError("Too few samples for t-SNE (need >= 3).")
        perplexity = min(tsne_perplexity, n_samples - 1)

        emb = TSNE(
            n_components=2,
            perplexity=perplexity,
            random_state=tsne_random_state,
        ).fit_transform(X)
        return pd.DataFrame(emb, columns=["Dim1", "Dim2"], index=raw_data.index.astype(str))

    raise ValueError(f"Unknown method: {method}")


def _max_z_gaussian_pvalue_from_permutations(
    z_single: float,
    z_two: float,
    single_perm_distances: np.ndarray,
    two_perm_distances: np.ndarray,
) -> Tuple[float, float]:
    max_z_original = np.nanmax([z_single, z_two])
    if not np.isfinite(max_z_original):
        return np.nan, np.nan

    mean_single = float(np.nanmean(single_perm_distances))
    std_single = float(np.nanstd(single_perm_distances, ddof=0))
    mean_two = float(np.nanmean(two_perm_distances))
    std_two = float(np.nanstd(two_perm_distances, ddof=0))

    if std_single <= 0 or std_two <= 0:
        return np.nan, np.nan

    z_single_perm = (single_perm_distances - mean_single) / std_single
    z_two_perm = (two_perm_distances - mean_two) / std_two
    max_z_perm = np.nanmax(np.vstack([z_single_perm, z_two_perm]), axis=0)

    mean_max = float(np.nanmean(max_z_perm))
    std_max = float(np.nanstd(max_z_perm, ddof=0))
    if std_max <= 0:
        return np.nan, np.nan

    z_max_gaus = (max_z_original - mean_max) / std_max
    pval_max_gaus = float(norm.sf(z_max_gaus)) if np.isfinite(z_max_gaus) else np.nan
    return float(z_max_gaus), float(pval_max_gaus)


def _compute_pair_stats_with_permutations(
    emb_df: pd.DataFrame,
    proximity_any: pd.Series,
    n_perm: int,
    rng: np.random.Generator,
) -> Dict[str, Any]:
    df = emb_df.copy()
    df["proximity"] = proximity_any.reindex(df.index).fillna(False).astype(bool).values

    num_proximal = int(df["proximity"].sum())
    if num_proximal == 0 or num_proximal == len(df):
        raise ValueError("Only one proximity state (all proximal or all non-proximal).")

    df_true = df[df["proximity"]].copy()
    df_false = df[~df["proximity"]].copy()

    centroid_all = [df["Dim1"].mean(), df["Dim2"].mean()]
    centroid_true = [df_true["Dim1"].mean(), df_true["Dim2"].mean()]
    centroid_false = [df_false["Dim1"].mean(), df_false["Dim2"].mean()]

    distances_true = euclidean_distances(df_true[["Dim1", "Dim2"]], [centroid_all]).flatten()
    original_single = float(distances_true.sum())

    distances_all = euclidean_distances(df[["Dim1", "Dim2"]], [centroid_all]).flatten()
    original_all_single = float(distances_all.sum())

    original_two = float(euclidean_distances([centroid_true], [centroid_false])[0, 0])

    single_perm_list: List[float] = []
    two_perm_list: List[float] = []

    idx_all = df.index.to_numpy()
    for _ in range(int(n_perm)):
        sampled_idx = rng.choice(idx_all, size=num_proximal, replace=False)
        sampled_df = df.loc[sampled_idx]

        d_perm_single = euclidean_distances(sampled_df[["Dim1", "Dim2"]], [centroid_all]).flatten().sum()
        single_perm_list.append(float(d_perm_single))

        centroid_true_perm = [sampled_df["Dim1"].mean(), sampled_df["Dim2"].mean()]
        false_df_perm = df.loc[~df.index.isin(sampled_idx)]
        centroid_false_perm = [false_df_perm["Dim1"].mean(), false_df_perm["Dim2"].mean()]
        d_perm_two = float(euclidean_distances([centroid_true_perm], [centroid_false_perm])[0, 0])
        two_perm_list.append(float(d_perm_two))

    single_perm = np.asarray(single_perm_list, dtype=float)
    two_perm = np.asarray(two_perm_list, dtype=float)

    mean_single = float(np.mean(single_perm))
    std_single = float(np.std(single_perm, ddof=0))
    z_single = (original_single - mean_single) / std_single if std_single > 0 else np.nan

    mean_two = float(np.mean(two_perm))
    std_two = float(np.std(two_perm, ddof=0))
    z_two = (original_two - mean_two) / std_two if std_two > 0 else np.nan

    max_z_original = float(np.nanmax([z_single, z_two])) if (np.isfinite(z_single) or np.isfinite(z_two)) else np.nan

    max_z_gaus, max_z_pval_gaus = _max_z_gaussian_pvalue_from_permutations(
        z_single=float(z_single) if np.isfinite(z_single) else np.nan,
        z_two=float(z_two) if np.isfinite(z_two) else np.nan,
        single_perm_distances=single_perm,
        two_perm_distances=two_perm,
    )

    relative_original = original_single / original_all_single if original_all_single > 0 else np.nan
    rel_perm_vec = single_perm / original_all_single if original_all_single > 0 else np.full(len(single_perm), np.nan)
    mean_rel_perm = float(np.nanmean(rel_perm_vec))
    std_rel_perm = float(np.nanstd(rel_perm_vec, ddof=0))

    return {
        "num_cells_ct1_used": int(df.shape[0]),
        "num_proximal": int(num_proximal),

        "original_dist_two_centroids": float(original_two),
        "original_dist_single_centroid": float(original_single),
        "original_dist_all_to_centroid": float(original_all_single),

        "relative_original_dist_single_centroid": float(relative_original),
        "mean_relative_perm_dist_single_centroid": float(mean_rel_perm),
        "std_relative_perm_dist_single_centroid": float(std_rel_perm),

        "z_single": float(z_single) if np.isfinite(z_single) else np.nan,
        "z_two": float(z_two) if np.isfinite(z_two) else np.nan,

        "max_z_original": float(max_z_original) if np.isfinite(max_z_original) else np.nan,
        "max_z_gaus": float(max_z_gaus) if np.isfinite(max_z_gaus) else np.nan,
        "max_z_pval_gaus": float(max_z_pval_gaus) if np.isfinite(max_z_pval_gaus) else np.nan,
    }


# =========================
# Stage A: Full pairs computation (PCA space)
# =========================

def compute_all_pairs_pca(
    labels_full: pd.DataFrame,
    distance_dir: Optional[Path],
    genes_names: List[str],
    tissue: str,
    dist_threshold: float,
    n_perm: int,
    rng: np.random.Generator,
    primary_cell_types: Optional[List[str]] = None,
    neighbor_cell_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    labels = labels_full[labels_full["tissue"].astype(str) == str(tissue)].copy()
    labels.index = labels.index.astype(str)
    if labels.empty:
        return pd.DataFrame()

    dm: Optional[pd.DataFrame] = None
    if distance_dir is not None:
        dm_path = distance_dir / f"distance_matrix_{str(tissue)}.csv"
        if dm_path.exists():
            dm = load_distance_matrix(distance_dir, tissue)
        else:
            warnings.warn(f"[tissue {tissue}] distance matrix not found at: {dm_path}. Falling back to coordinates.")


    cell_types = labels["cell_type"].dropna().unique().tolist()
    cell_types = [ct for ct in cell_types if pd.notna(ct)]

    primary_cell_types = _normalize_cell_type_list(primary_cell_types)
    neighbor_cell_types = _normalize_cell_type_list(neighbor_cell_types)

    cell_types_primary = _filter_cell_types(cell_types, primary_cell_types)
    cell_types_neighbor = _filter_cell_types(cell_types, neighbor_cell_types)

    genes_missing = [g for g in genes_names if g not in labels.columns]
    if len(genes_missing) > 0:
        warnings.warn(f"[tissue {tissue}] {len(genes_missing)} genes from genes_names were not found in labels columns.")

    genes_existing = [g for g in genes_names if g in labels.columns]
    if len(genes_existing) == 0:
        raise ValueError(f"[tissue {tissue}] No overlap between genes_names and labels gene columns.")

    out_rows: List[Dict[str, Any]] = []

    for ct1 in cell_types_primary:
        ct1_cells = labels.index[labels["cell_type"] == ct1].astype(str)

        raw_data = labels.loc[ct1_cells, genes_existing].copy()
        raw_data = raw_data.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")
        if raw_data.shape[0] < 2:
            continue

        try:
            emb_df = _compute_embedding_2d(raw_data=raw_data, method="pca", tsne_perplexity=0, tsne_random_state=None)
        except Exception as e:
            warnings.warn(f"[tissue {tissue}] PCA failed for ct1={ct1}: {e}")
            continue

        for ct2 in cell_types_neighbor:
            if ct1 == ct2:
                continue

            ct2_cells = labels.index[labels["cell_type"] == ct2].astype(str)
            if len(ct2_cells) == 0:
                continue

            if dm is not None:
                try:
                    prox_matrix = get_proximity(dm, emb_df.index, ct2_cells, dist_threshold)
                except Exception:
                    continue
                if prox_matrix.shape[1] == 0:
                    continue
                proximity_any = prox_matrix.any(axis=1)
            else:
                proximity_any = get_proximity_any_from_coordinates(
                    labels_tissue=labels,
                    ct1_index=emb_df.index,
                    ct2_index=ct2_cells,
                    dist_threshold=dist_threshold,
                    x_col="X_space",
                    y_col="Y_space",
                    z_col="Z_space",
                )


            try:
                stats = _compute_pair_stats_with_permutations(
                    emb_df=emb_df,
                    proximity_any=proximity_any.astype(bool),
                    n_perm=n_perm,
                    rng=rng,
                )
            except Exception:
                continue

            out_rows.append({
                "tissue": str(tissue),
                "cell_type_1": str(ct1),
                "cell_type_2": str(ct2),
                **stats,
            })

    return pd.DataFrame(out_rows)


# =========================
# Stage B: Shuffled proximity sanity check (PCA space, all pairs)
# =========================

def compute_all_pairs_pca_shuffled_proximity(
    labels_full: pd.DataFrame,
    distance_dir: Optional[Path],
    genes_names: List[str],
    tissue: str,
    dist_threshold: float,
    n_perm: int,
    rng: np.random.Generator,
    primary_cell_types: Optional[List[str]] = None,
    neighbor_cell_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    labels = labels_full[labels_full["tissue"].astype(str) == str(tissue)].copy()
    labels.index = labels.index.astype(str)
    if labels.empty:
        return pd.DataFrame()

    dm: Optional[pd.DataFrame] = None
    if distance_dir is not None:
        dm_path = distance_dir / f"distance_matrix_{str(tissue)}.csv"
        if dm_path.exists():
            dm = load_distance_matrix(distance_dir, tissue)
        else:
            warnings.warn(f"[tissue {tissue}] distance matrix not found at: {dm_path}. Falling back to coordinates.")


    cell_types = labels["cell_type"].dropna().unique().tolist()
    cell_types = [ct for ct in cell_types if pd.notna(ct)]

    primary_cell_types = _normalize_cell_type_list(primary_cell_types)
    neighbor_cell_types = _normalize_cell_type_list(neighbor_cell_types)

    cell_types_primary = _filter_cell_types(cell_types, primary_cell_types)
    cell_types_neighbor = _filter_cell_types(cell_types, neighbor_cell_types)

    genes_existing = [g for g in genes_names if g in labels.columns]
    if len(genes_existing) == 0:
        return pd.DataFrame()

    out_rows: List[Dict[str, Any]] = []

    for ct1 in cell_types_primary:
        ct1_cells = labels.index[labels["cell_type"] == ct1].astype(str)

        raw_data = labels.loc[ct1_cells, genes_existing].copy()
        raw_data = raw_data.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")
        if raw_data.shape[0] < 2:
            continue

        try:
            emb_df = _compute_embedding_2d(raw_data=raw_data, method="pca", tsne_perplexity=0, tsne_random_state=None)
        except Exception:
            continue

        for ct2 in cell_types_neighbor:
            if ct1 == ct2:
                continue

            ct2_cells = labels.index[labels["cell_type"] == ct2].astype(str)
            if len(ct2_cells) == 0:
                continue

            if dm is not None:
                try:
                    prox_matrix = get_proximity(dm, emb_df.index, ct2_cells, dist_threshold)
                except Exception:
                    continue
                if prox_matrix.shape[1] == 0:
                    continue
                proximity_any = prox_matrix.any(axis=1)
            else:
                proximity_any = get_proximity_any_from_coordinates(
                    labels_tissue=labels,
                    ct1_index=emb_df.index,
                    ct2_index=ct2_cells,
                    dist_threshold=dist_threshold,
                    x_col="X_space",
                    y_col="Y_space",
                    z_col="Z_space",
                )


            shuffled = proximity_any.to_numpy().copy()
            rng.shuffle(shuffled)
            proximity_any_shuf = pd.Series(shuffled, index=proximity_any.index)

            try:
                stats = _compute_pair_stats_with_permutations(
                    emb_df=emb_df,
                    proximity_any=proximity_any_shuf.astype(bool),
                    n_perm=n_perm,
                    rng=rng,
                )
            except Exception:
                continue

            out_rows.append({
                "tissue": str(tissue),
                "cell_type_1": str(ct1),
                "cell_type_2": str(ct2),
                **stats,
            })

    return pd.DataFrame(out_rows)

# =========================
# Shuffled proximity sanity check (tSNE space, all pairs)
# =========================

def compute_all_pairs_tsne_shuffled_proximity(
    labels_full: pd.DataFrame,
    distance_dir: Optional[Path],
    genes_names: List[str],
    tissue: str,
    dist_threshold: float,
    n_perm: int,
    rng: np.random.Generator,
    tsne_perplexity: int,
    tsne_random_state: Optional[int],
    primary_cell_types: Optional[List[str]] = None,
    neighbor_cell_types: Optional[List[str]] = None,
) -> pd.DataFrame:
    labels = labels_full[labels_full["tissue"].astype(str) == str(tissue)].copy()
    labels.index = labels.index.astype(str)
    if labels.empty:
        return pd.DataFrame()

    dm: Optional[pd.DataFrame] = None
    if distance_dir is not None:
        dm_path = distance_dir / f"distance_matrix_{str(tissue)}.csv"
        if dm_path.exists():
            dm = load_distance_matrix(distance_dir, tissue)
        else:
            warnings.warn(f"[tissue {tissue}] distance matrix not found at: {dm_path}. Falling back to coordinates.")


    cell_types = labels["cell_type"].dropna().unique().tolist()
    cell_types = [ct for ct in cell_types if pd.notna(ct)]

    primary_cell_types = _normalize_cell_type_list(primary_cell_types)
    neighbor_cell_types = _normalize_cell_type_list(neighbor_cell_types)

    cell_types_primary = _filter_cell_types(cell_types, primary_cell_types)
    cell_types_neighbor = _filter_cell_types(cell_types, neighbor_cell_types)

    genes_existing = [g for g in genes_names if g in labels.columns]
    if len(genes_existing) < 2:
        return pd.DataFrame()

    out_rows: List[Dict[str, Any]] = []

    for ct1 in cell_types_primary:
        ct1_cells = labels.index[labels["cell_type"] == ct1].astype(str)

        raw_data = labels.loc[ct1_cells, genes_existing].copy()
        raw_data = raw_data.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")
        if raw_data.shape[0] < 3:  # t-SNE needs >=3
            continue

        # compute t-SNE ONCE per ct1
        try:
            emb_df = _compute_embedding_2d(
                raw_data=raw_data,
                method="tsne",
                tsne_perplexity=tsne_perplexity,
                tsne_random_state=tsne_random_state,
            )
        except Exception as e:
            warnings.warn(f"[tissue {tissue}] t-SNE failed for ct1={ct1}: {e}")
            continue

        for ct2 in cell_types_neighbor:
            if ct1 == ct2:
                continue

            ct2_cells = labels.index[labels["cell_type"] == ct2].astype(str)
            if len(ct2_cells) == 0:
                continue

            if dm is not None:
                try:
                    prox_matrix = get_proximity(dm, emb_df.index, ct2_cells, dist_threshold)
                except Exception:
                    continue
                if prox_matrix.shape[1] == 0:
                    continue
                proximity_any = prox_matrix.any(axis=1)
            else:
                proximity_any = get_proximity_any_from_coordinates(
                    labels_tissue=labels,
                    ct1_index=emb_df.index,
                    ct2_index=ct2_cells,
                    dist_threshold=dist_threshold,
                    x_col="X_space",
                    y_col="Y_space",
                    z_col="Z_space",
                )


            # ===== SHUFFLE proximity labels across ct1 cells (keep same #True/#False) =====
            arr = proximity_any.to_numpy().copy()
            rng.shuffle(arr)
            proximity_any_shuf = pd.Series(arr, index=proximity_any.index)

            try:
                stats = _compute_pair_stats_with_permutations(
                    emb_df=emb_df,
                    proximity_any=proximity_any_shuf.astype(bool),
                    n_perm=n_perm,
                    rng=rng,
                )
            except Exception:
                continue

            out_rows.append({
                "tissue": str(tissue),
                "cell_type_1": str(ct1),
                "cell_type_2": str(ct2),
                **stats,
            })

    return pd.DataFrame(out_rows)

# =========================
# Stage C: Marker-removal sensitivity (only significant pairs)
# =========================

def compute_sig_pairs_pca_marker_removed(
    sig_pairs_df: pd.DataFrame,
    labels_full: pd.DataFrame,
    distance_dir: Optional[Path],
    genes_names: List[str],
    marker_genes_by_tissue: Optional[Mapping[str, str]],
    dist_threshold: float,
    n_perm: int,
    rng: np.random.Generator,
    skip_if_no_markers: bool = True,
) -> pd.DataFrame:
    if sig_pairs_df.empty:
        return pd.DataFrame()

    sig = sig_pairs_df.copy()
    sig["tissue"] = sig["tissue"].astype(str).str.strip()
    sig["cell_type_1"] = sig["cell_type_1"].astype(str).str.strip()
    sig["cell_type_2"] = sig["cell_type_2"].astype(str).str.strip()

    dm_cache: Dict[str, Optional[pd.DataFrame]] = {}
    marker_cache: Dict[str, Optional[pd.DataFrame]] = {}
    out_rows: List[Dict[str, Any]] = []

    for _, row in sig.iterrows():
        tissue = str(row["tissue"])
        ct1 = str(row["cell_type_1"])
        ct2 = str(row["cell_type_2"])

        if tissue not in marker_cache:
            marker_cache[tissue] = load_marker_genes_for_tissue(marker_genes_by_tissue, tissue=tissue)
        mg = marker_cache[tissue]
        if mg is None and skip_if_no_markers:
            continue

        labels = labels_full[labels_full["tissue"].astype(str) == str(tissue)].copy()
        if labels.empty:
            continue
        labels.index = labels.index.astype(str)


        # Load distance matrix for this tissue if available; otherwise fall back to coordinates
        if tissue not in dm_cache:
            dm = None

            if distance_dir is not None:
                dm_path = distance_dir / f"distance_matrix_{str(tissue)}.csv"
                if dm_path.exists():
                    dm = load_distance_matrix(distance_dir, tissue)
                else:
                    warnings.warn(
                        f"[tissue {tissue}] distance matrix not found at: {dm_path}. "
                        f"Falling back to coordinate-based proximity."
                    )

            dm_cache[tissue] = dm

        dm = dm_cache[tissue]
        # NOTE: dm may legitimately be None here → handled later by coordinate-based fallback


        ct1_cells = labels.index[labels["cell_type"] == ct1].astype(str)
        ct2_cells = labels.index[labels["cell_type"] == ct2].astype(str)
        if len(ct1_cells) < 2 or len(ct2_cells) < 1:
            continue

        genes_existing = [g for g in genes_names if g in labels.columns]
        if len(genes_existing) < 2:
            continue

        if mg is not None:
            markers_to_remove = (
                mg.loc[mg["CellType"] == ct2, "Marker"]
                .dropna().astype(str).unique().tolist()
            )
            markers_set = set(markers_to_remove)
        else:
            markers_set = set()

        genes_filtered = [g for g in genes_existing if g not in markers_set]

        if len(genes_filtered) < 2:
            continue

        raw_data = labels.loc[ct1_cells, genes_filtered].copy()
        raw_data = raw_data.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")
        if raw_data.shape[0] < 2:
            continue

        try:
            emb_df = _compute_embedding_2d(raw_data=raw_data, method="pca", tsne_perplexity=0, tsne_random_state=None)
        except Exception:
            continue

        if dm is not None:
            # Distance-matrix–based proximity
            try:
                prox_matrix = get_proximity(dm, emb_df.index, ct2_cells, dist_threshold)
            except Exception:
                continue

            if prox_matrix.shape[1] == 0:
                continue

            proximity_any = prox_matrix.any(axis=1)

        else:
            # Coordinate-based fallback: min distance to ct2
            proximity_any = get_proximity_any_from_coordinates(
                labels_tissue=labels,        # already filtered to this tissue
                ct1_index=emb_df.index,
                ct2_index=ct2_cells,
                dist_threshold=dist_threshold,
                x_col="X_space",
                y_col="Y_space",
                z_col="Z_space",
            )

         
        try:
            stats = _compute_pair_stats_with_permutations(
                emb_df=emb_df,
                proximity_any=proximity_any.astype(bool),
                n_perm=n_perm,
                rng=rng,
            )
        except Exception:
            continue

        out_rows.append({
            "tissue": tissue,
            "cell_type_1": ct1,
            "cell_type_2": ct2,
            "num_genes_used_after_marker_removal": int(len(genes_filtered)),
            **stats,
        })

    return pd.DataFrame(out_rows)


# =========================
# Stage D: t-SNE sensitivity (only significant pairs)
# =========================

def compute_sig_pairs_tsne(
    sig_pairs_df: pd.DataFrame,
    labels_full: pd.DataFrame,
    distance_dir: Optional[Path],
    genes_names: List[str],
    dist_threshold: float,
    n_perm: int,
    rng: np.random.Generator,
    tsne_perplexity: int,
    tsne_random_state: Optional[int],
    shuffle_proximity_before_stats: bool = False,
) -> pd.DataFrame:
    if sig_pairs_df.empty:
        return pd.DataFrame()

    sig = sig_pairs_df.copy()
    sig["tissue"] = sig["tissue"].astype(str).str.strip()

    dm_cache: Dict[str, Optional[pd.DataFrame]] = {}
    out_rows: List[Dict[str, Any]] = []

    for _, row in sig.iterrows():
        tissue = str(row["tissue"])
        ct1 = str(row["cell_type_1"])
        ct2 = str(row["cell_type_2"])

        labels = labels_full[labels_full["tissue"].astype(str) == str(tissue)].copy()
        if labels.empty:
            continue
        labels.index = labels.index.astype(str)

        # Load distance matrix for this tissue if available; otherwise fall back to coordinates
        if tissue not in dm_cache:
            dm: Optional[pd.DataFrame] = None
            if distance_dir is not None:
                dm_path = distance_dir / f"distance_matrix_{str(tissue)}.csv"
                if dm_path.exists():
                    dm = load_distance_matrix(distance_dir, tissue)
                else:
                    warnings.warn(
                        f"[tissue {tissue}] distance matrix not found at: {dm_path}. "
                        f"Falling back to coordinate-based proximity."
                    )
            dm_cache[tissue] = dm

        dm = dm_cache[tissue]
        # NOTE: dm may be None here; handled by coordinate-based fallback below


        ct1_cells = labels.index[labels["cell_type"] == ct1].astype(str)
        ct2_cells = labels.index[labels["cell_type"] == ct2].astype(str)
        if len(ct1_cells) < 2 or len(ct2_cells) < 1:
            continue

        genes_existing = [g for g in genes_names if g in labels.columns]
        if len(genes_existing) < 2:
            continue

        raw_data = labels.loc[ct1_cells, genes_existing].copy()
        raw_data = raw_data.apply(pd.to_numeric, errors="coerce").dropna(axis=0, how="any")
        if raw_data.shape[0] < 2:
            continue

        try:
            emb_df = _compute_embedding_2d(
                raw_data=raw_data,
                method="tsne",
                tsne_perplexity=tsne_perplexity,
                tsne_random_state=tsne_random_state,
            )
        except Exception:
            continue


        if dm is not None:
            try:
                prox_matrix = get_proximity(dm, emb_df.index, ct2_cells, dist_threshold)
            except Exception:
                continue
            if prox_matrix.shape[1] == 0:
                continue
            proximity_any = prox_matrix.any(axis=1)
        else:
            proximity_any = get_proximity_any_from_coordinates(
                labels_tissue=labels,
                ct1_index=emb_df.index,
                ct2_index=ct2_cells,
                dist_threshold=dist_threshold,
                x_col="X_space",
                y_col="Y_space",
                z_col="Z_space",
            )


        if shuffle_proximity_before_stats:
            arr = proximity_any.to_numpy().copy()
            rng.shuffle(arr)
            proximity_any = pd.Series(arr, index=proximity_any.index)

        try:
            stats = _compute_pair_stats_with_permutations(
                emb_df=emb_df,
                proximity_any=proximity_any.astype(bool),
                n_perm=n_perm,
                rng=rng,
            )
        except Exception:
            continue

        out_rows.append({
            "tissue": tissue,
            "cell_type_1": ct1,
            "cell_type_2": ct2,
            **stats,
        })

    return pd.DataFrame(out_rows)

# =========================
#   plot significant pairs
# =========================
def plot_sig_pairs_pca_panels(
    sig_pairs_df: pd.DataFrame,
    labels_full: pd.DataFrame,
    distance_dir: Optional[Path],
    genes_names: List[str],
    dist_threshold: float,
    out_path: Path,
    max_cols: int = 4,                 # unused, kept for compatibility
    max_panels: Optional[int] = None,  # cap number of figures
    panel_size: float = 4.0,
    point_size: float = 18.0,
) -> None:
    """
    Create ONE standard matplotlib figure per significant pair (PNG).

    - Single figure, single Axes
    - Standard ax.set_title (no suptitle)
    - Standard ax.legend (no bbox_to_anchor)
    - Saved under: <out_path.parent>/plots/{tissue}_{ct1}_{ct2}.png
    """

    sig = sig_pairs_df.copy()
    if sig.empty:
        print("No significant pairs to plot.")
        return

    if max_panels is not None:
        sig = sig.iloc[: int(max_panels)].copy()

    sig["tissue"] = sig["tissue"].astype(str).str.strip()
    sig["cell_type_1"] = sig["cell_type_1"].astype(str).str.strip()
    sig["cell_type_2"] = sig["cell_type_2"].astype(str).str.strip()

    plots_dir = Path(out_path).parent / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    dm_cache: Dict[str, Optional[pd.DataFrame]] = {}

    def _safe(s: str) -> str:
        return (
            str(s)
            .replace(" ", "_")
            .replace("/", "_")
            .replace("\\", "_")
            .replace(":", "_")
        )

    saved = 0

    for _, row in sig.iterrows():
        tissue = row["tissue"]
        ct1 = row["cell_type_1"]
        ct2 = row["cell_type_2"]

        labels = labels_full[labels_full["tissue"].astype(str) == tissue].copy()
        labels.index = labels.index.astype(str)
        if labels.empty:
            continue

        ct1_cells = labels.index[labels["cell_type"] == ct1]
        ct2_cells = labels.index[labels["cell_type"] == ct2]
        if len(ct1_cells) < 2 or len(ct2_cells) < 1:
            continue

        genes_existing = [g for g in genes_names if g in labels.columns]
        if len(genes_existing) < 2:
            continue

        raw_data = (
            labels.loc[ct1_cells, genes_existing]
            .apply(pd.to_numeric, errors="coerce")
            .dropna(axis=0, how="any")
        )
        if raw_data.shape[0] < 2:
            continue

        # PCA
        X = StandardScaler().fit_transform(raw_data.values)
        pcs = PCA(n_components=2).fit_transform(X)
        pca_df = pd.DataFrame(pcs, columns=["PC1", "PC2"], index=raw_data.index)

        # distance matrix (cached per tissue)
        if tissue not in dm_cache:
            dm = None
            if distance_dir is not None:
                dm_path = Path(distance_dir) / f"distance_matrix_{tissue}.csv"
                if dm_path.exists():
                    dm = load_distance_matrix(Path(distance_dir), tissue)
            dm_cache[tissue] = dm
        dm = dm_cache[tissue]

        # proximity labels
        if dm is not None:
            try:
                prox_matrix = get_proximity(
                    dm, pca_df.index, ct2_cells, float(dist_threshold)
                )
                proximity_any = prox_matrix.any(axis=1).reindex(pca_df.index).fillna(False)
            except Exception:
                continue
        else:
            try:
                proximity_any = get_proximity_any_from_coordinates(
                    labels_tissue=labels,
                    ct1_index=pca_df.index,
                    ct2_index=ct2_cells,
                    dist_threshold=float(dist_threshold),
                )
            except Exception:
                continue

        if proximity_any.nunique() < 2:
            continue

        is_prox = proximity_any.astype(bool).values

        # ----------------
        # Standard figure
        # ----------------
        fig, ax = plt.subplots(figsize=(panel_size, panel_size))

        ax.scatter(
            pca_df.loc[~is_prox, "PC1"],
            pca_df.loc[~is_prox, "PC2"],
            s=point_size,
            c="blue",
            edgecolors="black",
            linewidths=0.4,
            alpha=0.75,
            label="Distant cell",
        )
        ax.scatter(
            pca_df.loc[is_prox, "PC1"],
            pca_df.loc[is_prox, "PC2"],
            s=point_size,
            c="red",
            edgecolors="black",
            linewidths=0.4,
            alpha=0.75,
            label="Proximal cell",
        )

        cen_all = pca_df[["PC1", "PC2"]].mean().values
        cen_true = pca_df.loc[is_prox, ["PC1", "PC2"]].mean().values
        cen_false = pca_df.loc[~is_prox, ["PC1", "PC2"]].mean().values

        ax.scatter(
            cen_all[0], cen_all[1],
            s=55, c="yellow", edgecolors="black",
            linewidths=0.6, zorder=5,
            label="Overall centroid",
        )
        ax.scatter(
            [cen_true[0], cen_false[0]],
            [cen_true[1], cen_false[1]],
            s=55, c="green", edgecolors="black",
            linewidths=0.6, zorder=5,
            label="Proximal & distant centroids",
        )

        ax.set_title(f"Tissue: {tissue} | {ct1} vs {ct2}", fontsize=14, pad=10)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.grid(True, alpha=0.6)

        ax.legend(
            loc="upper left",
            bbox_to_anchor=(1.02, 1.0),
            borderaxespad=0.0,
            fontsize=11,
        )


        if hasattr(ax, "set_box_aspect"):
            ax.set_box_aspect(1)

        out_png = plots_dir / f"{_safe(tissue)}_{_safe(ct1)}_{_safe(ct2)}.png"
        fig.savefig(out_png, dpi=300, bbox_inches="tight")
        plt.close(fig)

        saved += 1

    print(f"[plot_sig_pairs_pca_panels] saved {saved} figures to {plots_dir}")

# =========================
# Orchestration: run across tissues + save outputs
# =========================

def run_dispersion_in_pca_space(
    summary_table_path: Path,
    genes_names_path: Path,
    distance_dir: Optional[Path],
    out_dir: Path,
    params: DispersionPCARunParams,
    marker_genes_by_tissue_json: Optional[Path] = None,
    rng_seed: Optional[int] = None,
) -> Dict[str, Any]:
    out_dir.mkdir(parents=True, exist_ok=True)

    genes_names = read_genes_list(genes_names_path)
    labels_full = load_labels(summary_table_path)
    marker_map = load_marker_map(marker_genes_by_tissue_json)

    # Enforce coordinate columns when no distance matrices are provided
    if distance_dir is None:
        missing = [c for c in ["X_space", "Y_space"] if c not in labels_full.columns]
        if missing:
            raise ValueError(
                "When --distance-dir is not provided, proximity is computed from coordinates and the summary table "
                "must include 'X_space' and 'Y_space' columns "
                f"(missing: {missing})."
            )


    # tissues: default all; if params.tissue_ids provided => filter
    all_tissues = list(labels_full["tissue"].unique().astype(str))
    if params.tissue_ids is not None and len(params.tissue_ids) > 0:
        allowed = set([str(x).strip() for x in params.tissue_ids])
        selected_tissues = [t for t in all_tissues if t in allowed]
    else:
        selected_tissues = all_tissues

    rng = np.random.default_rng(rng_seed)

    print(f"[START] Dispersion-PCA | tissues={selected_tissues} | primary={params.primary_cell_types} | neighbor={params.neighbor_cell_types}")
    print(f"[PARAMS] dist_threshold={params.dist_threshold} | n_perm={params.n_perm} | fdr_alpha={params.fdr_alpha}")

    # ---- 1) Full PCA all-pairs ----
    all_pairs_list: List[pd.DataFrame] = []
    for i, tissue in enumerate(selected_tissues, start=1):
        print(f"[STAGE 1/6] PCA all-pairs | tissue {tissue} ({i}/{len(selected_tissues)})")
        df_t = compute_all_pairs_pca(
            labels_full=labels_full,
            distance_dir=distance_dir,
            genes_names=genes_names,
            tissue=str(tissue),
            dist_threshold=params.dist_threshold,
            n_perm=params.n_perm,
            rng=rng,
            primary_cell_types=params.primary_cell_types,
            neighbor_cell_types=params.neighbor_cell_types,
        )
        if not df_t.empty:
            all_pairs_list.append(df_t.reset_index(drop=True))
        print(f"  -> pairs computed: {0 if df_t.empty else int(df_t.shape[0])}")

    df_all_pairs = pd.concat(all_pairs_list, ignore_index=True) if all_pairs_list else pd.DataFrame()

    if not df_all_pairs.empty and "max_z_pval_gaus" in df_all_pairs.columns:
        df_all_pairs["pval_gaus_fdr_bh"] = multipletests(df_all_pairs["max_z_pval_gaus"], method="fdr_bh")[1]
        df_all_pairs = df_all_pairs.sort_values("pval_gaus_fdr_bh", ascending=True, ignore_index=True)

    all_pairs_csv = out_dir / "pca_pairs_all_results.csv"
    df_all_pairs.reset_index(drop=True).to_csv(all_pairs_csv, index=False)

    if not df_all_pairs.empty and "pval_gaus_fdr_bh" in df_all_pairs.columns:
        df_sig = df_all_pairs[df_all_pairs["pval_gaus_fdr_bh"] <= params.fdr_alpha].copy()
    else:
        df_sig = pd.DataFrame()

    sig_csv = out_dir / "pca_pairs_significant.csv"
    df_sig.reset_index(drop=True).to_csv(sig_csv, index=False)

    print(f"[RESULT] PCA all-pairs: total={int(df_all_pairs.shape[0]) if not df_all_pairs.empty else 0} | sig(FDR)={int(df_sig.shape[0]) if not df_sig.empty else 0}")

    # ---- 2) Plot significant PCA panels ----
    print("[STAGE 2/6] Plot significant pairs (PCA panels)")
    fig_path = out_dir / "pca_significant_pairs_panels.tif"
    plot_sig_pairs_pca_panels(
        sig_pairs_df=df_sig,
        labels_full=labels_full,
        distance_dir=distance_dir,
        genes_names=genes_names,
        dist_threshold=params.dist_threshold,
        out_path=fig_path,
        max_cols=params.plot_n_cols,
        max_panels=params.plot_max_panels,
        panel_size=float(np.mean(params.figsize_per_panel))
    )

    # ---- 3) Shuffled proximity all-pairs (PCA) ----
    print("[STAGE 3/6] PCA shuffled-proximity (all pairs)")
    shuf_list: List[pd.DataFrame] = []
    for i, tissue in enumerate(selected_tissues, start=1):
        print(f"  tissue {tissue} ({i}/{len(selected_tissues)})")
        df_t = compute_all_pairs_pca_shuffled_proximity(
            labels_full=labels_full,
            distance_dir=distance_dir,
            genes_names=genes_names,
            tissue=str(tissue),
            dist_threshold=params.dist_threshold,
            n_perm=params.n_perm,
            rng=rng,
            primary_cell_types=params.primary_cell_types,
            neighbor_cell_types=params.neighbor_cell_types,
        )
        if not df_t.empty:
            shuf_list.append(df_t.reset_index(drop=True))
        print(f"    -> pairs computed: {0 if df_t.empty else int(df_t.shape[0])}")

    df_shuf = pd.concat(shuf_list, ignore_index=True) if shuf_list else pd.DataFrame()
    if not df_shuf.empty and "max_z_pval_gaus" in df_shuf.columns:
        df_shuf["pval_gaus_fdr_bh"] = multipletests(df_shuf["max_z_pval_gaus"], method="fdr_bh")[1]
        df_shuf = df_shuf.sort_values("pval_gaus_fdr_bh", ascending=True, ignore_index=True)

    shuf_csv = out_dir / "pca_pairs_all_results_shuffled_proximity.csv"
    df_shuf.reset_index(drop=True).to_csv(shuf_csv, index=False)

    n_sig_shuf = int(df_shuf[df_shuf["pval_gaus_fdr_bh"] <= params.fdr_alpha].shape[0]) if (not df_shuf.empty and "pval_gaus_fdr_bh" in df_shuf.columns) else 0
    print(f"[RESULT] PCA shuffled: total={int(df_shuf.shape[0]) if not df_shuf.empty else 0} | sig(FDR)={n_sig_shuf}")

    # ---- 4) Marker-removal sensitivity on significant pairs ----
    print("[STAGE 4/6] Marker-removal sensitivity (on significant pairs)")
    df_marker = compute_sig_pairs_pca_marker_removed(
        sig_pairs_df=df_sig,
        labels_full=labels_full,
        distance_dir=distance_dir,
        genes_names=genes_names,
        marker_genes_by_tissue=marker_map,
        dist_threshold=params.dist_threshold,
        n_perm=params.n_perm,
        rng=rng,
        skip_if_no_markers=params.skip_if_no_markers,
    )
    if not df_marker.empty and "max_z_pval_gaus" in df_marker.columns:
        df_marker["pval_gaus_fdr_bh"] = multipletests(df_marker["max_z_pval_gaus"], method="fdr_bh")[1]
        df_marker = df_marker.sort_values("pval_gaus_fdr_bh", ascending=True, ignore_index=True)

    marker_csv = out_dir / "pca_pairs_sig_marker_removal_sensitivity.csv"
    df_marker.reset_index(drop=True).to_csv(marker_csv, index=False)

    n_sig_marker = int(df_marker[df_marker["pval_gaus_fdr_bh"] <= params.fdr_alpha].shape[0]) if (not df_marker.empty and "pval_gaus_fdr_bh" in df_marker.columns) else 0
    print(f"[RESULT] Marker-removal: total={int(df_marker.shape[0]) if not df_marker.empty else 0} | sig(FDR)={n_sig_marker}")

    # ---- 5) t-SNE sensitivity on significant pairs ----
    print("[STAGE 5/6] t-SNE sensitivity (on significant pairs)")
    df_tsne = compute_sig_pairs_tsne(
        sig_pairs_df=df_sig,
        labels_full=labels_full,
        distance_dir=distance_dir,
        genes_names=genes_names,
        dist_threshold=params.dist_threshold,
        n_perm=params.n_perm,
        rng=rng,
        tsne_perplexity=params.tsne_perplexity,
        tsne_random_state=params.tsne_random_state,
        shuffle_proximity_before_stats=False,
    )
    if not df_tsne.empty and "max_z_pval_gaus" in df_tsne.columns:
        df_tsne["pval_gaus_fdr_bh"] = multipletests(df_tsne["max_z_pval_gaus"], method="fdr_bh")[1]
        df_tsne = df_tsne.sort_values("pval_gaus_fdr_bh", ascending=True, ignore_index=True)

    tsne_csv = out_dir / "tsne_pairs_sig_sensitivity.csv"
    df_tsne.reset_index(drop=True).to_csv(tsne_csv, index=False)

    n_sig_tsne = int(df_tsne[df_tsne["pval_gaus_fdr_bh"] <= params.fdr_alpha].shape[0]) if (not df_tsne.empty and "pval_gaus_fdr_bh" in df_tsne.columns) else 0
    print(f"[RESULT] t-SNE: total={int(df_tsne.shape[0]) if not df_tsne.empty else 0} | sig(FDR)={n_sig_tsne}")

    # ---- 6) t-SNE + shuffled proximity (ALL PAIRS) ----
    print("[STAGE 6/6] t-SNE + shuffled proximity (all pairs)")
    tsne_shuf_list: List[pd.DataFrame] = []
    for i, tissue in enumerate(selected_tissues, start=1):
        print(f"  tissue {tissue} ({i}/{len(selected_tissues)})")
        df_t = compute_all_pairs_tsne_shuffled_proximity(
            labels_full=labels_full,
            distance_dir=distance_dir,
            genes_names=genes_names,
            tissue=str(tissue),
            dist_threshold=params.dist_threshold,
            n_perm=params.n_perm,
            rng=rng,
            tsne_perplexity=params.tsne_perplexity,
            tsne_random_state=params.tsne_random_state,
            primary_cell_types=params.primary_cell_types,
            neighbor_cell_types=params.neighbor_cell_types,
        )
        if not df_t.empty:
            tsne_shuf_list.append(df_t.reset_index(drop=True))
        print(f"    -> pairs computed: {0 if df_t.empty else int(df_t.shape[0])}")

    df_tsne_shuf = pd.concat(tsne_shuf_list, ignore_index=True) if tsne_shuf_list else pd.DataFrame()
    if not df_tsne_shuf.empty and "max_z_pval_gaus" in df_tsne_shuf.columns:
        df_tsne_shuf["pval_gaus_fdr_bh"] = multipletests(df_tsne_shuf["max_z_pval_gaus"], method="fdr_bh")[1]
        df_tsne_shuf = df_tsne_shuf.sort_values("pval_gaus_fdr_bh", ascending=True, ignore_index=True)

    tsne_shuf_csv = out_dir / "tsne_pairs_all_sensitivity_shuffled_proximity.csv"
    df_tsne_shuf.reset_index(drop=True).to_csv(tsne_shuf_csv, index=False)

    n_sig_tsne_shuf = int(df_tsne_shuf[df_tsne_shuf["pval_gaus_fdr_bh"] <= params.fdr_alpha].shape[0]) if (not df_tsne_shuf.empty and "pval_gaus_fdr_bh" in df_tsne_shuf.columns) else 0
    print(f"[RESULT] t-SNE+shuffle (all pairs): total={int(df_tsne_shuf.shape[0]) if not df_tsne_shuf.empty else 0} | sig(FDR)={n_sig_tsne_shuf}")

    # ---- Summary counts table ----
    summary = {
        "n_total_pairs_pca_all": int(df_all_pairs.shape[0]) if not df_all_pairs.empty else 0,
        "n_sig_pairs_pca_all_fdr": int(df_sig.shape[0]) if not df_sig.empty else 0,
        "n_sig_pairs_shuffled_fdr": int(n_sig_shuf),
        "n_sig_marker_removal_fdr": int(n_sig_marker),
        "n_sig_tsne_fdr": int(n_sig_tsne),
        "n_sig_tsne_shuffle_allpairs_fdr": int(n_sig_tsne_shuf),
    }
    summary_df = pd.DataFrame([summary])
    summary_csv = out_dir / "dispersion_pca_summary_counts.csv"
    summary_df.to_csv(summary_csv, index=False)

    print(f"[DONE] Outputs written to: {out_dir}")
    print(f"  - all pairs: {all_pairs_csv}")
    print(f"  - significant pairs: {sig_csv}")
    print(f"  - significant plot: {fig_path}")
    print(f"  - shuffled all pairs: {shuf_csv}")
    print(f"  - marker sensitivity: {marker_csv}")
    print(f"  - tSNE sensitivity: {tsne_csv}")
    print(f"  - tSNE shuffled (all_pairs) sensitivity: {tsne_shuf_csv}")
    print(f"  - summary: {summary_csv}")

    return {
        "all_pairs_csv": str(all_pairs_csv),
        "sig_pairs_csv": str(sig_csv),
        "sig_plot_path": str(fig_path),
        "shuffled_pairs_csv": str(shuf_csv),
        "marker_sensitivity_csv": str(marker_csv),
        "tsne_sensitivity_csv": str(tsne_csv),
        "tsne_sensitivity_shuffled_csv": str(tsne_shuf_csv),
        "summary_csv": str(summary_csv),
        "n_all_pairs": int(df_all_pairs.shape[0]) if not df_all_pairs.empty else 0,
        "n_sig_pairs": int(df_sig.shape[0]) if not df_sig.empty else 0,
    }
