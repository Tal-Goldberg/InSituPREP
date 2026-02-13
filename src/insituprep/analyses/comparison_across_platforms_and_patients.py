from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Dict, Any, Tuple

import warnings
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import euclidean_distances


# =========================
# Params
# =========================

@dataclass
class ComparisonAcrossPlatformsRunParams:
    tissues: List[str]
    non_tumor_cell_types: List[str]

    # required in CLI (you want it required there)
    tumor_cell_type: str = "Epithelial"
    platform_prefix: str = "MERFISH_"

    # receptor status labelling
    hr_positive_tissues: Optional[List[str]] = None  # tissue names as they appear in pca_df["Tissue"]

    # stats / thresholds
    n_perm_pca: int = 10000
    n_perm_dist: int = 10000
    min_pv_for_tissue_filter: float = 0.01
    pc_pvalue_alpha: float = 0.05

    # progress printing
    verbose: bool = True
    perm_progress_every: int = 1000  # print every N permutations

    # read p-values from final_genes_* files
    # NOTE: you later add "*.csv" in glob, so template should NOT include ".csv"
    final_genes_prefix_template: str = "final_genes_{tissue}_{primary_cell_type}_vs_{neighbor_cell_type}_LFC"
    pval_column: str = "pvalue"  # raw pvalue
    gene_column: str = "gene"

    # plotting
    figsize: Tuple[float, float] = (12, 4)


# =========================
# Helpers
# =========================

def _try_import_adjust_text():
    try:
        from adjustText import adjust_text  # type: ignore
        return adjust_text
    except Exception:
        return None


def _normalize_cell_type_for_filename(s: str) -> str:
    # robust: collapses multiple spaces, trims, and replaces spaces with underscores
    return "_".join(str(s).strip().split())


def _ordered_base_tissues(tissues: List[str], platform_prefix: str) -> List[str]:
    base_tissues: List[str] = []
    seen = set()
    for t in tissues:
        base = str(t).replace(platform_prefix, "")
        if base not in seen:
            seen.add(base)
            base_tissues.append(base)
    return base_tissues


def _tissue_color_map(tissues: List[str], platform_prefix: str) -> Dict[str, Any]:
    base_tissues = _ordered_base_tissues(tissues, platform_prefix)
    base_colors = sns.color_palette("tab10", len(base_tissues))
    base_color_map = dict(zip(base_tissues, base_colors))
    return {str(t): base_color_map[str(t).replace(platform_prefix, "")] for t in tissues}


def _rescale_to_original_range(x: pd.Series, orig_min: float, orig_max: float) -> pd.Series:
    x_min, x_max = float(x.min()), float(x.max())
    if x_max == x_min:
        return pd.Series(np.full(shape=len(x), fill_value=(orig_min + orig_max) / 2), index=x.index)
    return ((x - x_min) / (x_max - x_min)) * (orig_max - orig_min) + orig_min


def _pc_permutation_pvalues(
    scaled_matrix: np.ndarray,
    real_var: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
    *,
    verbose: bool = False,
    progress_every: int = 1000,
    label: str = "",
) -> np.ndarray:
    null_variance_ratios = np.zeros((n_perm, 2))
    for i in range(n_perm):
        perm = np.copy(scaled_matrix)
        for col in range(perm.shape[1]):
            rng.shuffle(perm[:, col])
        pca_perm = PCA(n_components=2)
        pca_perm.fit(perm)
        null_variance_ratios[i] = pca_perm.explained_variance_ratio_

        if verbose and progress_every > 0 and (i + 1) % progress_every == 0:
            print(f"[perm PC] {label}  {i+1}/{n_perm}")

    pvals = (np.sum(null_variance_ratios >= real_var, axis=0) + 1) / n_perm
    return pvals


def _load_pval_matrix_from_final_genes(
    results_dir: Path,
    tissues: List[str],
    primary_cell_type: str,
    neighbor_cell_type: str,
    *,
    file_prefix_template: str,
    gene_column: str,
    pval_column: str,
) -> Optional[pd.DataFrame]:
    """
    Loads p-values from files named like:
      final_genes_{tissue}_{primary_cell_type}_vs_{neighbor_cell_type}_LFC*.csv
    Returns:
      DataFrame index=gene, columns=tissues, INNER join across tissues (shared genes only).
    """
    pval_df: Optional[pd.DataFrame] = None

    for tissue in tissues:
        prefix = file_prefix_template.format(
            tissue=str(tissue),
            primary_cell_type=primary_cell_type,
            neighbor_cell_type=neighbor_cell_type,
        )
        matches = sorted(results_dir.glob(prefix + "*.csv"))
        if len(matches) == 0:
            continue

        file_path = matches[0]  # deterministic first match

        df = pd.read_csv(file_path)
        if gene_column not in df.columns:
            raise ValueError(f"Missing gene column '{gene_column}' in {file_path}")
        if pval_column not in df.columns:
            raise ValueError(
                f"Missing p-value column '{pval_column}' in {file_path}. "
                f"Available columns: {list(df.columns)}"
            )

        df_sub = df[[gene_column, pval_column]].dropna()
        df_sub = df_sub.set_index(gene_column)
        df_sub.columns = [str(tissue)]

        if pval_df is None:
            pval_df = df_sub
        else:
            pval_df = pval_df.join(df_sub, how="inner")

    return pval_df


def _platform_pair_distance_pvalues(
    pval_df_n: pd.DataFrame,  # genes x tissues  (-log10)
    pca_df: pd.DataFrame,     # Tissue, PC1, PC2
    n_perms: int,
    rng: np.random.Generator,
    *,
    verbose: bool = False,
    progress_every: int = 1000,
    label: str = "",
) -> pd.DataFrame:
    """
    p-value for ExSeq vs MERFISH distance (paired tissue id) in PCA space.
    p = P(perm_dist <= original_dist)
    """
    pca_df = pca_df.copy()
    pca_df["Tissue_stripped"] = pca_df["Tissue"].astype(str).str.extract(r"(\d+)$")[0]
    grouped = pca_df.groupby("Tissue_stripped")

    pc1_min, pc1_max = float(pca_df["PC1"].min()), float(pca_df["PC1"].max())
    pc2_min, pc2_max = float(pca_df["PC2"].min()), float(pca_df["PC2"].max())

    pairs_pv_df = pd.DataFrame(columns=["Tissue_ID", "p_value_dist"])
    pval_df_for_perm = pval_df_n.T  # tissues x genes

    for tissue_id, group in grouped:
        if len(group) != 2:
            continue

        tissue1 = str(group["Tissue"].iloc[0])
        tissue2 = str(group["Tissue"].iloc[1])

        original_coords = group[["PC1", "PC2"]].values
        original_dist = float(euclidean_distances([original_coords[0]], [original_coords[1]])[0, 0])

        if verbose:
            print(f"[dist] {label}  tissue_pair={tissue_id}  ({tissue1} vs {tissue2})  perms={n_perms}")

        perm_distances: List[float] = []
        for k in range(n_perms):
            permuted = pval_df_for_perm.copy()
            for col in permuted.columns:
                permuted[col] = rng.permutation(permuted[col].values)

            scaled = StandardScaler().fit_transform(permuted)
            pca = PCA(n_components=2)
            coords = pca.fit_transform(scaled)
            pca_df_perm = pd.DataFrame(coords, columns=["PC1", "PC2"], index=permuted.index)

            pca_df_perm["PC1"] = _rescale_to_original_range(pca_df_perm["PC1"], pc1_min, pc1_max)
            pca_df_perm["PC2"] = _rescale_to_original_range(pca_df_perm["PC2"], pc2_min, pc2_max)

            perm_dist = float(euclidean_distances([pca_df_perm.loc[tissue1]], [pca_df_perm.loc[tissue2]])[0, 0])
            perm_distances.append(perm_dist)

            if verbose and progress_every > 0 and (k + 1) % progress_every == 0:
                print(f"[perm dist] {label}  tissue_pair={tissue_id}  {k+1}/{n_perms}")

        perm_distances_arr = np.asarray(perm_distances)
        p_value = (np.sum(perm_distances_arr <= original_dist) + 1) / n_perms

        pairs_pv_df = pd.concat(
            [pairs_pv_df, pd.DataFrame([{"Tissue_ID": tissue_id, "p_value_dist": float(p_value)}])],
            ignore_index=True,
        )

    return pairs_pv_df


def _exseq_hr_proximity_pvalue(
    pval_df_n: pd.DataFrame,   # genes x tissues (-log10)
    pca_df: pd.DataFrame,      # Tissue, PC1, PC2, ReceptorStatus
    min_max_df: pd.DataFrame,  # min/max of raw pvalues per tissue
    n_perms: int,
    platform_prefix: str,
    rng: np.random.Generator,
    *,
    verbose: bool = False,
    progress_every: int = 1000,
    label: str = "",
) -> Optional[float]:
    """
    Mean pairwise PCA distance among ExSeq HR+/HER2- tissues.
    p = P(perm_mean_dist <= original_mean_dist)
    """
    significant_tissues = min_max_df[min_max_df["min"] <= 0.01].index
    pca_df_filtered = pca_df[pca_df["Tissue"].isin(significant_tissues)]

    pca_df_receptor = pca_df_filtered[pca_df_filtered["ReceptorStatus"] == "HR+/HER2-"]
    exseq_tissues = [t for t in pca_df_receptor["Tissue"] if not str(t).startswith(platform_prefix)]
    if len(exseq_tissues) <= 1:
        return None

    if verbose:
        print(f"[receptor] {label}  exseq_hr_tissues={len(exseq_tissues)}  perms={n_perms}")

    exseq_coords = pca_df[pca_df["Tissue"].isin(exseq_tissues)][["PC1", "PC2"]].reset_index(drop=True)

    from itertools import combinations
    pairs = list(combinations(range(len(exseq_coords)), 2))
    original_dists = [
        float(euclidean_distances([exseq_coords.iloc[i]], [exseq_coords.iloc[j]])[0, 0]) for i, j in pairs
    ]
    original_mean = float(np.mean(original_dists))

    pc1_min, pc1_max = float(pca_df["PC1"].min()), float(pca_df["PC1"].max())
    pc2_min, pc2_max = float(pca_df["PC2"].min()), float(pca_df["PC2"].max())

    pval_df_for_perm = pval_df_n.T  # tissues x genes

    perm_means: List[float] = []
    for k in range(n_perms):
        permuted = pval_df_for_perm.copy()
        for col in permuted.columns:
            permuted[col] = rng.permutation(permuted[col].values)

        scaled = StandardScaler().fit_transform(permuted)
        pca = PCA(n_components=2)
        coords = pca.fit_transform(scaled)
        pca_df_perm = pd.DataFrame(coords, columns=["PC1", "PC2"], index=permuted.index)

        pca_df_perm["PC1"] = _rescale_to_original_range(pca_df_perm["PC1"], pc1_min, pc1_max)
        pca_df_perm["PC2"] = _rescale_to_original_range(pca_df_perm["PC2"], pc2_min, pc2_max)

        perm_exseq_coords = pca_df_perm.loc[exseq_tissues][["PC1", "PC2"]].reset_index(drop=True)
        perm_dists = [
            float(euclidean_distances([perm_exseq_coords.iloc[i]], [perm_exseq_coords.iloc[j]])[0, 0])
            for i, j in pairs
        ]
        perm_means.append(float(np.mean(perm_dists)))

        if verbose and progress_every > 0 and (k + 1) % progress_every == 0:
            print(f"[perm receptor] {label}  {k+1}/{n_perms}")

    perm_means_arr = np.asarray(perm_means)
    p_value = (np.sum(perm_means_arr <= original_mean) + 1) / n_perms
    return float(p_value)


def _format_pvals_per_tissue(pairs_pv_df: pd.DataFrame) -> str:
    if pairs_pv_df is None or pairs_pv_df.empty:
        return ""
    parts = []
    for _, r in pairs_pv_df.iterrows():
        parts.append(f"Tissue {r['Tissue_ID']}: p-value={float(r['p_value_dist']):.3g}")
    return " ; ".join(parts)


def _format_val(val: float) -> str:
    return f"{val:.3e}" if abs(val) < 0.01 else f"{val:.3f}"


# =========================
# Core runner for one direction
# =========================

def _run_one_direction(
    *,
    results_dir: Path,
    out_dir: Path,
    tissues: List[str],
    cell_types: List[str],
    primary_label_fn,   # callable(cell_type)->str
    neighbor_label_fn,  # callable(cell_type)->str
    hr_positive_tissues: List[str],
    tissue_colors: Dict[str, Any],
    params: ComparisonAcrossPlatformsRunParams,
    rng: np.random.Generator,
    direction_tag: str,
    mode: str,  # "not_tumor_to_tumor" or "tumor_to_not_tumor"
) -> pd.DataFrame:
    adjust_text = _try_import_adjust_text()
    rows: List[Dict[str, Any]] = []

    plots_dir = out_dir / "plots"
    plots_dir.mkdir(parents=True, exist_ok=True)

    total_pairs = len(cell_types)
    if params.verbose:
        print(f"\n=== Running {mode} : {total_pairs} pairs ===")

    for idx, cell_type in enumerate(cell_types, start=1):
        # raw labels (for printing)
        if mode == "not_tumor_to_tumor":
            primary_ct_raw = str(cell_type)
            neighbor_ct_raw = str(params.tumor_cell_type)
        elif mode == "tumor_to_not_tumor":
            primary_ct_raw = str(params.tumor_cell_type)
            neighbor_ct_raw = str(cell_type)
        else:
            raise ValueError(f"Unknown mode: {mode}")

        # normalized for filename matching
        primary_ct = _normalize_cell_type_for_filename(primary_ct_raw)
        neighbor_ct = _normalize_cell_type_for_filename(neighbor_ct_raw)

        pair_label = f"{mode} ({idx}/{total_pairs})  {primary_ct_raw} vs {neighbor_ct_raw}"
        if params.verbose:
            print(f"\n--- {pair_label} ---")
            print("[stage] loading p-values")

        pval_df = _load_pval_matrix_from_final_genes(
            results_dir=results_dir,
            tissues=tissues,
            primary_cell_type=primary_ct,
            neighbor_cell_type=neighbor_ct,
            file_prefix_template=params.final_genes_prefix_template,
            gene_column=params.gene_column,
            pval_column=params.pval_column,
        )

        if pval_df is None or pval_df.empty:
            if params.verbose:
                print("[skip] no matching files / empty pval matrix")
            continue

        if params.verbose:
            print(f"[stage] pval matrix loaded: genes={pval_df.shape[0]}, tissues_found={pval_df.shape[1]}")

        min_max_df = pd.DataFrame({"min": pval_df.min(), "max": pval_df.max()})

        # replace zeros with 0.8 * global min (avoid -log10(0))
        min_val = float(pval_df.min().min())
        replacement_value = min_val * 0.8 if min_val > 0 else 1e-300
        pval_df_replaced = pval_df.replace(0, replacement_value)

        # -log10
        pval_df_n = -np.log10(pval_df_replaced)

        # standardize tissues x genes
        if params.verbose:
            print("[stage] scaling + PCA (real)")

        pval_df_scaled = StandardScaler().fit_transform(pval_df_n.T)

        pca = PCA(n_components=2)
        pca_results = pca.fit_transform(pval_df_scaled)
        real_var = pca.explained_variance_ratio_

        # permutation p-values for explained variance
        if params.verbose:
            print(f"[stage] PC permutation p-values (n={params.n_perm_pca})")

        p_values_PCs = _pc_permutation_pvalues(
            scaled_matrix=pval_df_scaled,
            real_var=real_var,
            n_perm=params.n_perm_pca,
            rng=rng,
            verbose=params.verbose,
            progress_every=params.perm_progress_every,
            label=pair_label,
        )

        if not bool(np.all(p_values_PCs <= params.pc_pvalue_alpha)):
            if params.verbose:
                print(f"[skip] PC p-values not significant: PC1={p_values_PCs[0]:.3g}, PC2={p_values_PCs[1]:.3g}")
            continue

        # PCA df
        pca_df = pd.DataFrame(pca_results, columns=["PC1", "PC2"])
        pca_df["Tissue"] = list(pval_df_n.columns)

        pca_df["ReceptorStatus"] = pca_df["Tissue"].apply(
            lambda x: "HR+/HER2-" if str(x) in set(hr_positive_tissues) else "Other"
        )

        # platform distance p-values
        if params.verbose:
            print(f"[stage] platform distance permutation (n={params.n_perm_dist})")

        pairs_pv_df = _platform_pair_distance_pvalues(
            pval_df_n=pval_df_n,
            pca_df=pca_df,
            n_perms=params.n_perm_dist,
            rng=rng,
            verbose=params.verbose,
            progress_every=params.perm_progress_every,
            label=pair_label,
        )

        # filter significant tissues and paired ids
        significant_tissues = min_max_df[min_max_df["min"] <= params.min_pv_for_tissue_filter].index
        numeric_tissues = [t for t in significant_tissues if str(t).isnumeric()]
        paired_tissues = [t for t in numeric_tissues if f"{params.platform_prefix}{t}" in significant_tissues]
        pairs_pv_df = pairs_pv_df[pairs_pv_df["Tissue_ID"].isin(paired_tissues)].copy()

        # receptor-status proximity p-value
        if params.verbose:
            print(f"[stage] receptor-status proximity permutation (n={params.n_perm_dist})")

        receptor_pv = _exseq_hr_proximity_pvalue(
            pval_df_n=pval_df_n,
            pca_df=pca_df,
            min_max_df=min_max_df,
            n_perms=params.n_perm_dist,
            platform_prefix=params.platform_prefix,
            rng=rng,
            verbose=params.verbose,
            progress_every=params.perm_progress_every,
            label=pair_label,
        )

        # plotting
        if params.verbose:
            print("[stage] plotting")

        pca_df_filtered = pca_df[pca_df["Tissue"].isin(significant_tissues)].copy()

        marker_styles = {"HR+/HER2-": "*", "Other": "o"}
        plt.figure(figsize=params.figsize)

        for _, r in pca_df_filtered.iterrows():
            tissue = str(r["Tissue"])
            receptor_status = str(r["ReceptorStatus"])
            color = tissue_colors.get(tissue, (0.2, 0.2, 0.2))
            marker = marker_styles.get(receptor_status, "o")
            size = 300 if marker == "*" else 100
            plt.scatter(
                float(r["PC1"]),
                float(r["PC2"]),
                color=color,
                marker=marker,
                s=size,
                edgecolor="black",
                alpha=0.7,
            )

        texts = []
        for _, r in pca_df_filtered.iterrows():
            display_name = str(r["Tissue"]).replace(params.platform_prefix, "M_")
            t = plt.text(float(r["PC1"]), float(r["PC2"]), display_name, fontsize=10, ha="center", va="center")
            texts.append(t)

        if adjust_text is not None:
            adjust_text(
                texts,
                expand=(2.7, 2.7),
                arrowprops=dict(arrowstyle="-", lw=0.5),
                force_text=(3.5, 3.5),
                force_static=(3.5, 3.5),
                force_explode=(3.5, 3.5),
                avoid_self=True,
                prevent_crossings=True,
                expand_axes=True,
                time_lim=5,
            )
        else:
            warnings.warn("adjustText is not installed; labels may overlap. Install with: pip install adjustText")

        legend_handles_color = []
        for tissue in pca_df_filtered["Tissue"]:
            tissue = str(tissue)
            color = tissue_colors.get(tissue, (0.2, 0.2, 0.2))
            min_v = _format_val(float(min_max_df.loc[tissue, "min"])) if tissue in min_max_df.index else ""
            max_v = _format_val(float(min_max_df.loc[tissue, "max"])) if tissue in min_max_df.index else ""
            display_name = tissue.replace(params.platform_prefix, "M_")
            label_txt = f"{display_name} ({min_v}, {max_v})"
            legend_handles_color.append(mpatches.Patch(color=color, label=label_txt))

        shape_handles = [
            plt.Line2D([0], [0], marker="*", color="w", label="HR+/HER2-", markerfacecolor="gray",
                       markersize=12, markeredgecolor="black"),
            plt.Line2D([0], [0], marker="o", color="w", label="Other receptor status", markerfacecolor="gray",
                       markersize=10, markeredgecolor="black"),
        ]

        legend1 = plt.legend(
            handles=legend_handles_color,
            title="Tissue (min, max)",
            bbox_to_anchor=(1.02, 1.0),
            loc="upper left",
            borderaxespad=0.5,
            fontsize=8,
            title_fontsize=10,
            handlelength=1.2,
            handletextpad=0.6,
            labelspacing=0.6,
            borderpad=0.7,
        )
        plt.gca().add_artist(legend1)

        plt.legend(
            handles=shape_handles,
            title="Biopsy receptor status",
            bbox_to_anchor=(1.02, 0.3),
            loc="upper left",
            borderaxespad=0.5,
            fontsize=8,
            title_fontsize=10,
            handlelength=1.2,
            handletextpad=0.6,
            labelspacing=0.6,
            borderpad=0.7,
        )

        primary_label = primary_label_fn(cell_type)
        neighbor_label = neighbor_label_fn(cell_type)

        plt.title(f"PCA of DeSeq P-Values: {primary_label} proximal to {neighbor_label}")
        plt.xlabel(f"PC1 (p-value = {float(p_values_PCs[0]):.3g})")
        plt.ylabel(f"PC2 (p-value = {float(p_values_PCs[1]):.3g})")
        plt.grid(True)
        plt.tight_layout()

        fig_name = f"pca_{direction_tag}_primary_{primary_label}_neighbor_{neighbor_label}.tif"
        fig_path = plots_dir / fig_name.replace(" ", "_")
        plt.savefig(fig_path, format="tiff", dpi=300, bbox_inches="tight")
        plt.close()

        if params.verbose:
            print(f"[done] saved plot: {fig_path}")

        rows.append(
            {
                "Primary cell type": primary_label,
                "Neighbor cell type": neighbor_label,
                "ExSeq-MERFISH similarity p-value": _format_pvals_per_tissue(pairs_pv_df),
                "Receptor status similarity p-value": (np.nan if receptor_pv is None else float(receptor_pv)),
                # explained variance (NEW)
                "PC1 explained variance": float(real_var[0]),
                "PC2 explained variance": float(real_var[1]),

                "PC1 p-value": float(p_values_PCs[0]),
                "PC2 p-value": float(p_values_PCs[1]),
                "Plot path": str(fig_path),
            }
        )

    return pd.DataFrame(rows)


# =========================
# Public API
# =========================

def run_comparison_across_platforms_and_patients(
    *,
    results_dir: Path,
    out_dir: Path,
    params: ComparisonAcrossPlatformsRunParams,
    rng_seed: Optional[int] = None,
) -> Dict[str, Any]:
    results_dir = Path(results_dir)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(rng_seed)

    hr_positive_tissues = params.hr_positive_tissues or []
    tissue_colors = _tissue_color_map(params.tissues, params.platform_prefix)

    if params.verbose:
        print("\n==============================")
        print("Platform comparison started")
        print(f"results_dir: {results_dir}")
        print(f"out_dir:     {out_dir}")
        print(f"tissues:     {len(params.tissues)}")
        print(f"cell types:  {len(params.non_tumor_cell_types)} (+ tumor={params.tumor_cell_type})")
        print(f"n_perm_pca:  {params.n_perm_pca}")
        print(f"n_perm_dist: {params.n_perm_dist}")
        print("==============================\n")

    df_not_tumor_to_tumor = _run_one_direction(
        results_dir=results_dir,
        out_dir=out_dir / "not_tumor_to_tumor",
        tissues=params.tissues,
        cell_types=params.non_tumor_cell_types,
        primary_label_fn=lambda ct: str(ct),
        neighbor_label_fn=lambda ct: str(params.tumor_cell_type),
        hr_positive_tissues=hr_positive_tissues,
        tissue_colors=tissue_colors,
        params=params,
        rng=rng,
        direction_tag="not_tumor_to_tumor",
        mode="not_tumor_to_tumor",
    )

    df_tumor_to_not_tumor = _run_one_direction(
        results_dir=results_dir,
        out_dir=out_dir / "tumor_to_not_tumor",
        tissues=params.tissues,
        cell_types=params.non_tumor_cell_types,
        primary_label_fn=lambda ct: str(params.tumor_cell_type),
        neighbor_label_fn=lambda ct: str(ct),
        hr_positive_tissues=hr_positive_tissues,
        tissue_colors=tissue_colors,
        params=params,
        rng=rng,
        direction_tag="tumor_to_not_tumor",
        mode="tumor_to_not_tumor",
    )

    combined = pd.concat(
        [
            df_not_tumor_to_tumor.assign(Direction="NotTumor→Tumor"),
            df_tumor_to_not_tumor.assign(Direction="Tumor→NotTumor"),
        ],
        ignore_index=True,
    )

    summary_csv = out_dir / "comparison_across_platforms_and_patients_summary.csv"
    combined.to_csv(summary_csv, index=False)

    if params.verbose:
        print("\n==============================")
        print("Platform comparison finished")
        print(f"Summary CSV: {summary_csv}")
        print(f"Rows: {combined.shape[0]}")
        print("==============================\n")

    return {
        "summary_csv": str(summary_csv),
        "n_rows": int(combined.shape[0]),
        "out_dir": str(out_dir),
    }
