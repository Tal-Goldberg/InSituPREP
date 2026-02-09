"""
expression_as_function_of_distance.py

Implements the "Expression as a Function of Distance" analysis.

GOAL (conceptual):
For each tissue and each ordered pair of cell types (primary_cell_type -> neighbor_cell_type),
we compute the minimum Euclidean distance from each primary cell to the nearest neighbor cell.
We then perform linear regression to identify genes whose expression correlates with distance.

We quantify this using:
    - Linear regression slope, R², and p-values
    - Permutation-based statistical significance testing
    - Global FDR correction
    - Gaussian-smoothed R² filtering

This module also supports:
    1) Low Expression Gene (LEG) filtering
    2) Gaussian-smoothed regression filtering
    3) Heatmap visualization of results
    4) Multi-panel expression-distance plots
"""

from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Any, List, Tuple, Mapping
import warnings
import json
import pandas as pd
import numpy as np
import seaborn as sns
from scipy.stats import linregress, norm
from scipy.ndimage import gaussian_filter1d
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests
from sklearn.preprocessing import StandardScaler
from dataclasses import asdict, is_dataclass
from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# =========================
# Parameters / configuration
# =========================

@dataclass
class ExpressionDistanceRunParams:
        sigma_param: float = 5.0
        r2_threshold: float = 0.3
        maximum_distance: float = 145.0
        leg_threshold: float = 10.0
        cell_count_threshold: int = 20
        num_iterations: int = 100
        fdr_thresh: float = 0.05

        # filters
        tissue_ids: Optional[List[str]] = None
        primary_cell_types: Optional[List[str]] = None
        neighbor_cell_types: Optional[List[str]] = None


# =========================
# IO helpers
# =========================

def read_genes_list(genes_names_path: Path) -> List[str]:
        """Read gene names from a file."""
        genes: List[str] = []
        with open(genes_names_path, "r", encoding="utf-8") as f:
                for line in f:
                        g = line.strip()
                        if g:
                                genes.append(g)
        return genes


def load_summary_table(summary_table_path: Path) -> pd.DataFrame:
        """Load summary table with cell data, coordinates, and gene expressions."""
        return pd.read_csv(summary_table_path, dtype={"Var1": str, "tissue": str})


def save_run_metadata(
        output_path: Path,
        summary_table_path: Path,
        genes_names_path: Path,
        params: ExpressionDistanceRunParams,
) -> Path:
        """
        Save a JSON file describing the inputs, output path and chosen parameters.
        Returns the path to the written file.
        """
        meta = {
                "timestamp_utc": datetime.utcnow().isoformat() + "Z",
                "summary_table_path": str(summary_table_path),
                "genes_names_path": str(genes_names_path),
                "output_path": str(output_path),
                "params": asdict(params) if is_dataclass(params) else dict(params),
        }

        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)
        out_file = output_path / "run_parameters.json"
        with open(out_file, "w", encoding="utf-8") as fh:
                json.dump(meta, fh, indent=2, sort_keys=True)
        return out_file


# =========================
# Distance computation helpers
# =========================

def create_min_dist_vector_from_coordinates(
        tissue_list: List[str],
        cell_type_1: str,
        cell_type_2: str,
        summary_table: pd.DataFrame,
        maximum_distance: float,
) -> pd.Series:
        """
        Compute minimum Euclidean distance from each cell in cell_type_1
        to the nearest cell in cell_type_2 across tissues.
        
        Uses X_space, Y_space, and optionally Z_space coordinates.
        Distances above maximum_distance are replaced with NaN.
        """
        min_dist_vector = pd.Series(dtype=float)

        for tissue in tissue_list:
                print(f"Processing tissue: {tissue}")
                
                df_tissue = summary_table[summary_table['tissue'] == tissue].copy()
                df_ct1 = df_tissue[df_tissue['cell_type'] == cell_type_1].copy()
                df_ct2 = df_tissue[df_tissue['cell_type'] == cell_type_2].copy()

                if df_ct1.empty or df_ct2.empty:
                        print(f"Skipping tissue {tissue}: one or both cell types are missing.")
                        continue

                # Get coordinates (2D or 3D)
                if 'Z_space' in df_tissue.columns:
                        coords_ct1 = df_ct1[['X_space', 'Y_space', 'Z_space']].values
                        coords_ct2 = df_ct2[['X_space', 'Y_space', 'Z_space']].values
                        print("Using 3D coordinates (X, Y, Z).")
                else:
                        coords_ct1 = df_ct1[['X_space', 'Y_space']].values
                        coords_ct2 = df_ct2[['X_space', 'Y_space']].values
                        print("Using 2D coordinates (X, Y).")
                
                tree_ct2 = cKDTree(coords_ct2)
                
                if cell_type_1 == cell_type_2:
                        # Exclude self-distance
                        distances, _ = tree_ct2.query(coords_ct1, k=2)
                        distances = distances[:, 1]
                else:
                        distances, _ = tree_ct2.query(coords_ct1, k=1)

                # Apply threshold
                distances = np.where(distances > maximum_distance, np.nan, distances)
                
                min_dist_vector = pd.concat([
                        min_dist_vector,
                        pd.Series(distances, index=df_ct1["Var1"])
                ])

                if not min_dist_vector.index.is_unique:
                        print("Warning: min_dist_vector index is not unique.")

                print(f"Initial cells: {len(df_ct1)}, within max distance: {np.sum(~np.isnan(distances))}")
                if np.all(np.isnan(distances)):
                        print("All distances are NaN after thresholding.")
                else:
                        print(f"Min: {np.nanmin(distances):.2f}, Max: {np.nanmax(distances):.2f}")


        return min_dist_vector


# =========================
# Filtering helpers
# =========================

def analyze_expression_by_distance_filter_leg(
        tissue_list: List[str],
        cell_type_1: str,
        cell_type_2: str,
        summary_table: pd.DataFrame,
        min_dist: pd.Series,
        exp_mat: pd.DataFrame,
        leg_threshold: float,
        cell_count_threshold: int,
        output_path: Path,
) -> pd.DataFrame:
        """
        Filter expression matrix by Low Expression Genes (LEG) threshold
        and minimum cell count per gene and cell type.
        These genes are set to NaN in the summary table for the given cell type and tissue.
        """
        summary_table_leg = summary_table.copy()
        min_dist.index = min_dist.index.astype(str)

        for tissue in tissue_list:
                mask = (summary_table_leg.tissue == tissue) & (summary_table_leg.cell_type == cell_type_1)
                tissue_df = summary_table_leg[mask].copy()
                tissue_df["Var1"] = tissue_df["Var1"].astype(str)

                min_dist_drop = min_dist.dropna()
                valid_cells = tissue_df["Var1"].isin(min_dist_drop.index)
                
                if valid_cells.sum() == 0:
                        continue

                filtered_df = tissue_df[valid_cells].copy()
                expr = filtered_df.loc[:, exp_mat.columns]

                valid_counts = expr.notna().sum(axis=0)
                max_expr = pd.Series(np.percentile(expr, 98, axis=0), index=expr.columns)
                low_expr_or_few_cells = (valid_counts <= cell_count_threshold) | (max_expr <= leg_threshold)

                low_expression_genes = low_expr_or_few_cells[low_expr_or_few_cells].index
                high_expression_genes = low_expr_or_few_cells[~low_expr_or_few_cells].index

                print(f"Tissue {tissue}:")
                print(f"  Genes filtered: {len(low_expression_genes)}, retained: {len(high_expression_genes)}")

                # Save LEG filter plot
                output_folder = Path(output_path) / "LEGfilter"
                output_folder.mkdir(parents=True, exist_ok=True)
                output_file = output_folder / f"max_expr_tissue_{tissue}_pct_{cell_type_1}_nct_{cell_type_2}.png"
                
                plt.figure(figsize=(4, 4))
                sns.histplot(max_expr.dropna(), bins=50, kde=False, color="skyblue")
                plt.axvline(leg_threshold, color='red', linestyle='--', label=f'Threshold = {leg_threshold}')
                plt.title(f"Max Expression — Tissue {tissue}, {cell_type_1} vs {cell_type_2}, Retained genes: {len(high_expression_genes)}")
                plt.xlabel("Max Expression")
                plt.ylabel("Number of Genes")
                plt.legend()
                plt.savefig(output_file, bbox_inches='tight')
                plt.close()

                # Replace with NaN
                summary_table_leg.loc[mask, low_expression_genes] = np.nan

        return summary_table_leg


# =========================
# Linear Regression Analysis
# =========================

def analyze_expression_by_distance_linear_regression(
        tissue_list: List[str],
        cell_type_1: str,
        cell_type_2: str,
        summary_table: pd.DataFrame,
        num_iterations: int,
        min_dist: pd.Series,
        exp_mat: pd.DataFrame,
        output_path: Path,
) -> pd.DataFrame:
        """
        Perform linear regression of gene expression on distance,
        with permutation-based significance testing.
        """
        all_results = []

        for tissue in tissue_list:
                print(f"Tissue: {tissue}")
                
                summary_table_tissue = summary_table[summary_table.tissue == tissue].copy()
                summary_table_ct1 = summary_table_tissue[summary_table_tissue.cell_type == cell_type_1].copy()
                
                summary_table_ct1["Var1"] = summary_table_ct1["Var1"].astype(str)
                min_dist.index = min_dist.index.astype(str)
                
                filtered_df = summary_table_ct1[summary_table_ct1["Var1"].isin(min_dist.index)].copy()
                filtered_df["min_dist"] = filtered_df["Var1"].map(min_dist)
                
                expression_with_dist = filtered_df.loc[:, ["Var1"] + list(exp_mat.columns) + ["min_dist"]]
                expression_with_dist_filt = expression_with_dist.dropna(subset=["min_dist"])

                if expression_with_dist_filt.empty:
                        print(f"No valid data for tissue {tissue}, {cell_type_1} vs {cell_type_2}.")
                        continue

                results = []

                # Linear regression for each gene
                for gene in expression_with_dist_filt.columns:
                        if gene in ["Var1", "min_dist"]:
                                continue
                        
                        y = expression_with_dist_filt[gene]
                        x = expression_with_dist_filt["min_dist"]
                        valid = (~y.isna()) & (~x.isna())
                        
                        if valid.sum() < 3:
                                continue

                        slope, intercept, r_value, p_value, std_err = linregress(x[valid], y[valid])

                        results.append({
                                "tissue": tissue,
                                "primary_cell_type": cell_type_1,
                                "neighbor_cell_type": cell_type_2,
                                "gene": gene,
                                "slope": slope,
                                "r2_orig": r_value ** 2,
                                "p_value_orig": p_value,
                                "cells_num": x[valid].shape[0],
                                "p98_expression": np.percentile(y[valid], 98),
                                "p98_distance": np.percentile(x[valid], 98)
                        })

                regression_results = pd.DataFrame(results)

                # Permutation test
                if not regression_results.empty and num_iterations > 0:
                        print(f"Running {num_iterations} permutations...")
                        genes = expression_with_dist_filt.columns
                        perm_pval_store = {gene: [] for gene in genes}
                        perm_r2_store = {gene: [] for gene in genes}

                        for iter in range(num_iterations):
                                # print every 50 iterations
                                if (iter + 1) % 50 == 0:
                                        print(f"  Permutation {iter + 1}/{num_iterations}")
                                
                                summary_perm = summary_table_ct1.copy()
                                summary_perm["Var1"] = np.random.permutation(summary_perm["Var1"].values)
                                filtered_perm = summary_perm[summary_perm["Var1"].isin(min_dist.index)].copy()
                                filtered_perm["min_dist"] = filtered_perm["Var1"].map(min_dist)
                                
                                expr_perm = filtered_perm.loc[:, ["Var1"] + list(exp_mat.columns) + ["min_dist"]]
                                expr_perm = expr_perm.dropna(subset=["min_dist"])

                                for gene in expr_perm.columns:
                                        if gene in ["Var1", "min_dist"]:
                                                continue
                                        
                                        y = expr_perm[gene]
                                        x = expr_perm["min_dist"]
                                        valid = (~y.isna()) & (~x.isna())
                                        
                                        if valid.sum() < 3:
                                                continue
                                        
                                        slope, intercept, r_value, p_value, std_err = linregress(x[valid], y[valid])
                                        perm_pval_store[gene].append(p_value)
                                        perm_r2_store[gene].append(r_value ** 2)

                        # Compute permutation-based p-values
                        perm_pvals = []
                        for gene in regression_results["gene"]:
                                real_pv = regression_results.loc[regression_results["gene"] == gene, "p_value_orig"].values[0]
                                rand_pvs = np.array(perm_pval_store[gene])
                                real_z_pv = (real_pv - rand_pvs.mean()) / (rand_pvs.std(ddof=1) + 1e-10)
                                pv_perm_pv = norm.cdf(real_z_pv)
                                real_r2 = regression_results.loc[regression_results["gene"] == gene, "r2_orig"].values[0]
                                rand_r2 = np.array(perm_r2_store[gene])
                                real_z_r2 = (real_r2 - rand_r2.mean()) / (rand_r2.std(ddof=1) + 1e-10)
                                if np.isfinite(real_z_r2):
                                        pv_perm_r2 = max(1 - norm.cdf(real_z_r2), 1e-16)
                                else:
                                        pv_perm_r2 = np.nan
                                perm_pvals.append(pv_perm_r2)

                        regression_results["p_value_perm"] = perm_pvals
                        regression_results["r2_perm_mean"] = regression_results["gene"].map(lambda g: np.mean(perm_r2_store[g]))


                all_results.append(regression_results)

        if len(all_results) > 0:
                combined_df = pd.concat(all_results, ignore_index=True)
        else:
                return pd.DataFrame()
        
        #save
        output_folder = Path(output_path) / "LinearRegression"
        output_folder.mkdir(parents=True, exist_ok=True)
        output_file = output_folder / f"linear_regression_results_{cell_type_1}_{cell_type_2}.csv"
        combined_df.to_csv(output_file, index=False)

        return combined_df


# =========================
# FDR Filtering
# =========================

def apply_global_fdr_and_filter(
        results_by_pairs: Dict[str, pd.DataFrame],
        fdr_thresh: float,
        output_path: Path,
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
        """
        Apply global FDR correction across all genes, cell type pairs, and tissues.
        
        Parameters:
        results_by_pairs (dict): {celltype_pair: DataFrame} from linear regression
        FDR_thresh (float): FDR significance threshold
        output_path (str): Path to output folder

        Returns:
        filtered_results (dict): Updated dict with NaN-filtered DataFrames
        summary_df (DataFrame): Initial and significant gene counts per pair of cell types
        """
        all_pvals = []

        for pair, df in results_by_pairs.items():
                if df.empty or "gene" not in df.columns or "tissue" not in df.columns or "p_value_orig" not in df.columns:
                        print(f"Skipping empty or missing column DataFrame for pair {pair}")
                        continue
                
                temp = df[["gene", "tissue", "p_value_orig"]].copy()
                temp["pair_cell_type"] = pair
                all_pvals.append(temp)

        if not all_pvals:
                return {}, pd.DataFrame()

        all_pvals_df = pd.concat(all_pvals, ignore_index=True)
        _, fdr_vals, _, _ = multipletests(all_pvals_df["p_value_orig"], method='fdr_bh')
        all_pvals_df["p_value_fdr_global"] = fdr_vals

        filtered_results = {}
        summary_records = []

        for pair, df in results_by_pairs.items():
                if df.empty:
                        continue

                df_pair = df.copy()
                df_pair = df_pair.merge(
                        all_pvals_df[all_pvals_df["pair_cell_type"] == pair][["gene", "tissue", "p_value_fdr_global"]],
                        on=["gene", "tissue"],
                        how="left"
                )

                tissue_counts_before = df_pair.groupby("tissue")["gene"].nunique()
                df_pair.loc[df_pair["p_value_fdr_global"] > fdr_thresh, :] = np.nan
                df_cleaned = df_pair.dropna(how="all")
                tissue_counts_after = df_cleaned.groupby("tissue")["gene"].nunique()

                filtered_results[pair] = df_cleaned

                for tissue in tissue_counts_before.index:
                        summary_records.append({
                                "pair_cell_type": pair,
                                "tissue": tissue,
                                "initial_genes": tissue_counts_before.get(tissue, 0),
                                "significant_genes": tissue_counts_after.get(tissue, 0)
                        })

        summary_df = pd.DataFrame(summary_records)

        # Save results
        output_folder = Path(output_path) / "FDRfilter"
        output_folder.mkdir(parents=True, exist_ok=True)
        summary_df.to_csv(output_folder / "FDR_filter_summary.csv", index=False)

        for pair, df in filtered_results.items():
                df.to_csv(output_folder / f"{pair}_FDR_filtered.csv", index=False)

        # Plot
        if not summary_df.empty:
                plt.figure(figsize=(14, 6))
                plot_df = summary_df.set_index(["pair_cell_type", "tissue"])
                plot_df[["initial_genes", "significant_genes"]].plot(kind="bar", alpha=0.8)
                plt.title("Initial vs. Significant Genes per Cell Type Pair and Tissue")
                plt.ylabel("Number of Genes")
                plt.xlabel("Cell Type Pair & Tissue")
                plt.xticks(rotation=90)
                plt.savefig(output_folder / "FDR_filter_summary.png", bbox_inches='tight')
                plt.close()

        return filtered_results, summary_df


# =========================
# Gaussian R² Filtering
# =========================

def apply_gaussian_r2_filter(
        results_by_pairs: Dict[str, pd.DataFrame],
        summary_table_by_pairs: Dict[str, pd.DataFrame],
        min_dist_by_pair: Dict[str, pd.Series],
        r2_threshold: float,
        sigma: float,
        exp_mat: pd.DataFrame,
        output_path: Path,
) -> Tuple[Dict[str, pd.DataFrame], pd.DataFrame]:
        """
        Apply Gaussian-smoothed R² filtering to results.
        """
        filtered_results = {}
        summary_records = []
        r2_by_pair_tissue_gene = {}

        for pair, df in results_by_pairs.items():
                if df.empty or "gene" not in df.columns or "tissue" not in df.columns:
                        print(f"Skipping empty or missing column DataFrame for pair {pair}")
                        continue

                if pair not in min_dist_by_pair or pair not in summary_table_by_pairs:
                        print(f"Missing min_dist or summary_table for pair {pair}, skipping Gaussian R² filter.")
                        continue

                min_dist = min_dist_by_pair[pair].copy()
                summary_table = summary_table_by_pairs[pair].copy()
                cell_type = df["primary_cell_type"].iloc[0]
                tissues = df["tissue"].dropna().unique()

                for tissue in tissues:
                        tissue = str(tissue)
                        subset = summary_table[(summary_table.tissue == tissue) & (summary_table.cell_type == cell_type)].copy()
                        subset["Var1"] = subset["Var1"].astype(str)
                        min_dist.index = min_dist.index.astype(str)
                        subset = subset[subset["Var1"].isin(min_dist.index)].copy()
                        subset["min_dist"] = subset["Var1"].map(min_dist)
                        gene_expr_df = subset.dropna(subset=["min_dist"])

                        failing_genes = []

                        for gene in exp_mat.columns:
                                gene_df = gene_expr_df[["min_dist", gene]].dropna().sort_values("min_dist")
                                if gene_df.shape[0] < 3:
                                        continue
                                
                                x = gene_df["min_dist"].values
                                y = gene_df[gene].values
                                y_smoothed = gaussian_filter1d(y, sigma=sigma)
                                slope, intercept, r_value, _, _ = linregress(x, y_smoothed)
                                r2 = r_value ** 2
                                r2_by_pair_tissue_gene[(pair, tissue, gene)] = r2
                                
                                if r2 < r2_threshold:
                                        failing_genes.append(gene)

                        tissue_mask = (df["tissue"].astype(str) == tissue)
                        before = df[tissue_mask]["gene"].nunique()
                        
                        if failing_genes:
                                df.loc[tissue_mask & df["gene"].isin(failing_genes), :] = np.nan
                        
                        after = df[tissue_mask].dropna(how="all")["gene"].nunique()

                        summary_records.append({
                                "pair_cell_type": pair,
                                "tissue": tissue,
                                "initial_genes": before,
                                "significant_genes": after,
                                "removed": before - after
                        })

                filtered_results[pair] = df.dropna(how="all")

        # Add r2_gaussian column
        for pair, df in filtered_results.items():
                if df.empty:
                        continue
                
                r2_values = []
                for idx, row in df.iterrows():
                        key = (pair, row["tissue"], row["gene"])
                        r2_val = r2_by_pair_tissue_gene.get(key, np.nan)
                        r2_values.append(r2_val)

                df = df.dropna(how="all").copy()
                df["r2_gaussian"] = r2_values
                filtered_results[pair] = df

        # Save results
        output_folder = Path(output_path) / "GaussianR2filter"
        output_folder.mkdir(parents=True, exist_ok=True)

        for pair, df in filtered_results.items():
                df.to_csv(output_folder / f"{pair}_gaussian_r2_filtered.csv", index=False)

        summary_df = pd.DataFrame(summary_records)

        if not summary_df.empty:
                summary_df.to_csv(output_folder / "Gaussian_R2_filter_summary.csv", index=False)

                plt.figure(figsize=(14, 6))
                plot_df = summary_df.copy()
                plot_df["label"] = plot_df["pair_cell_type"] + " | " + plot_df["tissue"].astype(str)
                plot_df = plot_df.sort_values("removed", ascending=False)
                plt.bar(plot_df["label"], plot_df["removed"], color="tomato")
                plt.xticks(rotation=90)
                plt.title(f"Number of Genes Removed (R² < {r2_threshold}) per Pair and Tissue")
                plt.ylabel("Number of Genes Removed")
                plt.xlabel("Cell Type Pair | Tissue")
                plt.savefig(output_folder / "Gaussian_R2_filter_summary.png", bbox_inches='tight')
                plt.close()

        return filtered_results, summary_df


# =========================
# Plotting
# =========================

def plot_multi_feature_heatmaps(
        df: pd.DataFrame,
        pair_label: str,
        output_path: Path,
) -> None:
        """
        Plot heatmaps of slope, -log10(FDR), and R² Gaussian.
       
        Parameters:
        -----------
        df : pd.DataFrame
                Filtered dataframe with columns ['gene', 'tissue', 'slope', 'p_value_fdr_global', 'r2_gaussian']
        pair_label : str
                Label for the plot title (e.g., "T_cells_Tumor")
        """
        if df.empty:
                print(f"[{pair_label}] Skipping plot — no data.")
                return

        required_cols = ["gene", "tissue", "slope", "p_value_fdr_global", "r2_gaussian"]
        if not all(col in df.columns for col in required_cols):
                print(f"[{pair_label}] Skipping plot — missing columns.")
                return

        df_pair = df.copy()
        melted = pd.melt(
                df_pair,
                id_vars=["gene", "tissue"],
                value_vars=["slope", "p_value_fdr_global", "r2_gaussian"],
                var_name="feature",
                value_name="value"
        )

        melted.loc[melted["feature"] == "p_value_fdr_global", "value"] = (
                -np.log10(melted.loc[melted["feature"] == "p_value_fdr_global", "value"].clip(lower=1e-16))
        )

        melted["col"] = melted["tissue"].astype(str) + "_" + melted["feature"]
        heatmap_data = melted.pivot(index="gene", columns="col", values="value")

        features = ["slope", "p_value_fdr_global", "r2_gaussian"]
        col_groups = {f: [c for c in heatmap_data.columns if c.endswith(f)] for f in features}

        fig_height = max(6, 0.4 * len(heatmap_data))
        fig, axes = plt.subplots(1, 3, figsize=(20, fig_height), sharey=True, gridspec_kw={'wspace': 0.1})

        cmaps = {"slope": "vlag", "p_value_fdr_global": "Greens", "r2_gaussian": "Blues"}
        cbar_labels = {"slope": "Slope", "p_value_fdr_global": "-log10(FDR)", "r2_gaussian": "R² Gaussian"}

        for ax, feature in zip(axes, features):
                cols = col_groups[feature]

                if not cols:
                        ax.axis('off')
                        ax.set_title(f"{feature} (no data)")
                        continue

                feature_data = heatmap_data[cols].dropna(how='all')

                if feature_data.empty:
                        ax.axis('off')
                        ax.set_title(f"{feature} (no valid values)")
                        continue

                sns.heatmap(
                        feature_data,
                        ax=ax,
                        cmap=cmaps[feature],
                        center=0 if feature == "slope" else None,
                        linewidths=0.5,
                        linecolor='white',
                        cbar_kws={"label": cbar_labels[feature], "location": "bottom", "pad": 0.05}
                )
                ax.set_title(cbar_labels[feature], fontsize=14, pad=10)
                ax.set_xlabel("Tissue", fontsize=12)
                if ax == axes[0]:
                        ax.set_ylabel("Gene", fontsize=12)

        plt.suptitle(f"Gene-Level Regression Features by Tissue — {pair_label}", fontsize=18, y=0.98)
        
        output_folder = Path(output_path) / "Heatmaps"
        output_folder.mkdir(parents=True, exist_ok=True)
        plt.savefig(output_folder / f"combined_heatmap_{pair_label}.png", bbox_inches='tight')
        plt.close()


def plot_expression_distance_panels(
        pair_name: str,
        tissue: str,
        gene: str,
        filtered_results_dict: Dict[str, pd.DataFrame],
        expression_dict: Dict[str, pd.DataFrame],
        min_dist_dict: Dict[str, pd.Series],
        sigma: float,
        output_path: Path,
) -> None:
        """
        Plot expression vs. distance panels for a given pair, gene and tissue.
        Parameters:
        pair_name (str): Cell type pair name (e.g., "T_cells_Tumor")
        tissue (str): Tissue ID to plot
        gene (str): Gene name to plot
        filtered_results_dict (dict): Dictionary of filtered results DataFrames by pair name
        expression_dict (dict): Dictionary of expression DataFrames by pair name
        min_dist_dict (dict): Dictionary of minimum distance Series by pair name
        sigma (float): Smoothing parameter for Gaussian filter
        output_path (Path): Path to save the plot
        """
        df_results = filtered_results_dict[pair_name]
        df_expression = expression_dict[pair_name]
        min_dist = min_dist_dict[pair_name]
        cell_type = df_results["primary_cell_type"].iloc[0]
        cell_type2 = df_results["neighbor_cell_type"].iloc[0]

        df_expression = df_expression[(df_expression["tissue"].astype(str) == str(tissue)) & (df_expression["cell_type"] == cell_type)].copy()
        df_expression["Var1"] = df_expression["Var1"].astype(str)
        min_dist.index = min_dist.index.astype(str)
        df_expression = df_expression[df_expression["Var1"].isin(min_dist.index)]
        df_expression["min_dist"] = df_expression["Var1"].map(min_dist)

        df_gene = df_expression[["min_dist", gene]].dropna().sort_values("min_dist")

        if df_gene.shape[0] < 3:
                print(f"Not enough data for {gene} in tissue {tissue}.")
                return

        x = df_gene["min_dist"].values
        y = df_gene[gene].values
        y_perm = np.random.permutation(y)
        y_smooth = gaussian_filter1d(y, sigma=sigma)
        y_perm_smooth = gaussian_filter1d(y_perm, sigma=sigma)

        output_folder = Path(output_path) / "Plots"
        output_folder.mkdir(parents=True, exist_ok=True)

        fig1, axes = plt.subplots(2, 3, figsize=(15, 8))
        fdr_val = df_results[df_results['gene'] == gene]['p_value_fdr_global'].values
        fdr_str = f"{fdr_val[0]:.1e}" if len(fdr_val) > 0 else "N/A"
        fig1.suptitle(fr"$\it{{{gene}}}$ expression in {cell_type} vs {cell_type2} | FDR: {fdr_str}", fontsize=22)

        raw_ymax = np.max(np.concatenate([y, y_perm])) * 1.1

        for i, (yy, label, color) in enumerate(zip([y, y_perm], ["Raw", "Shuffled y"], ["red", "gray"])):
                ax = axes[i, 0]
                slope, intercept, r, p, _ = linregress(x, yy)
                ax.scatter(x, yy, alpha=0.6, s=35)
                ax.plot(x, intercept + slope * x, color=color)
                ax.set_ylim(-1, raw_ymax)
                ax.set_title(label, fontsize=22)
                ax.set_ylabel("Expression", fontsize=18)
                ax.set_xlabel("Distance [μm]", fontsize=18)
                ax.text(0.05, 0.95, f"r²={r**2:.2f}\nslope={slope:.3f}\np={p:.2e}",
                                transform=ax.transAxes, verticalalignment='top',
                                bbox=dict(boxstyle="round", facecolor='white', alpha=0.7), fontsize=14)
                ax.tick_params(axis='both', which='major', labelsize=18)

        for i, (yy, label, color) in enumerate(zip([y_smooth, y_perm_smooth], [f"Smoothed σ={sigma}", f"Shuffled σ={sigma}"], ["red", "gray"])):
                ax = axes[i, 1]
                slope, intercept, r, p, _ = linregress(x, yy)
                ax.scatter(x, yy, alpha=0.6, s=35)
                ax.plot(x, intercept + slope * x, color=color)
                ax.set_ylim(-1, raw_ymax)
                ax.set_title(label, fontsize=22)
                ax.set_xlabel("Distance [μm]", fontsize=18)
                ax.text(0.05, 0.95, f"r²={r**2:.2f}\nslope={slope:.3f}\np={p:.2e}",
                                transform=ax.transAxes, verticalalignment='top',
                                bbox=dict(boxstyle="round", facecolor='white', alpha=0.7), fontsize=14)
                ax.tick_params(axis='both', which='major', labelsize=18)

        zoom_ymax = np.max(np.concatenate([y_smooth, y_perm_smooth])) * 1.1
        for i, (yy, label, color) in enumerate(zip([y_smooth, y_perm_smooth], ["Zoomed", "Zoomed"], ["red", "gray"])):
                ax = axes[i, 2]
                slope, intercept, r, p, _ = linregress(x, yy)
                ax.scatter(x, yy, alpha=0.6, s=35)
                ax.plot(x, intercept + slope * x, color=color)
                ax.set_ylim(-1, zoom_ymax)
                ax.set_title(label, fontsize=22)
                ax.set_xlabel("Distance [μm]", fontsize=18)
                ax.tick_params(axis='both', which='major', labelsize=18)

        fig1.tight_layout(rect=[0, 0, 1, 0.95])
        fig1.savefig(output_folder / f"All_panels_{gene}_{cell_type}_{cell_type2}.png", bbox_inches='tight')
        plt.close(fig1)

        # Raw panel only
        fig2, ax2 = plt.subplots(figsize=(6, 4))
        slope, intercept, r, p, _ = linregress(x, y)
        ax2.scatter(x, y, alpha=0.6, s=35)
        ax2.plot(x, intercept + slope * x, color="red")
        ax2.set_ylim(-1, raw_ymax)
        ax2.set_xlabel("Distance [μm]", fontsize=20)
        ax2.set_ylabel(fr"$\it{{{gene}}}$ Expression", fontsize=20)
        ax2.tick_params(axis='both', which='major', labelsize=18)
        fig2.savefig(output_folder / f"Raw_panel_{gene}_{cell_type}_{cell_type2}.png", bbox_inches='tight')
        plt.close(fig2)


# =========================
# Main Orchestration
# =========================

def run_expression_as_function_of_distance(
        summary_table_path: Path,
        genes_names_path: Path,
        output_path: Path,
        params: ExpressionDistanceRunParams,
) -> Dict[str, Any]:
        """
        Main function to run the full expression-distance analysis pipeline.
        """
        output_path = Path(output_path)
        output_path.mkdir(parents=True, exist_ok=True)

        # Load data
        summary_table = load_summary_table(summary_table_path)
        genes_names = read_genes_list(genes_names_path)

        # Get tissues and cell types
        tissue_list = summary_table["tissue"].astype(str).unique().tolist()
        if params.tissue_ids is not None and len(params.tissue_ids) > 0:
                wanted = set(params.tissue_ids)
                tissue_list = [t for t in tissue_list if t in wanted]

        if params.primary_cell_types is not None and len(params.primary_cell_types) > 0:
            primary_cell_types_list = params.primary_cell_types
        else:
            primary_cell_types_list = summary_table["cell_type"].unique().tolist()

        if params.neighbor_cell_types is not None and len(params.neighbor_cell_types) > 0:
            neighbor_cell_types_list = params.neighbor_cell_types
        else:
            neighbor_cell_types_list = summary_table["cell_type"].unique().tolist()

        if params.maximum_distance == 0:
                params.maximum_distance = np.inf

        # Filter to common genes
        common_genes = [g for g in genes_names if g in summary_table.columns]
        missing_genes = set(genes_names) - set(common_genes)
        if missing_genes:
                print(f"Warning: {len(missing_genes)} genes from the provided list are missing in the summary table.")
        if not common_genes:
                raise ValueError("No common genes found between the provided gene list and the summary table.")
        exp_mat = summary_table.loc[:, common_genes]
        exp_mat.index = summary_table["Var1"].astype(str)

        print(f"[START] Expression-Distance | tissues={tissue_list} | genes={len(common_genes)}")
        print(f"[PARAMS] leg_threshold={params.leg_threshold} | r2_threshold={params.r2_threshold} | fdr_thresh={params.fdr_thresh}")

        # Save run metadata
        save_run_metadata(output_path, summary_table_path, genes_names_path, params)

        # --- Stage 1: Distance computation ---
        print("[STAGE 1/6] Computing distance vectors...")
        dict_min_dist_by_pair = {}
        for pct in primary_cell_types_list:
                for nct in neighbor_cell_types_list:
                        key = f"{pct}_{nct}"
                        print(f"  Computing distance vector for {key}...")
                        min_dist = create_min_dist_vector_from_coordinates(
                                tissue_list, pct, nct, summary_table, params.maximum_distance
                        )
                        dict_min_dist_by_pair[key] = min_dist

        # --- Stage 2: LEG filtering ---
        print("[STAGE 2/6] Applying LEG filter...")
        dict_summary_table_leg_by_pair = {}
        for pct in primary_cell_types_list:
                for nct in neighbor_cell_types_list:
                        key = f"{pct}_{nct}"
                        print(f"  Filtering LEG for {key}...")
                        summary_table_leg = analyze_expression_by_distance_filter_leg(
                                tissue_list, pct, nct, summary_table,
                                dict_min_dist_by_pair[key], exp_mat,
                                params.leg_threshold, params.cell_count_threshold, output_path
                        )
                        dict_summary_table_leg_by_pair[key] = summary_table_leg

        # --- Stage 3: Linear regression ---
        print("[STAGE 3/6] Running linear regression...")
        dict_linear_regression_by_pair = {}
        for pct in primary_cell_types_list:
                for nct in neighbor_cell_types_list:
                        key = f"{pct}_{nct}"
                        print(f"  Regressing {key}...")
                        linear_regression = analyze_expression_by_distance_linear_regression(
                                tissue_list, pct, nct,
                                dict_summary_table_leg_by_pair[key],
                                params.num_iterations,
                                dict_min_dist_by_pair[key],
                                exp_mat, output_path
                        )
                        dict_linear_regression_by_pair[key] = linear_regression

        # --- Stage 4: FDR filtering ---
        print("[STAGE 4/6] Applying global FDR correction...")
        dict_filtered_global_fdr_by_pair, summary_fdr = apply_global_fdr_and_filter(
                dict_linear_regression_by_pair, params.fdr_thresh, output_path
        )

        # --- Stage 5: Gaussian R² filtering ---
        print("[STAGE 5/6] Applying Gaussian R² filter...")
        dict_filtered_gaussian_r2_by_pair, summary_gaussian_r2 = apply_gaussian_r2_filter(
                dict_filtered_global_fdr_by_pair,
                dict_summary_table_leg_by_pair,
                dict_min_dist_by_pair,
                params.r2_threshold, params.sigma_param, exp_mat,
                output_path
        )

        # --- Stage 6: Final results and plots ---
        print("[STAGE 6/6] Generating plots and final results...")
        dfs = [df for df in dict_filtered_gaussian_r2_by_pair.values() if df is not None and not df.empty]
        if dfs:
            final_combined_results = pd.concat(dfs, ignore_index=True)
        else:
            final_combined_results = pd.DataFrame()
        final_combined_results.to_csv(output_path / "Final_combined_results.csv", index=False)

        for pct in primary_cell_types_list:
                for nct in neighbor_cell_types_list:
                        key = f"{pct}_{nct}"
                        df_heatmap = dict_filtered_gaussian_r2_by_pair.get(key)
                        if df_heatmap is None or df_heatmap.empty:
                                continue
                        print(f"  Plotting heatmaps for {key}...")
                        plot_multi_feature_heatmaps(df_heatmap, key, output_path)

        for _, row in final_combined_results.iterrows():
                pair_name = f"{row['primary_cell_type']}_{row['neighbor_cell_type']}"
                tissue = str(row['tissue'])
                gene = row['gene']
                
                try:
                        plot_expression_distance_panels(
                                pair_name=pair_name,
                                tissue=tissue,
                                gene=gene,
                                filtered_results_dict=dict_filtered_gaussian_r2_by_pair,
                                expression_dict=dict_summary_table_leg_by_pair,
                                min_dist_dict=dict_min_dist_by_pair,
                                sigma=params.sigma_param,
                                output_path=output_path
                        )
                except Exception as e:
                        print(f"Failed to plot {pair_name}, {tissue}, {gene}: {e}")

        print(f"[DONE] Outputs written to: {output_path}")

        return {
                "final_results_csv": str(output_path / "Final_combined_results.csv"),
                "n_results": len(final_combined_results),
                "output_directory": str(output_path),
        }
