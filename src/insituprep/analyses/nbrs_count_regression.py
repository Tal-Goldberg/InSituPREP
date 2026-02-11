"""
Neighbor Count Regression Analysis.

Analyzes the relationship between gene expression and the number of neighboring
cells of a specific type. Implements three regression methods:
1. Standard linear regression
2. Sampling-based regression (2a: slope averaging, 2b: line-on-line)
3. Weighted least squares regression (3.1: unnormalized, 3.2: normalized)
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Tuple

import numpy as np
import pandas as pd
import warnings
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap


@dataclass
class NbrsCountRegressionParams:
    dist_threshold: float = 15.0
    min_expression: int = 10
    variance_quantile: float = 0.20
    min_neighbor_bins: int = 4
    min_cells_per_bin: int = 10
    n_iterations: int = 1000
    sample_size: int = 10
    fdr_alpha: float = 0.05


def read_genes_list(genes_names_path: Path) -> List[str]:
    """Read gene names from a text file (one gene per line)."""
    genes = []
    with open(genes_names_path, "r", encoding="utf-8") as f:
        for line in f:
            g = line.strip()
            if g:
                genes.append(g)
    return genes


def load_labels(labels_path: Path) -> pd.DataFrame:
    """Load the summary table CSV. The Var1 column (cell IDs) is read as string
    to preserve IDs like '1.9' vs '1.90' vs '1.900' (which collapse as float)."""
    # Determine which column is Var1 so we can force it to string dtype
    cols = pd.read_csv(labels_path, nrows=0).columns.tolist()
    dtype = {}
    if "Var1" in cols:
        dtype["Var1"] = str
    elif "cell_id" in cols:
        dtype["cell_id"] = str
    
    df = pd.read_csv(labels_path, dtype=dtype)
    
    if "Var1" in df.columns:
        df = df.rename(columns={"Var1": "cell_id"})
    elif "cell_id" not in df.columns:
        raise ValueError(
            "Summary table must contain a 'Var1' or 'cell_id' column with cell identifiers."
        )
    return df


def load_distance_matrix(distance_dir: Path, tissue) -> pd.DataFrame:
    """Load pre-computed distance matrix for a specific tissue.
    Uses row/column labels exactly as stored in the CSV (cell IDs),
    without adding any tissue prefix.
    """
    fp = distance_dir / f"distance_matrix_{tissue}.csv"
    dist = pd.read_csv(fp, index_col=0, dtype={0: str})

    # Ensure IDs are strings and consistent (no whitespace)
    dist.index = dist.index.map(lambda s: str(s).strip())
    dist.columns = [str(c).strip() for c in dist.columns]

    return dist

def calculate_neighbor_counts(
    sub_primary: pd.DataFrame,
    sub_neighbor: pd.DataFrame,
    dist_threshold: float,
    distance_matrix: Optional[pd.DataFrame] = None,
) -> np.ndarray:
    """
    Calculate how many neighbor cells are within dist_threshold of each primary cell.
    If distance_matrix provided, uses it (index-based lookup).
    Otherwise calculates from coordinates (cKDTree).
    """
    if distance_matrix is not None:
        dist_sub = distance_matrix.loc[sub_primary.index, sub_neighbor.index]
        return (dist_sub.values <= dist_threshold).sum(axis=1)
    
    # Calculate from coordinates using cKDTree
    coord_cols = ["X_space", "Y_space"]
    if "Z_space" in sub_primary.columns and sub_primary["Z_space"].notna().any():
        coord_cols.append("Z_space")
    
    coords_primary = sub_primary[coord_cols].values
    coords_neighbor = sub_neighbor[coord_cols].values
    
    tree = cKDTree(coords_neighbor)
    counts = tree.query_ball_point(coords_primary, r=dist_threshold, return_length=True)
    
    return np.array(counts)


def get_cell_type_pairs(
    labels: pd.DataFrame,
    tissue,
    primary_cell_types: Optional[List[str]] = None,
    neighbor_cell_types: Optional[List[str]] = None,
) -> List[Tuple[str, str]]:
    """
    Get all cell type pairs to analyze (including self-pairs), dynamically from data.
    If specific types provided, filters to those. Otherwise uses all types in data.
    """
    labels_tissue = labels[labels["tissue"].astype(str) == str(tissue)]
    all_cell_types = labels_tissue["cell_type"].dropna().unique().tolist()
    
    if primary_cell_types is not None:
        primary_types = [ct for ct in primary_cell_types if ct in all_cell_types]
    else:
        primary_types = all_cell_types
    
    if neighbor_cell_types is not None:
        neighbor_types = [ct for ct in neighbor_cell_types if ct in all_cell_types]
    else:
        neighbor_types = all_cell_types
    
    return [(p, n) for p in primary_types for n in neighbor_types]


def get_tissues(labels: pd.DataFrame) -> List[str]:
    """Get all unique tissue IDs from the data (as strings)."""
    return [str(t) for t in labels["tissue"].dropna().unique().tolist()]


def get_gene_columns(labels: pd.DataFrame, genes_names_path: Path) -> List[str]:
    """
    Get gene columns from external file, filtered to those that exist in data.
    Replaces METADATA_COLS = 9 approach.
    """
    requested_genes = read_genes_list(genes_names_path)
    available_genes = [g for g in requested_genes if g in labels.columns]
    
    if len(available_genes) < len(requested_genes):
        missing = set(requested_genes) - set(available_genes)
        warnings.warn(f"{len(missing)} genes not found in data: {missing}")
    
    return available_genes


# Color scheme for heatmap (from original config.py)
HEATMAP_COLORS = ['#000080', '#0066CC', '#00CCFF', '#FFFF00', '#FFD700']
HEATMAP_CMAP = LinearSegmentedColormap.from_list('blue_yellow', HEATMAP_COLORS, N=256)


def save_results_csv(
    results_df: pd.DataFrame,
    out_dir: Path,
    method: str,
    target_cell: str,
    neighbor_cell: str,
    dist_threshold: float,
    fdr_alpha: float = 0.05,
) -> None:
    """
    Save results to CSV files with descriptive names.
    
    File naming pattern:
    - {method}_{target}_vs_{neighbor}_dist{threshold}_all.csv
    - {method}_{target}_vs_{neighbor}_dist{threshold}_significant.csv
    
    Parameters:
        results_df: DataFrame with all results
        out_dir: Output directory
        method: Analysis method (1, 2a, 2b, 3.1, 3.2)
        target_cell, neighbor_cell: Cell type names
        dist_threshold: Distance threshold used
        fdr_alpha: Significance threshold for filtering
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # Create descriptive base name
    # Example: "method1_B_vs_Epithelial_dist15"
    base_name = f"method{method}_{target_cell}_vs_{neighbor_cell}_dist{int(dist_threshold)}"
    
    # Drop plot-data columns for CSV (large arrays that break Excel).
    # Plots are saved separately as PNGs; the in-memory DataFrame is not modified.
    plot_cols = [c for c in results_df.columns if c.startswith("plot_")]
    csv_df = results_df.drop(columns=plot_cols, errors="ignore")
    
    # Save all results
    all_results_path = out_dir / f"{base_name}_all.csv"
    csv_df.to_csv(all_results_path, index=False)
    print(f"Saved {len(csv_df)} total results to {all_results_path}")
    
    # Save significant results only
    if "q_value" in csv_df.columns:
        sig_df = csv_df[csv_df["q_value"] < fdr_alpha]
        sig_path = out_dir / f"{base_name}_significant.csv"
        sig_df.to_csv(sig_path, index=False)
        print(f"Saved {len(sig_df)} significant results to {sig_path}")


def create_gene_plot(
    gene_name: str,
    neighbor_counts: np.ndarray,
    expression_values: np.ndarray,
    slope: float,
    intercept: float,
    r_value: float,
    tissue,
    target_cell: str,
    neighbor_cell: str,
    original_neighbor_counts: Optional[np.ndarray] = None,
    original_expression_values: Optional[np.ndarray] = None,
) -> plt.Figure:
    """
    Create the 4-panel gene plot (same as your GUI).
    
    Panels:
    - Top-left: Expression vs Neighbors scatter + regression line
    - Top-right: Neighbor count histogram
    - Bottom-left: Mean expression ± std vs neighbors
    - Bottom-right: Fit on raw data (if Method 2/3)
    
    Parameters:
        gene_name: Name of the gene
        neighbor_counts: Array of neighbor counts (after bin filtering)
        expression_values: Array of expression values (after bin filtering)
        slope, intercept, r_value: Regression results
        tissue: Tissue ID
        target_cell, neighbor_cell: Cell type names
        original_*: Raw data before filtering (for Method 2/3)
        
    Returns:
        matplotlib Figure object
    """
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle(
        f"Gene Analysis: {gene_name} (Tissue {tissue})\n{target_cell} vs. {neighbor_cell}",
        fontsize=12
    )
    
    # --- Panel 1: Expression vs Neighbors (top-left) ---
    ax1 = axes[0, 0]
    plot_df = pd.DataFrame({'neighbors': neighbor_counts, 'expression': expression_values})
    plot_df['cell_count'] = plot_df.groupby(['neighbors', 'expression'])['neighbors'].transform('count')
    
    scatter = ax1.scatter(
        plot_df['neighbors'], plot_df['expression'],
        c=plot_df['cell_count'], cmap='viridis', alpha=0.9, s=25
    )
    cbar = fig.colorbar(scatter, ax=ax1)
    cbar.set_label('Cell Count')
    
    # Add regression line
    x_lims = ax1.get_xlim()
    x_vals = np.linspace(x_lims[0], x_lims[1], 100)
    ax1.plot(x_vals, intercept + slope * x_vals, color='red', linestyle='--')
    
    ax1.set_title(f"Expression vs. Neighbors (slope={slope:.3f}, r={r_value:.3f})", fontsize=10)
    ax1.set_xlabel("Number of Neighbors")
    ax1.set_ylabel("Gene Expression")
    ax1.grid(True, linestyle='--', linewidth=0.5)
    
    # --- Panel 2: Neighbor Count Histogram (top-right) ---
    ax2 = axes[0, 1]
    if len(neighbor_counts) > 0:
        ax2.hist(neighbor_counts, bins=range(int(neighbor_counts.max()) + 2),
                 alpha=0.7, edgecolor='black', color='lightblue')
        ax2.axvline(np.mean(neighbor_counts), color='red', linestyle='--',
                    label=f'Mean: {np.mean(neighbor_counts):.1f}')
        ax2.legend()
    ax2.set_title("Neighbor Count Distribution (Analysis Data)", fontsize=10)
    ax2.set_xlabel("Number of Neighbors")
    ax2.set_ylabel("Number of Target Cells")
    ax2.grid(True, linestyle='--', linewidth=0.5)
    
    # --- Panel 3: Mean Expression ± Std (bottom-left) ---
    ax3 = axes[1, 0]
    mean_expr = plot_df.groupby('neighbors')['expression'].mean()
    std_expr = plot_df.groupby('neighbors')['expression'].std()
    ax3.errorbar(mean_expr.index, mean_expr.values, yerr=std_expr.values, fmt='-o', capsize=5)
    ax3.set_title("Mean Expression vs. Neighbors", fontsize=10)
    ax3.set_xlabel("Number of Neighbor Cells")
    ax3.set_ylabel("Mean Expression ± Std Dev")
    ax3.grid(True, linestyle='--', linewidth=0.5)
    
    # --- Panel 4: Fit on Raw Data (bottom-right) ---
    ax4 = axes[1, 1]
    if original_neighbor_counts is not None and original_expression_values is not None:
        orig_plot_df = pd.DataFrame({
            'neighbors': original_neighbor_counts,
            'expression': original_expression_values
        })
        orig_plot_df['cell_count'] = orig_plot_df.groupby(
            ['neighbors', 'expression']
        )['neighbors'].transform('count')
        
        ax4.scatter(
            orig_plot_df['neighbors'], orig_plot_df['expression'],
            c=orig_plot_df['cell_count'], cmap='viridis', alpha=0.9, s=25
        )
        
        x_lims4 = ax4.get_xlim()
        x_vals4 = np.linspace(x_lims4[0], x_lims4[1], 100)
        ax4.plot(x_vals4, intercept + slope * x_vals4, color='red', linestyle='--')
        
        ax4.set_title("Fit on Raw Data", fontsize=10)
        ax4.set_xlabel("Number of Neighbors")
        ax4.set_ylabel("Gene Expression")
        ax4.grid(True, linestyle='--', linewidth=0.5)
    else:
        ax4.axis('off')
    
    fig.tight_layout()
    plt.subplots_adjust(left=0.08, right=0.95, top=0.9, bottom=0.08, hspace=0.35)
    
    return fig


def create_heatmap(
    results_df: pd.DataFrame,
    target_cell: str,
    neighbor_cell: str,
    fdr_alpha: float = 0.05,
    tissue_list: Optional[List[int]] = None,
) -> Optional[plt.Figure]:
    """
    Create heatmap showing -log10(q-value) for significant genes across tissues.
    
    Same as your GUI's heatmap.
    
    Parameters:
        results_df: DataFrame with results (must have 'gene', 'tissue', 'q_value')
        target_cell, neighbor_cell: Cell type names (for title)
        fdr_alpha: Significance threshold
        tissue_list: Optional list of tissues to show (default: all in data)
        
    Returns:
        matplotlib Figure, or None if no significant genes
    """
    # Filter to significant only
    sig_df = results_df[results_df["q_value"] < fdr_alpha].copy()
    
    if sig_df.empty:
        return None
    
    # Create pivot table: genes (rows) x tissues (columns)
    sig_df['neg_log_q'] = -np.log10(sig_df['q_value'].replace(0, 1e-300))
    heatmap_data = sig_df.pivot_table(index='gene', columns='tissue', values='neg_log_q')
    
    # Add missing tissues if tissue_list provided
    if tissue_list is not None:
        for tissue in tissue_list:
            if tissue not in heatmap_data.columns:
                heatmap_data[tissue] = np.nan
        heatmap_data = heatmap_data[sorted(tissue_list)]
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, max(6, len(heatmap_data.index) * 0.4)))
    im = ax.imshow(heatmap_data, cmap=HEATMAP_CMAP, aspect='auto')
    
    # Labels
    ax.set_xticks(np.arange(len(heatmap_data.columns)))
    ax.set_yticks(np.arange(len(heatmap_data.index)))
    ax.set_xticklabels(heatmap_data.columns)
    ax.set_yticklabels(heatmap_data.index)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
    
    # Grid
    ax.set_xticks(np.arange(heatmap_data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(heatmap_data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="black", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", size=0)
    
    # Colorbar
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label('-log10(q-value)', rotation=-90, va="bottom")
    
    ax.set_title(f"{target_cell} vs. {neighbor_cell}", fontsize=12)
    fig.tight_layout()
    plt.subplots_adjust(top=0.92, bottom=0.15, left=0.1, right=0.98)
    
    return fig


def save_plots(
    results_df: pd.DataFrame,
    out_dir: Path,
    method: str,
    target_cell: str,
    neighbor_cell: str,
    dist_threshold: float,
    fdr_alpha: float = 0.05,
    tissue_list: Optional[List[int]] = None,
) -> None:
    """
    Save all plots for significant genes with descriptive names.
    
    Creates:
    - plots/heatmap_{method}_{target}_vs_{neighbor}_dist{threshold}.png
    - plots/genes/{gene}_T{tissue}_{method}_{target}_vs_{neighbor}.png
    
    Parameters:
        results_df: DataFrame with results (must include plot data columns)
        out_dir: Output directory
        method: Analysis method (1, 2a, 2b, 3.1, 3.2)
        target_cell, neighbor_cell: Cell type names
        dist_threshold: Distance threshold used
        fdr_alpha: Significance threshold
        tissue_list: Optional list of tissues for heatmap
    """
    plots_dir = out_dir / "plots"
    genes_dir = plots_dir / "genes"
    plots_dir.mkdir(parents=True, exist_ok=True)
    genes_dir.mkdir(parents=True, exist_ok=True)
    
    # Descriptive base name
    base_name = f"method{method}_{target_cell}_vs_{neighbor_cell}_dist{int(dist_threshold)}"
    
    # Save heatmap
    heatmap_fig = create_heatmap(results_df, target_cell, neighbor_cell, fdr_alpha, tissue_list)
    if heatmap_fig is not None:
        heatmap_path = plots_dir / f"heatmap_{base_name}.png"
        heatmap_fig.savefig(heatmap_path, dpi=300, bbox_inches='tight')
        plt.close(heatmap_fig)
        print(f"Saved heatmap to {heatmap_path}")
    
    # Save individual gene plots
    sig_df = results_df[results_df["q_value"] < fdr_alpha]
    
    for _, row in sig_df.iterrows():
        if 'plot_neighbor_counts' not in row or 'plot_expression_values' not in row:
            continue
            
        gene_fig = create_gene_plot(
            gene_name=row['gene'],
            neighbor_counts=np.array(row['plot_neighbor_counts']),
            expression_values=np.array(row['plot_expression_values']),
            slope=row['slope'],
            intercept=row['intercept'],
            r_value=row['r_value'],
            tissue=row['tissue'],
            target_cell=target_cell,
            neighbor_cell=neighbor_cell,
            original_neighbor_counts=np.array(row.get('plot_original_neighbor_counts')) if 'plot_original_neighbor_counts' in row else None,
            original_expression_values=np.array(row.get('plot_original_expression_values')) if 'plot_original_expression_values' in row else None,
        )
        
        # Descriptive gene plot name
        gene_path = genes_dir / f"{row['gene']}_T{row['tissue']}_{base_name}.png"
        gene_fig.savefig(gene_path, dpi=300, bbox_inches='tight')
        plt.close(gene_fig)
    
    print(f"Saved {len(sig_df)} gene plots to {genes_dir}")


# =============================================================================
# GENE FILTERING FUNCTIONS (from original helper_functions.py)
# =============================================================================

def filter_genes_by_expression(
    expression_df: pd.DataFrame,
    genes: List[str],
    threshold: int,
) -> List[str]:
    """Keep genes whose MAX count > threshold."""
    expr = expression_df[genes]
    keep = expr.columns[expr.max(axis=0) > threshold].tolist()
    return keep


def filter_genes_by_variance(
    expression_df: pd.DataFrame,
    genes: List[str],
    quantile: float,
) -> List[str]:
    """Keep top X% most variable genes (by variance)."""
    expr = expression_df[genes]
    
    # Robust min-max normalization (same as original)
    def _robust_minmax(col: pd.Series) -> pd.Series:
        vmin = col.min()
        vmax = np.percentile(col, 95)
        if vmax == vmin:
            return pd.Series(np.nan, index=col.index)
        return (col - vmin) / (vmax - vmin)
    
    expr_norm = expr.apply(_robust_minmax, axis=0)
    var_ = expr_norm.var(axis=0, skipna=True)
    cutoff_quantile = 1.0 - quantile
    cutoff = var_.quantile(cutoff_quantile)
    
    keep = var_[var_ >= cutoff].index.tolist()
    return keep


def validate_bin_requirements(
    neighbor_counts: np.ndarray,
    min_bins: int,
    min_cells_per_bin: int,
) -> Tuple[bool, List[int]]:
    """Check if we have enough bins with enough cells for Method 2/3."""
    unique_counts, count_frequencies = np.unique(neighbor_counts, return_counts=True)
    
    valid_mask = count_frequencies >= min_cells_per_bin
    valid_counts = unique_counts[valid_mask]
    
    is_valid = len(valid_counts) >= min_bins
    return is_valid, sorted(valid_counts.tolist())


# =============================================================================
# FDR CORRECTION (from original helper_functions.py)
# =============================================================================

def apply_fdr_correction(p_values: np.ndarray, alpha: float = 0.05) -> np.ndarray:
    """Apply Benjamini-Hochberg FDR correction."""
    from statsmodels.stats.multitest import multipletests
    
    if len(p_values) == 0:
        return np.array([])
    
    _, q_values, _, _ = multipletests(p_values, alpha=alpha, method='fdr_bh')
    return q_values


# =============================================================================
# METHOD 1: LINEAR REGRESSION (from original analysis_methods.py)
# =============================================================================

def analyze_method1(
    expression_df: pd.DataFrame,
    neighbor_counts: np.ndarray,
    genes: List[str],
    tissue,
    target_cell: str,
    neighbor_cell: str,
    params: NbrsCountRegressionParams,
) -> List[dict]:
    """
    Method 1: Standard linear regression for each gene.
    
    For each gene: regress expression ~ neighbor_counts
    """
    from scipy import stats
    
    results = []
    
    for gene in genes:
        expression_values = expression_df[gene].values
        
        # Skip if no variance
        if np.all(expression_values == expression_values[0]):
            continue
        
        try:
            slope, intercept, r_value, p_value, _ = stats.linregress(
                neighbor_counts, expression_values
            )
            
            results.append({
                'tissue': tissue,
                'target_cell': target_cell,
                'neighbor_cell': neighbor_cell,
                'gene': gene,
                'method': '1',
                'p_value': p_value,
                'r_value': r_value,
                'r_squared': r_value ** 2,
                'slope': slope,
                'intercept': intercept,
                'n_cells': len(neighbor_counts),
                'plot_neighbor_counts': neighbor_counts.tolist(),
                'plot_expression_values': expression_values.tolist(),
            })
        except Exception:
            continue
    
    return results


# =============================================================================
# METHOD 2: SAMPLING-BASED REGRESSION (from original analysis_methods.py)
# =============================================================================

def analyze_method2(
    expression_df: pd.DataFrame,
    neighbor_counts: np.ndarray,
    genes: List[str],
    tissue,
    target_cell: str,
    neighbor_cell: str,
    params: NbrsCountRegressionParams,
) -> List[dict]:
    """
    Method 2: Sampling-based regression.
    
    Returns results for both 2a (slope averaging) and 2b (line-on-line).
    """
    from scipy import stats
    
    # Validate bin requirements
    is_valid, valid_bins = validate_bin_requirements(
        neighbor_counts, params.min_neighbor_bins, params.min_cells_per_bin
    )
    if not is_valid:
        return []
    
    # Filter to valid bins only
    valid_mask = np.isin(neighbor_counts, valid_bins)
    filtered_df = expression_df[valid_mask].copy()
    filtered_counts = neighbor_counts[valid_mask]
    
    if len(filtered_df) == 0:
        return []
    
    # Prepare sampling data
    sampling_data = {nc: np.where(filtered_counts == nc)[0] for nc in valid_bins}
    
    # Run sampling iterations
    slopes_dict = {gene: [] for gene in genes}
    intercepts_dict = {gene: [] for gene in genes}
    
    for _ in range(params.n_iterations):
        sampled_counts = []
        sampled_indices = []
        
        for nc in valid_bins:
            available = sampling_data[nc]
            if len(available) == 0:
                continue
            sampled_idx = np.random.choice(available, size=params.sample_size, replace=True)
            sampled_counts.extend([nc] * params.sample_size)
            sampled_indices.extend(sampled_idx)
        
        for gene in genes:
            try:
                expr_vals = filtered_df[gene].iloc[sampled_indices].values
                if np.all(expr_vals == expr_vals[0]):
                    continue
                slope, intercept, _, _, _ = stats.linregress(sampled_counts, expr_vals)
                slopes_dict[gene].append(slope)
                intercepts_dict[gene].append(intercept)
            except Exception:
                continue
    
    # Compute results for 2a and 2b
    results = []
    max_nc = int(np.max(filtered_counts))
    evaluation_range = np.arange(0, max_nc + 1)
    
    for gene in genes:
        slopes = slopes_dict.get(gene, [])
        intercepts = intercepts_dict.get(gene, [])
        
        if not slopes:
            continue
        
        expression_values = filtered_df[gene].values
        raw_expression = expression_df[gene].values
        
        # --- Method 2a: Slope averaging ---
        mean_slope = np.mean(slopes)
        optimal_intercept = np.mean(expression_values - mean_slope * filtered_counts)
        
        # Calculate correlation
        y_actual = expression_values
        x = filtered_counts
        numerator = np.sum((x - np.mean(x)) * (y_actual - np.mean(y_actual)))
        denominator = np.sqrt(np.sum((x - np.mean(x))**2) * np.sum((y_actual - np.mean(y_actual))**2))
        
        if denominator == 0:
            r_value_2a, p_value_2a = 0.0, 1.0
        else:
            r_value_2a = numerator / denominator
            n = len(filtered_counts)
            if abs(r_value_2a) == 1.0:
                p_value_2a = 0.0
            else:
                t_stat = r_value_2a * np.sqrt((n - 2) / (1 - r_value_2a**2))
                p_value_2a = 2 * stats.t.sf(abs(t_stat), df=n - 2)
        
        results.append({
            'tissue': tissue,
            'target_cell': target_cell,
            'neighbor_cell': neighbor_cell,
            'gene': gene,
            'method': '2a',
            'p_value': p_value_2a,
            'r_value': r_value_2a,
            'r_squared': r_value_2a ** 2,
            'slope': mean_slope,
            'intercept': optimal_intercept,
            'n_cells': len(filtered_counts),
            'plot_neighbor_counts': filtered_counts.tolist(),
            'plot_expression_values': expression_values.tolist(),
            'plot_original_neighbor_counts': neighbor_counts.tolist(),
            'plot_original_expression_values': raw_expression.tolist(),
        })
        
        # --- Method 2b: Line-on-line regression ---
        if len(slopes) == len(intercepts) and len(slopes) > 0:
            all_x = np.tile(evaluation_range, len(slopes))
            all_y = np.concatenate([s * evaluation_range + i for s, i in zip(slopes, intercepts)])
            
            if len(all_x) >= 2 and not np.all(all_y == all_y[0]):
                slope_2b, intercept_2b, r_value_2b, p_value_2b, _ = stats.linregress(all_x, all_y)
                
                results.append({
                    'tissue': tissue,
                    'target_cell': target_cell,
                    'neighbor_cell': neighbor_cell,
                    'gene': gene,
                    'method': '2b',
                    'p_value': p_value_2b,
                    'r_value': r_value_2b,
                    'r_squared': r_value_2b ** 2,
                    'slope': slope_2b,
                    'intercept': intercept_2b,
                    'n_cells': len(filtered_counts),
                    'plot_neighbor_counts': filtered_counts.tolist(),
                    'plot_expression_values': expression_values.tolist(),
                    'plot_original_neighbor_counts': neighbor_counts.tolist(),
                    'plot_original_expression_values': raw_expression.tolist(),
                })
    
    return results


# =============================================================================
# METHOD 3: WEIGHTED LEAST SQUARES (from original analysis_methods.py)
# =============================================================================

def analyze_method3(
    expression_df: pd.DataFrame,
    neighbor_counts: np.ndarray,
    genes: List[str],
    tissue,
    target_cell: str,
    neighbor_cell: str,
    params: NbrsCountRegressionParams,
) -> List[dict]:
    """
    Method 3: Weighted least squares regression.
    
    Returns results for both 3.1 (unnormalized) and 3.2 (normalized weights).
    """
    import statsmodels.api as sm
    
    # Validate bin requirements
    is_valid, valid_bins = validate_bin_requirements(
        neighbor_counts, params.min_neighbor_bins, params.min_cells_per_bin
    )
    if not is_valid:
        return []
    
    # Filter to valid bins only
    valid_mask = np.isin(neighbor_counts, valid_bins)
    filtered_df = expression_df[valid_mask].copy()
    filtered_counts = neighbor_counts[valid_mask]
    
    if len(filtered_df) == 0:
        return []
    
    # Compute weights
    def compute_weights(counts, normalize=False):
        unique, freq = np.unique(counts, return_counts=True)
        n_total = len(counts)
        probs = {int(c): float(f) / float(n_total) for c, f in zip(unique, freq)}
        weight_map = {c: 1.0 / p for c, p in probs.items()}
        weights = np.array([weight_map[int(c)] for c in counts], dtype=float)
        if normalize:
            s = weights.sum()
            if s > 0:
                weights /= s
        return weights
    
    weights_unnorm = compute_weights(filtered_counts, normalize=False)
    weights_norm = compute_weights(filtered_counts, normalize=True)
    
    results = []
    
    for gene in genes:
        expression_values = filtered_df[gene].values
        raw_expression = expression_df[gene].values
        
        X = sm.add_constant(filtered_counts)
        
        # --- Method 3.1: Unnormalized weights ---
        try:
            wls_model = sm.WLS(expression_values, X, weights=weights_unnorm).fit()
            
            results.append({
                'tissue': tissue,
                'target_cell': target_cell,
                'neighbor_cell': neighbor_cell,
                'gene': gene,
                'method': '3.1',
                'p_value': float(wls_model.pvalues[1]),
                'r_value': float(np.sign(wls_model.params[1]) * np.sqrt(max(wls_model.rsquared, 0.0))),
                'r_squared': float(max(wls_model.rsquared, 0.0)),
                'slope': float(wls_model.params[1]),
                'intercept': float(wls_model.params[0]),
                'n_cells': len(filtered_counts),
                'plot_neighbor_counts': filtered_counts.tolist(),
                'plot_expression_values': expression_values.tolist(),
                'plot_original_neighbor_counts': neighbor_counts.tolist(),
                'plot_original_expression_values': raw_expression.tolist(),
            })
        except Exception:
            pass
        
        # --- Method 3.2: Normalized weights ---
        try:
            wls_model = sm.WLS(expression_values, X, weights=weights_norm).fit()
            
            results.append({
                'tissue': tissue,
                'target_cell': target_cell,
                'neighbor_cell': neighbor_cell,
                'gene': gene,
                'method': '3.2',
                'p_value': float(wls_model.pvalues[1]),
                'r_value': float(np.sign(wls_model.params[1]) * np.sqrt(max(wls_model.rsquared, 0.0))),
                'r_squared': float(max(wls_model.rsquared, 0.0)),
                'slope': float(wls_model.params[1]),
                'intercept': float(wls_model.params[0]),
                'n_cells': len(filtered_counts),
                'plot_neighbor_counts': filtered_counts.tolist(),
                'plot_expression_values': expression_values.tolist(),
                'plot_original_neighbor_counts': neighbor_counts.tolist(),
                'plot_original_expression_values': raw_expression.tolist(),
            })
        except Exception:
            pass
    
    return results


# =============================================================================
# MAIN RUN FUNCTION
# =============================================================================

def run_analysis(
    labels: pd.DataFrame,
    genes: List[str],
    tissue,
    target_cell: str,
    neighbor_cell: str,
    method: str,
    params: NbrsCountRegressionParams,
    distance_matrix: Optional[pd.DataFrame] = None,
) -> pd.DataFrame:
    """
    Run analysis for one cell type pair.
    
    Parameters:
        labels: Full data DataFrame
        genes: List of gene names to analyze
        tissue: Tissue ID
        target_cell, neighbor_cell: Cell type names
        method: '1', '2a', '2b', '3.1', or '3.2' (or '2' for both 2a/2b, '3' for both 3.1/3.2)
        params: Analysis parameters
        distance_matrix: Optional pre-computed distance matrix
        
    Returns:
        DataFrame with results
    """
    # Filter to tissue
    labels_tissue = labels[labels["tissue"].astype(str) == str(tissue)].copy()
    labels_tissue = labels_tissue.reset_index(drop=True)
    
    # Get cells by type
    primary_mask = labels_tissue["cell_type"] == target_cell
    neighbor_mask = labels_tissue["cell_type"] == neighbor_cell
    
    sub_primary = labels_tissue[primary_mask].copy()
    sub_neighbor = labels_tissue[neighbor_mask].copy()
    
    if len(sub_primary) == 0 or len(sub_neighbor) == 0:
        return pd.DataFrame()
    
    # For distance-matrix mode, set cell_id as index so .loc lookups work
    if distance_matrix is not None:
        sub_primary_dm = sub_primary.set_index("cell_id")
        sub_neighbor_dm = sub_neighbor.set_index("cell_id")
    else:
        sub_primary_dm = sub_primary
        sub_neighbor_dm = sub_neighbor

    # Inform user which proximity source is used
    if distance_matrix is not None:
        print(f"[INFO] Tissue {tissue}: using distance matrix for neighbor counting.")
    else:
        print(f"[INFO] Tissue {tissue}: using coordinate-based distances (cKDTree).")

    # Calculate neighbor counts
    neighbor_counts = calculate_neighbor_counts(
        sub_primary_dm,
        sub_neighbor_dm,
        params.dist_threshold,
        distance_matrix,
    )
    
    if len(np.unique(neighbor_counts)) <= 1:
        return pd.DataFrame()
    
    # Filter genes
    filtered_genes = filter_genes_by_expression(sub_primary, genes, params.min_expression)
    filtered_genes = filter_genes_by_variance(sub_primary, filtered_genes, params.variance_quantile)
    
    if not filtered_genes:
        return pd.DataFrame()
    
    # Run appropriate method
    if method == '1':
        results = analyze_method1(
            sub_primary, neighbor_counts, filtered_genes,
            tissue, target_cell, neighbor_cell, params
        )
    elif method in ['2', '2a', '2b']:
        results = analyze_method2(
            sub_primary, neighbor_counts, filtered_genes,
            tissue, target_cell, neighbor_cell, params
        )
        # Filter to specific sub-method if requested
        if method in ['2a', '2b']:
            results = [r for r in results if r['method'] == method]
    elif method in ['3', '3.1', '3.2']:
        results = analyze_method3(
            sub_primary, neighbor_counts, filtered_genes,
            tissue, target_cell, neighbor_cell, params
        )
        # Filter to specific sub-method if requested
        if method in ['3.1', '3.2']:
            results = [r for r in results if r['method'] == method]
    else:
        raise ValueError(f"Unknown method: {method}")
    
    if not results:
        return pd.DataFrame()
    
    # Convert to DataFrame and apply FDR correction
    results_df = pd.DataFrame(results)
    
    if 'p_value' in results_df.columns and len(results_df) > 0:
        results_df['q_value'] = apply_fdr_correction(
            results_df['p_value'].values, params.fdr_alpha
        )
    
    return results_df
