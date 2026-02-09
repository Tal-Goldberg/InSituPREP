from __future__ import annotations

from pathlib import Path
from typing import List, Optional

import typer
import json
from typing import Literal

ProximityMode = Literal["continuous", "binary"]


from insituprep.analyses.rna_velocity_stage1 import (
    stage1_gamma_grid,
    stage1_gamma_plots_all_genes,
    stage1_find_best_l_fscore,
    stage1_t_preview,
)

from insituprep.analyses.rna_velocity_stage2 import run_velocity_stage2

app = typer.Typer(add_completion=True, help="RNA velocity analysis (stages 1+2).")


@app.command("stage1-a", help="Scan K and quantile grids and compute per-gene gamma values (slope u_knn vs s_knn).")
def stage1_a_gamma_grid(
    pc_origin_path: Path = typer.Option(
        ...,
        "--pc-origin-path",
        exists=True,
        help=(
            "CSV with PCA embedding of cells. "
            "Index must be cell IDs. "
            "Required columns: PC1, PC2, PC3. "
            "Optional column: cell_type (used for coloring gamma plots)."
        ),
    ),
    spliced_counts_path: Path = typer.Option(
        ...,
        "--spliced-counts-path",
        exists=True,
        help=(
            "CSV with cytoplasmic (spliced) RNA counts. "
            "Rows correspond to individual cells. "
            "The first column must be 'tissue' (tissue ID), "
            "followed by gene expression count columns (cells x genes)."
        ),
    ),
    unspliced_counts_path: Path = typer.Option(
        ...,
        "--unspliced-counts-path",
        exists=True,
        help=(
            "CSV with nuclear (unspliced) RNA counts. "
            "Rows correspond to individual cells. "
            "The first column must be 'tissue' (tissue ID), "
            "followed by gene expression count columns (cells x genes)."
        ),
    ),
    out: Path = typer.Option(..., "--out", help="Output directory for RNA velocity results."),
    k_values: str = typer.Option(..., "--k-values", help='JSON list of k nearest neighbors values to scan, e.g. "[5,10,20,30,40,50]".'),
    quantiles: str = typer.Option(..., "--quantiles", help='JSON list of quantile outliers values to scan, e.g. "[0.05,0.075,0.1]".')
):

    print("\n[RNA-velocity] stage1-a START")
    print(f"  out      : {out}")
    
    k_values_list = json.loads(k_values)
    quantiles_list = json.loads(quantiles)

    # basic sanity
    if not isinstance(k_values_list, list) or not all(isinstance(x, int) for x in k_values_list):
        raise typer.BadParameter("Expected --k-values to be a list of ints, e.g. [10,20,30]")
    if not isinstance(quantiles_list, list) or not all(isinstance(x, (int, float)) for x in quantiles_list):
        raise typer.BadParameter("Expected --quantiles to be a list of floats, e.g. [0.05,0.075]")


    print(f"  k_values (parsed) : {k_values_list}")
    print(f"  quantiles(parsed) : {quantiles_list}")
    
    stage1_gamma_grid(
        pc_origin_path=pc_origin_path,
        spliced_counts_path=spliced_counts_path,
        unspliced_counts_path=unspliced_counts_path,
        out_dir=out,
        k_values=k_values_list,
        quantiles=quantiles_list,
        l_fixed=0.038,
    )

    print("[RNA-velocity] stage1-a DONE\n")

@app.command("stage1-b", help="Generate per-gene gamma portrait plots (scatter by cell type + fitted gamma line) for chosen K and quantile.")
def stage1_b_gamma_plots(
    pc_origin_path: Path = typer.Option(
        ...,
        "--pc-origin-path",
        exists=True,
        help=(
            "CSV with PCA embedding of cells. "
            "Index must be cell IDs. "
            "Required columns: PC1, PC2, PC3. "
            "Optional column: cell_type (used for coloring gamma plots)."
        ),
    ),
    spliced_counts_path: Path = typer.Option(
        ...,
        "--spliced-counts-path",
        exists=True,
        help=(
            "CSV with cytoplasmic (spliced) RNA counts. "
            "Rows correspond to individual cells. "
            "The first column must be 'tissue' (tissue ID), "
            "followed by gene expression count columns (cells x genes)."
        ),
    ),
    unspliced_counts_path: Path = typer.Option(
        ...,
        "--unspliced-counts-path",
        exists=True,
        help=(
            "CSV with nuclear (unspliced) RNA counts. "
            "Rows correspond to individual cells. "
            "The first column must be 'tissue' (tissue ID), "
            "followed by gene expression count columns (cells x genes)."
        ),
    ),
    out: Path = typer.Option(..., "--out", help="Output directory for RNA velocity results."),
    k: int = typer.Option(..., "--k", help = "Chosen K (number of PCA nearest neighbors) used for KNN pooling and gamma/velocity calculations."),
    quantile: float = typer.Option(..., "--quantile", help="Chosen quantile outliers used for gamma fitting (extreme-point regression).")
):
    print("\n[RNA-velocity] stage1-b START")
    print(f"  out      : {out}")
    print(f"  k        : {k}")
    print(f"  quantile : {quantile}")


    stage1_gamma_plots_all_genes(
        pc_origin_path=pc_origin_path,
        spliced_counts_path=spliced_counts_path,
        unspliced_counts_path=unspliced_counts_path,
        out_dir=out,
        k=k,
        quantile=quantile,
        l_fixed=0.038,
    )

    print("[RNA-velocity] stage1-b DONE\n")


@app.command("stage1-c", help="Find best_l by scanning expression thresholds and maximizing precision/recall F-score vs ground-truth genes. Writes stage1_selected_params.txt.")
def stage1_c_best_l(
    pc_origin_path: Path = typer.Option(
        ...,
        "--pc-origin-path",
        exists=True,
        help=(
            "CSV with PCA embedding of cells. "
            "Index must be cell IDs. "
            "Required columns: PC1, PC2, PC3. "
            "Optional column: cell_type (used for coloring gamma plots)."
        ),
    ),
    spliced_counts_path: Path = typer.Option(
        ...,
        "--spliced-counts-path",
        exists=True,
        help=(
            "CSV with cytoplasmic (spliced) RNA counts. "
            "Rows correspond to individual cells. "
            "The first column must be 'tissue' (tissue ID), "
            "followed by gene expression count columns (cells x genes)."
        ),
    ),
    unspliced_counts_path: Path = typer.Option(
        ...,
        "--unspliced-counts-path",
        exists=True,
        help=(
            "CSV with nuclear (unspliced) RNA counts. "
            "Rows correspond to individual cells. "
            "The first column must be 'tissue' (tissue ID), "
            "followed by gene expression count columns (cells x genes)."
        ),
    ),
    out: Path = typer.Option(..., "--out", help="Output directory for RNA velocity results."),
    k: int = typer.Option(..., "--k", help = "Chosen K (number of PCA nearest neighbors) used for KNN pooling and gamma/velocity calculations."),
    quantile: float = typer.Option(..., "--quantile", help="Chosen quantile outliers used for gamma fitting (extreme-point regression)."),
    ground_truth_path: Path = typer.Option(..., "--ground-truth-path", help="Path to ground_truth_genes.txt (one gene per line). "
                                                                            "Reference genes used for parameter tuning in stage1-c. "
                                                                            "These genes were selected by visually inspecting the phase portrait of each gene "
                                                                            "and identifying genes with the expected shape and a clear linear trend "
                                                                            "in the spliced–unspliced RNA velocity space."),
):
    print("\n[RNA-velocity] stage1-c START (best_l)")
    print(f"  out               : {out}")
    print(f"  k                 : {k}")
    print(f"  quantile          : {quantile}")
    print(f"  ground_truth_path : {ground_truth_path}")

    best_l = stage1_find_best_l_fscore(
        pc_origin_path=pc_origin_path,
        spliced_counts_path=spliced_counts_path,
        unspliced_counts_path=unspliced_counts_path,
        out_dir=out,
        k=k,
        quantile=quantile,
        ground_truth_txt=ground_truth_path,
    )
    typer.echo(f"best_l = {best_l} (saved to {out/'stage1_selected_params.txt'})")
    print("[RNA-velocity] stage1-c DONE\n")

@app.command("stage1-d", help="Compute velocity and predicted future spliced states s(t) across a grid of t values. Saves spliced_values_in_multiple_t.pkl and PCA arrow plots.")
def stage1_d_t_preview(
    pc_origin_path: Path = typer.Option(
        ...,
        "--pc-origin-path",
        exists=True,
        help=(
            "CSV with PCA embedding of cells. "
            "Index must be cell IDs. "
            "Required columns: PC1, PC2, PC3. "
            "Optional column: cell_type (used for coloring gamma plots)."
        ),
    ),
    spliced_counts_path: Path = typer.Option(
        ...,
        "--spliced-counts-path",
        exists=True,
        help=(
            "CSV with cytoplasmic (spliced) RNA counts. "
            "Rows correspond to individual cells. "
            "The first column must be 'tissue' (tissue ID), "
            "followed by gene expression count columns (cells x genes)."
        ),
    ),
    unspliced_counts_path: Path = typer.Option(
        ...,
        "--unspliced-counts-path",
        exists=True,
        help=(
            "CSV with nuclear (unspliced) RNA counts. "
            "Rows correspond to individual cells. "
            "The first column must be 'tissue' (tissue ID), "
            "followed by gene expression count columns (cells x genes)."
        ),
    ),
    out: Path = typer.Option(..., "--out", help="Output directory for RNA velocity results."),
    k: int = typer.Option(..., "--k", help = "Chosen K (number of PCA nearest neighbors) used for KNN pooling and gamma/velocity calculations."),
    quantile: float = typer.Option(..., "--quantile", help="Chosen quantile outliers used for gamma fitting (extreme-point regression)."),
    ground_truth_path: Path = typer.Option(..., "--ground-truth-path", help="Path to ground_truth_genes.txt (one gene per line). "
                                                                            "Reference genes used for parameter tuning in stage1-c. "
                                                                            "These genes were selected by visually inspecting the phase portrait of each gene "
                                                                            "and identifying genes with the expected shape and a clear linear trend "
                                                                            "in the spliced–unspliced RNA velocity space."),
    best_l: float = typer.Option(..., "--best-l", help="Normalized-expression threshold chosen in stage1-c, defined as the expression cutoff "
                                                        "that maximizes the precision–recall F1 score. "
                                                        "Controls which genes are retained in the filtered velocity and s(t) calculations."),
    n_cells: int = typer.Option(500, "--n-cells", help= "Number of cells to sample for PCA arrow plots (t_figures). Use e.g. 500 for preview."),
 ):
    print("\n[RNA-velocity] stage1-d START (t preview)")
    print(f"  out               : {out}")
    print(f"  k                 : {k}")
    print(f"  quantile          : {quantile}")
    print(f"  best_l            : {best_l}")
    print(f"  n_cells (arrows)  : {n_cells}")
    print(f"  ground_truth_path : {ground_truth_path}")


    stage1_t_preview(
        pc_origin_path=pc_origin_path,
        spliced_counts_path=spliced_counts_path,
        unspliced_counts_path=unspliced_counts_path,
        out_dir=out,
        k=k,
        quantile=quantile,
        ground_truth_txt=ground_truth_path,
        best_l=best_l,
        n_cells=n_cells,
    )
    print("[RNA-velocity] stage1-d DONE\n")


@app.command(
    "stage2",
    help=
    "Stage 2: loads the Stage1 outputs (predicted s(t)/velocity-derived features stored in spliced_values_in_multiple_t.pkl) and computes per-cell proximity to each cell type "
    "(continuous min-distance and a binary thresholded version), then tests primary→neighbor associations between proximity and the "
    "Stage1-derived velocity geometry features (magnitude/phase; optional permutation p-values + BH-FDR), optionally repeats the tests "
    "per tissue, and optionally runs a gene-level analysis correlating |velocity_gene| with proximity for significant primary→neighbor pairs."
)
def stage2_run(
    pc_origin_path: Path = typer.Option(
        ...,
        "--pc-origin-path",
        exists=True,
        help="CSV with PCA embedding of cells. Index must be cell IDs. Required: PC1, PC2. Optional: PC3, cell_type."
    ),

    summary_table_path: Path = typer.Option(
    ...,
    "--summary-table-path",
    exists=True,
        help=("CSV with per-cell metadata including cell centroid coordinates (µm) and gene expression.\n\n"
               "Required columns: 'Var1' (cell IDs), 'tissue', 'cell_type'.\n\n"
                "The table may also include gene expression columns "
                "(one column per gene, raw numeric count values) but these columns are ignored (since we use the 'spliced' and 'unspliced' counts matrices instead).\n\n"
                "Proximity source:\n\n"
                "  - If --distance-dir is provided: proximity is computed from the per-tissue distance matrices.\n\n"
                "  - If --distance-dir is NOT provided: the summary table MUST include centroid coordinates:\n\n"
                "      * Required: 'X_space', 'Y_space'\n\n"
                "      * Optional: 'Z_space' (if present, 3D distances are used; otherwise 2D).\n\n"),
    ),

    distance_dir: Optional[Path] = typer.Option(
        None,
        "--distance-dir",
        help="Optional directory with per-tissue cell-cell distance matrices "
            "('distance_matrix_{tissue}.csv'), containing pairwise cell–cell distances (µm). "
            "Rows and columns correspond to cell IDs. "
            "If omitted, proximity is computed from cell coordinates in the labels file "
            "(X_space, Y_space, and optionally Z_space). "
    ),
    gene_dict_pkl: Path = typer.Option(
        ...,
        "--gene-dict-pkl",
        exists=True,
        help="Pickle produced by stage1 (eg. spliced_values_in_multiple_t.pkl or a dict containing 'all_gene_level_dfs_per_t')."
    ),
    out: Path = typer.Option(..., "--out", help="Output directory for stage2 results."),
    t: float = typer.Option(3.0, "--t", help="Timepoint t to analyze (must exist as key in the gene dict; e.g. 3.0)."),
    proximity_mode: ProximityMode = typer.Option("continuous","--proximity-mode", case_sensitive=False, help="Which proximity table to analyze.", show_choices=True),
    distance_threshold: float = typer.Option(1.0, "--distance-threshold", help="Distance threshold (µm) used to build the binary proximity table."),
    maximum_distance_threshold: float = typer.Option(float("inf"), "--maximum-distance-threshold", help="For XY/Z mode only: distances greater than this threshold "
                                                     "are set to NaN and excluded from downstream analyses."),
    min_cells: int = typer.Option(10, "--min-cells", help="Minimum number of primary cells required to test a primary→neighbor pair."),
    n_perm: int = typer.Option(10000, "--n-perm", help="Number of permutations for pairwise permutation testing."),
    rng_seed: Optional[int] = typer.Option(None, "--rng-seed", help="Random seed for permutation tests. "
                                                                        "If provided (e.g. 0), results are reproducible. "
                                                                        "If omitted, permutations are randomized on each run."),
    fdr_alpha: float = typer.Option(0.05, "--fdr-alpha", help="FDR alpha (BH) for reporting q-values."),
    compute_pair_permutations: bool = typer.Option(True, "--compute-pair-permutations/--no-compute-pair-permutations", help="Run pairwise permutations (slower)."),
    compute_tissue: bool = typer.Option(True, "--compute-tissue/--no-compute-tissue", help="Also run tissue-stratified per-tissue pairwise analysis."),
    compute_gene_level: bool = typer.Option(True, "--compute-gene-level/--no-compute-gene-level", help="Run gene-level velocity↔proximity analysis for significant pairs."),
    gene_from_feature: str = typer.Option("magnitude_delta", "--gene-from-feature", help="Which pairwise feature to select significant pairs from for gene-level analysis."),
    gene_pval_col: str = typer.Option("p_value_perm_fdr", "--gene-pval-col", help="Which p-value column to use when selecting significant pairs for gene-level analysis."),
    gene_alpha: float = typer.Option(0.05, "--gene-alpha", help="Alpha threshold when selecting significant pairs for gene-level analysis."),
    run_gene_permutations: bool = typer.Option(False, "--run-gene-permutations/--no-run-gene-permutations", help="If enabled, compute gene-level permutation p-values (very slow)."),
    gene_n_perm: int = typer.Option(10000, "--gene-n-perm", help="Number of permutations for per-gene permutations (if enabled)."),
    gene_rng_seed: Optional[int] = typer.Option(None, "--gene-rng-seed", help="RNG seed for gene-level permutations."),
    verbose: bool = typer.Option(True, "--verbose/--no-verbose", help="Verbose logging."),
):
    """
    Run the Stage2 pairwise velocity↔proximity pipeline.

    Typical usage notes:
     - gene_dict_pkl should be the pickle produced by stage1 (spliced_values_in_multiple_t.pkl)
       or a dict that contains key 'all_gene_level_dfs_per_t'.
     - Make sure the selected t exists in the pickle keys (e.g. 3.0 or '3.0').
    """
    print("\n[RNA-velocity] stage2 START")
    print(f"  out: {out}")
    print(f"  t: {t} | proximity_mode: {proximity_mode} | distance_threshold: {distance_threshold}")
    print(f"  gene_dict_pkl: {gene_dict_pkl}")
    print(f"  summary_table_path: {summary_table_path} | distance_dir: {distance_dir}")

    # call wrapper
    run_velocity_stage2(
        pc_origin_path=pc_origin_path,
        summary_table_path=summary_table_path,
        distance_dir=distance_dir,
        gene_dict_pkl=gene_dict_pkl,
        out_dir=out,
        t=t,
        proximity_mode=proximity_mode,  # "continuous" or "binary"
        distance_threshold=distance_threshold,
        maximum_distance_threshold=maximum_distance_threshold,
        min_cells=min_cells,
        n_perm=n_perm,
        rng_seed=rng_seed,
        fdr_alpha=fdr_alpha,
        compute_pair_permutations=compute_pair_permutations,
        compute_tissue=compute_tissue,
        compute_gene_level=compute_gene_level,
        gene_from_feature=gene_from_feature,
        gene_pval_col=gene_pval_col,
        gene_alpha=gene_alpha,
        run_gene_permutations=run_gene_permutations,
        gene_n_perm=gene_n_perm,
        gene_rng_seed=gene_rng_seed,
        verbose=verbose,
    )
    print("[RNA-velocity] stage2 DONE\n")
