from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, List, Tuple, Dict
from multiprocessing import Pool, cpu_count

import numpy as np
import pandas as pd
import typer

from insituprep.analyses.nbrs_count_regression import (
    NbrsCountRegressionParams,
    load_labels,
    load_distance_matrix,
    get_gene_columns,
    get_cell_type_pairs,
    get_tissues,
    run_analysis,
    save_results_csv,
    save_plots,
    create_heatmap,
)

# Global variables for multiprocessing (set by _init_worker)
_mp_labels = None
_mp_genes = None
_mp_tissue = None
_mp_method = None
_mp_params = None
_mp_distance_matrix = None
_mp_out = None


def _init_worker(labels, genes, tissue, method, params, distance_matrix, out):
    global _mp_labels, _mp_genes, _mp_tissue, _mp_method, _mp_params, _mp_distance_matrix, _mp_out
    _mp_labels = labels
    _mp_genes = genes
    _mp_tissue = tissue
    _mp_method = method
    _mp_params = params
    _mp_distance_matrix = distance_matrix
    _mp_out = out


def _process_pair(pair: Tuple[str, str]):
    target_cell, neighbor_cell = pair

    results_df = run_analysis(
        labels=_mp_labels,
        genes=_mp_genes,
        tissue=_mp_tissue,
        target_cell=target_cell,
        neighbor_cell=neighbor_cell,
        method=_mp_method,
        params=_mp_params,
        distance_matrix=_mp_distance_matrix,
    )

    if len(results_df) > 0:
        save_results_csv(
            results_df, _mp_out, _mp_method, target_cell, neighbor_cell,
            _mp_params.dist_threshold, _mp_params.fdr_alpha
        )
        save_plots(
            results_df, _mp_out, _mp_method, target_cell, neighbor_cell,
            _mp_params.dist_threshold, _mp_params.fdr_alpha
        )
        return results_df
    return None


app = typer.Typer(
    add_completion=True,
    help=(
        "Neighbor Count Regression Analysis.\n"
        "Analyzes the relationship between gene expression and number of neighboring cells.\n"
        "Implements three regression methods: (1) linear, (2a/2b) sampling-based, (3.1/3.2) weighted."
    ),
)


def _parse_json_str_list(s: Optional[str], opt_name: str) -> Optional[List[str]]:
    """Parse JSON list of strings (cell types). Returns list[str] or None."""
    if s is None:
        return None
    try:
        x = json.loads(s)
    except Exception as e:
        raise typer.BadParameter(f"{opt_name} must be a valid JSON list.") from e
    if not isinstance(x, list):
        raise typer.BadParameter(f"{opt_name} must be a JSON list.")
    # enforce strings (but keep exact text)
    return [str(v) for v in x]


def _parse_json_list(s: Optional[str], opt_name: str) -> Optional[List[str]]:
    if s is None:
        return None
    try:
        x = json.loads(s)
    except Exception as e:
        raise typer.BadParameter(f"{opt_name} must be a valid JSON list.") from e
    if not isinstance(x, list):
        raise typer.BadParameter(f"{opt_name} must be a JSON list.")
    return [str(v) for v in x]


@app.command("run",
    help=(
        "Neighbor Count Regression Analysis.\n"
        "Analyzes the relationship between gene expression and number of neighboring cells.\n"
        "Implements three regression methods: (1) linear, (2a/2b) sampling-based, (3.1/3.2) weighted."
             ))
def run(
    summary_table_path: Path = typer.Option(
        ..., "--summary-table-path", "-s",
        exists=True,
        help=("CSV with per-cell metadata including cell centroid coordinates (µm) and gene expression.\n\n"
               "Required columns: 'Var1' (cell IDs), 'tissue', 'cell_type'.\n\n"
                "The table should also include gene expression columns "
                "(one column per gene, raw numeric count values).\n\n"
                "Proximity source:\n\n"
                "  - If --distance-dir is provided: proximity is computed from the per-tissue distance matrices.\n\n"
                "  - If --distance-dir is NOT provided: the summary table MUST include centroid coordinates:\n\n"
                "      * Required: 'X_space', 'Y_space'\n\n"
                "      * Optional: 'Z_space' (if present, 3D distances are used; otherwise 2D).\n\n"),
    ),
    genes_names_path: Path = typer.Option(
        ..., "--genes-names-path", "-g",
        exists=True,
        help="Text file containing one gene name per line. "
            "Used to select a subset of genes and extract the corresponding "
            "gene expression columns from the 'summary_table' CSV."
    ),
    out: Path = typer.Option(
        ..., "--out", "-o",
        help="Output directory for results and plots.",
    ),
    tissue: Optional[str] = typer.Option(
        None, "--tissue", "-t",
        help='JSON list of tissue IDs.\n\n Example: \'[100,313]\'.\n\n If omitted: run all tissues.',
    ),

    # Distance options
    distance_dir: Optional[Path] = typer.Option(
        None, "--distance-dir",
        exists=True,
        file_okay=False,
        help="Directory with distance matrices. If omitted, calculates from coordinates.",
    ),
    dist_threshold: float = typer.Option(
        40.0, "--dist-threshold",
        help="Distance threshold in microns.",
    ),

    # Method selection
    method: str = typer.Option(
        "1", "--method", "-m",
        help="Regression method: 1, 2a, 2b, 3.1, or 3.2.",
    ),

    # Gene filtering
    min_expression: int = typer.Option(
        10, "--min-expression",
        help="Keep genes with max expression > threshold.",
    ),
    variance_quantile: float = typer.Option(
        0.20, "--variance-quantile",
        help="Keep top X%% most variable genes (0.20 = top 20%%).",
    ),

    # Marker gene filtering
    marker_genes_by_tissue_json: Optional[Path] = typer.Option(
        None, "--marker-genes-by-tissue-json",
        exists=True,
        help="JSON mapping tissue ID to marker genes CSV path. "
             "Removes marker genes to reduce segmentation artifacts.",
    ),

    # Bin requirements
    min_neighbor_bins: int = typer.Option(
        4, "--min-neighbor-bins",
        help="Minimum neighbor count bins required (Method 2/3).",
    ),
    min_cells_per_bin: int = typer.Option(
        10, "--min-cells-per-bin",
        help="Minimum cells per bin (Method 2/3).",
    ),

    # Method 2 specific
    n_iterations: int = typer.Option(
        1000, "--n-iterations",
        help="Sampling iterations (Method 2).",
    ),
    sample_size: int = typer.Option(
        10, "--sample-size",
        help="Cells per bin per iteration (Method 2).",
    ),

    # Statistics
    fdr_alpha: float = typer.Option(
        0.05, "--fdr-alpha",
        help="FDR significance threshold.",
    ),

    # Cell type filtering
    primary_cell_types: Optional[str] = typer.Option(
        None, "--primary-cell-types",
        help='Primary cell types as JSON list.\n\n Example: \'["Endothelial","Smooth muscle"]\' (default: all).',
    ),
    neighbor_cell_types: Optional[str] = typer.Option(
        None, "--neighbor-cell-types",
        help='Neighbor cell types as JSON list.\n\n Example: \'["Endothelial","Smooth muscle"]\' (default: all).',
    ),

    # Reproducibility
    rng_seed: Optional[int] = typer.Option(
        None, "--rng-seed",
        help="Random seed for reproducibility (Method 2 sampling).",
    ),

    # Parallelization
    workers: int = typer.Option(
        1, "--workers", "-w",
        help="Parallel workers: positive=exact count, -1=all CPUs, -2=all-1, etc.",
    ),
):
    import matplotlib.pyplot as plt

    out.mkdir(parents=True, exist_ok=True)

    if rng_seed is not None:
        np.random.seed(rng_seed)
        typer.echo(f"Random seed: {rng_seed}")

    if distance_dir is None:
        df_head = pd.read_csv(summary_table_path, nrows=5)
        if "X_space" not in df_head.columns or "Y_space" not in df_head.columns:
            raise typer.BadParameter(
                "When --distance-dir is not provided, the summary table must include "
                "X_space and Y_space columns for coordinate-based distance calculation."
            )

    marker_dict: Optional[Dict[str, List[str]]] = None
    if marker_genes_by_tissue_json is not None:
        marker_dict = json.loads(marker_genes_by_tissue_json.read_text(encoding="utf-8"))
        typer.echo(f"Loaded marker genes for {len(marker_dict)} tissues")

    params = NbrsCountRegressionParams(
        dist_threshold=dist_threshold,
        min_expression=min_expression,
        variance_quantile=variance_quantile,
        min_neighbor_bins=min_neighbor_bins,
        min_cells_per_bin=min_cells_per_bin,
        n_iterations=n_iterations,
        sample_size=sample_size,
        fdr_alpha=fdr_alpha,
    )

    labels = load_labels(summary_table_path)
    genes = get_gene_columns(labels, genes_names_path)

    tissue_list = _parse_json_list(tissue, "--tissue")
    if tissue_list is None:
        tissue_list = get_tissues(labels)

    typer.echo(f"Tissues to analyze: {tissue_list}")
    typer.echo(f"Method: {method}, Distance: {dist_threshold}µm")
    typer.echo(f"Genes: {len(genes)}")

    primary_list = _parse_json_str_list(primary_cell_types, "--primary-cell-types")
    neighbor_list = _parse_json_str_list(neighbor_cell_types, "--neighbor-cell-types")

    if workers > 0:
        n_workers = workers
    else:
        n_workers = max(1, cpu_count() + workers + 1)

    all_results = []
    results_by_pair = {}

    for tissue_id in tissue_list:
        typer.echo(f"\n=== Tissue {tissue_id} ===")

        tissue_genes = genes
        if marker_dict is not None:
            tissue_key = str(tissue_id)
            if tissue_key in marker_dict:
                marker_genes = marker_dict[tissue_key]
                if isinstance(marker_genes, str):
                    marker_df = pd.read_csv(marker_genes)
                    marker_genes = marker_df.iloc[:, 0].tolist()
                tissue_genes = [g for g in genes if g not in marker_genes]
                typer.echo(f"  Filtered {len(genes) - len(tissue_genes)} marker genes")

        distance_matrix = None
        if distance_dir is not None:
            try:
                distance_matrix = load_distance_matrix(distance_dir, tissue_id)
            except FileNotFoundError:
                typer.echo(f"  Warning: No distance matrix for tissue {tissue_id}, using coordinates")

        pairs = get_cell_type_pairs(labels, tissue_id, primary_list, neighbor_list)
        if not pairs:
            typer.echo("  No cell type pairs found, skipping")
            continue

        typer.echo(f"  Cell pairs: {len(pairs)}")

        tissue_workers = min(n_workers, len(pairs))

        if tissue_workers > 1:
            with Pool(
                processes=tissue_workers,
                initializer=_init_worker,
                initargs=(labels, tissue_genes, tissue_id, method, params, distance_matrix, out),
            ) as pool:
                results = pool.map(_process_pair, pairs)
            tissue_results = [r for r in results if r is not None]
        else:
            tissue_results = []
            for target_cell, neighbor_cell in pairs:
                typer.echo(f"  Analyzing: {target_cell} vs {neighbor_cell}...")

                results_df = run_analysis(
                    labels=labels,
                    genes=tissue_genes,
                    tissue=tissue_id,
                    target_cell=target_cell,
                    neighbor_cell=neighbor_cell,
                    method=method,
                    params=params,
                    distance_matrix=distance_matrix,
                )

                if len(results_df) > 0:
                    tissue_results.append(results_df)
                    save_results_csv(
                        results_df, out, method, target_cell, neighbor_cell,
                        dist_threshold, params.fdr_alpha
                    )
                    save_plots(
                        results_df, out, method, target_cell, neighbor_cell,
                        dist_threshold, params.fdr_alpha
                    )

        for df in tissue_results:
            if len(df) == 0:
                continue
            all_results.append(df)
            key = (df["target_cell"].iloc[0], df["neighbor_cell"].iloc[0])
            results_by_pair.setdefault(key, []).append(df)

    if results_by_pair:
        typer.echo("\n=== Creating combined heatmaps ===")
        plots_dir = out / "plots"
        plots_dir.mkdir(parents=True, exist_ok=True)

        for (target_cell, neighbor_cell), dfs in results_by_pair.items():
            combined_df = pd.concat(dfs, ignore_index=True)

            heatmap_fig = create_heatmap(
                combined_df, target_cell, neighbor_cell,
                params.fdr_alpha, tissue_list
            )
            if heatmap_fig is not None:
                base_name = f"method{method}_{target_cell}_vs_{neighbor_cell}_dist{int(dist_threshold)}"
                heatmap_path = plots_dir / f"heatmap_{base_name}_ALL_TISSUES.png"
                heatmap_fig.savefig(heatmap_path, dpi=300, bbox_inches="tight")
                plt.close(heatmap_fig)
                typer.echo(f"  Saved: {heatmap_path.name}")

    if all_results:
        total_results = sum(len(df) for df in all_results)
        total_sig = sum(len(df[df["q_value"] < params.fdr_alpha]) for df in all_results)
        typer.echo(f"\n[DONE] Total: {total_results} results, {total_sig} significant (q < {params.fdr_alpha})")
        typer.echo(f"  output dir: {out}")
    else:
        typer.echo("\n[DONE] No results generated.")
        typer.echo(f"  output dir: {out}")
