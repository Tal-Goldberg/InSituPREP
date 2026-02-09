from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, List

import typer

from insituprep.analyses.dispersion_in_pca_space import (
    run_dispersion_in_pca_space,
    DispersionPCARunParams,
)

app = typer.Typer(
    help=(
        "Dispersion in PCA space:\n"
        "Tests whether primary-type cells proximal to neighbor-type occupy a distinct region in PCA space\n"
        "vs non-proximal primary cells. Includes permutation testing + BH-FDR and sensitivity analyses."
    )
)


def _parse_json_list_of_strings(s: Optional[str], opt_name: str) -> Optional[List[str]]:
    if s is None:
        return None
    try:
        x = json.loads(s)
    except Exception as e:
        raise typer.BadParameter(f"{opt_name} must be a valid JSON list of strings.") from e
    if not isinstance(x, list) or not all(isinstance(v, str) for v in x):
        raise typer.BadParameter(f"{opt_name} must be a JSON list of strings.")
    return x


@app.command("run", help= "Dispersion in PCA space:\n"
        "Tests whether primary-type cells proximal to neighbor-type occupy a distinct region in PCA space "
        "vs non-proximal primary cells. Includes permutation testing + BH-FDR and sensitivity analyses.")
def run(
    summary_table_path: Path = typer.Option(
        ...,
        "--summary-table-path",
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
    distance_dir: Optional[Path] = typer.Option(
        None,
        "--distance-dir",
        exists=True,
        file_okay=False,
        help="Optional directory with per-tissue cell-cell distance matrices "
            "('distance_matrix_{tissue}.csv'), containing pairwise cell–cell distances (µm). "
            "Rows and columns correspond to cell IDs. "
            "If omitted, proximity is computed from cell coordinates in the summary table file "
            "(X_space, Y_space, and optionally Z_space). "
    ),
    genes_names_path: Path = typer.Option(
        ...,
        "--genes-names-path",
        exists=True,
        help="Text file containing one gene name per line. "
            "Used to select a subset of genes and extract the corresponding "
            "gene expression columns from the 'summary_table' CSV."
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        help="Output directory.",
    ),
    dist_threshold: float = typer.Option(
        1.0,
        "--dist-threshold",
        help="Distance threshold (µm) for proximity labeling between primary-type cell and neighbor-type cell",
    ),
    n_perm: int = typer.Option(
        10000,
        "--n-perm",
        help="Number of permutations.",
    ),
    fdr_alpha: float = typer.Option(
        0.05,
        "--fdr-alpha",
        help="FDR cutoff for significant pairs (BH).",
    ),
    tissue: Optional[str] = typer.Option(
        None,
        "--tissue-ids",
        help='JSON list of tissue IDs as strings. Example: \'["100","313"]\'. If omitted: run all tissues.',
    ),
    primary: Optional[str] = typer.Option(
        None,
        "--primary",
        help='JSON list of primary cell types. Example: \'["Endothelial","T_CD3"]\'.'
        'If omitted, all cell types are treated as primary.'
    ),
    neighbor: Optional[str] = typer.Option(
        None,
        "--neighbor",
        help='JSON list of neighbor cell types. Example: \'["Epithelial","B"]\'.'
        'If omitted, all cell types are treated as neighbor.'
    ),
    marker_genes_by_tissue_json: Optional[Path] = typer.Option(
        None,
        "--marker-genes-by-tissue-json",
        exists=True,
        help="Optional JSON mapping tissue->marker CSV path (more info in readme file)."
        "Marker filtering removes genes that appear as markers for the neighbor cell type in order to reduce potential segmentation artifacts. "
    ),
    skip_if_no_markers: bool = typer.Option(
        True,
        "--skip-if-no-markers/--no-skip-if-no-markers",
        help="If True, skip marker-removal stage for tissues without markers.",
    ),
    tsne_perplexity: int = typer.Option(
        30,
        "--tsne-perplexity",
        help="t-SNE perplexity (sensitivity test).",
    ),
    tsne_random_state: Optional[int] = typer.Option(
        None,
        "--tsne-random-state",
        help="t-SNE random_state (None for stochastic).",
    ),
    plot_n_cols: int = typer.Option(
        4,
        "--plot-n-cols",
        help="Number of columns in the significant-pairs PCA panel plot.",
    ),
    plot_max_panels: Optional[int] = typer.Option(
        None,
        "--plot-max-panels",
        help="Optional max number of panels to plot.",
    ),
    rng_seed: Optional[int] = typer.Option(
        None,
        "--rng-seed",
        help="Random seed for reproducibility.",
    ),
):
    """
    Run dispersion-in-PCA analysis.
    """
    out.mkdir(parents=True, exist_ok=True)

    tissue_list = _parse_json_list_of_strings(tissue, "--tissue-ids")
    primary_list = _parse_json_list_of_strings(primary, "--primary")
    neighbor_list = _parse_json_list_of_strings(neighbor, "--neighbor")

    params = DispersionPCARunParams(
        dist_threshold=float(dist_threshold),
        n_perm=int(n_perm),
        fdr_alpha=float(fdr_alpha),
        tissue_ids=tissue_list,
        primary_cell_types=primary_list,
        neighbor_cell_types=neighbor_list,
        tsne_perplexity=int(tsne_perplexity),
        tsne_random_state=int(tsne_random_state) if tsne_random_state is not None else None,
        plot_n_cols=int(plot_n_cols),
        plot_max_panels=int(plot_max_panels) if plot_max_panels is not None else None,
        skip_if_no_markers=bool(skip_if_no_markers),
    )

    res = run_dispersion_in_pca_space(
        summary_table_path=summary_table_path,
        genes_names_path=genes_names_path,
        distance_dir=distance_dir,
        out_dir=out,
        params=params,
        marker_genes_by_tissue_json=marker_genes_by_tissue_json,
        rng_seed=rng_seed,
    )

    typer.echo(f"[DONE] Outputs written to: {out}")
    if isinstance(res, dict):
        if "summary_csv" in res:
            typer.echo(f"  summary: {res['summary_csv']}")
        if "sig_pairs_csv" in res:
            typer.echo(f"  significant pairs: {res['sig_pairs_csv']}")
        if "sig_plot_path" in res:
            typer.echo(f"  significant plot: {res['sig_plot_path']}")
