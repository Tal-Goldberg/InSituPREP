from __future__ import annotations
import json
from pathlib import Path
from typing import Optional, List
import typer

# filepath: expression_as_function_of_distance.py



from insituprep.analyses.expression_as_function_of_distance import (
    run_expression_as_function_of_distance,
    ExpressionDistanceRunParams,
)

app = typer.Typer(
    help=(
        "Expression as a Function of Distance:\n"
        "Analyzes how gene expression in primary-type cells correlates with distance to neighbor-type cells.\n"
        "Including linear regression, permutation testing, FDR correction, and Gaussian R² filtering."
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


@app.command("run",
    help=(
        "Expression as a Function of Distance:\n"
        "Analyzes how gene expression in primary-type cells correlates with distance to neighbor-type cells.\n"
        "Including linear regression, permutation testing, FDR correction, and Gaussian R² filtering."

))
def run(
    summary_table_path: Path = typer.Option(
        ...,
        "--summary-table-path",
        exists=True,
        help="CSV with per-cell metadata including coordinates (X_space, Y_space, Z_space) and gene expression.\n\n "
            "Required columns: 'Var1' (cell IDs), 'tissue', 'cell_type'.\n\n "
            "The table must also include columns: 'X_space', 'Y_space' (cell centroid coordinates in µm).\n\n "
            "The table should also include gene expression columns. "
            "(one column per gene, raw numeric counts values).\n\n "
            "Optional column 'Z_space' for 3D coordinates (defaults to 2D if absent). "
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
    sigma_param: float = typer.Option(
        5.0,
        "--sigma-param",
        help="Gaussian smoothing parameter for R² filtering."
        "Sigma value for Gaussian smoothing applied to expression data during regression analysis to reduce noise.",
    ),
    r2_threshold: float = typer.Option(
        0.3,
        "--r2-threshold",
        help="R² threshold for Gaussian-smoothed regression filtering."
        "Minimum R-squared threshold for retaining genes in the Gaussian-filtered regression; genes below this are excluded.",
    ),
    maximum_distance: float = typer.Option(
        145.0,
        "--maximum-distance",
        help="Maximum distance (µm) to consider. Set to 0 for unlimited."
        "Maximum Euclidean distance (in µm) between primary and neighbor cells to consider for interactions; set to 0 to include all distances. 145 is approximately the field of view size in ExSeq, equivalent to a square of ~100x100x25 µm.",
    ),
    leg_threshold: float = typer.Option(
        10.0,
        "--leg-threshold",
        help="Low Expression Gene (LEG) threshold. Genes with 98th percentile below this are filtered."
        "Minimum expression level required for a gene in a cell type to be considered for analysis (Low Expression Genes threshold).",
    ),
    cell_count_threshold: int = typer.Option(
        20,
        "--cell-count-threshold",
        help="Minimum number of cells expressing a gene to retain it."
        "Minimum number of cells required per cell type in a tissue to include that cell type in the analysis.",
    ),
    num_iterations: int = typer.Option(
        100,
        "--num-iterations",
        help="Number of permutations for statistical significance testing."
        "Number of iterations for the permutation test to assess statistical significance of regression results.",
    ),
    fdr_thresh: float = typer.Option(
        0.05,
        "--fdr-thresh",
        help="FDR cutoff for significant genes (BH correction)."
        "False Discovery Rate (FDR) threshold for correcting p-values across all tests; genes with FDR above this are filtered out.",
    ),
    tissue: Optional[str] = typer.Option(
        None,
        "--tissue",
        help='JSON list of tissue IDs as strings.\n\n Example: \'["100","313"]\'.\n\n If omitted: run all tissues.',
    ),
    primary: Optional[str] = typer.Option(
        None,
        "--primary",
        help='JSON list of primary cell types.\n\n Example: \'["Endothelial","T_CD3"]\'.\n\n If omitted: use all cell types.',
    ),
    neighbor: Optional[str] = typer.Option(
        None,
        "--neighbor",
        help='JSON list of neighbor cell types.\n\n Example: \'["Epithelial","B"]\'.\n\n If omitted: use all cell types.',
    ),
):
    """
    Run expression-as-function-of-distance analysis.
    """
    out.mkdir(parents=True, exist_ok=True)

    tissue_list = _parse_json_list_of_strings(tissue, "--tissue")
    primary_list = _parse_json_list_of_strings(primary, "--primary") or []
    neighbor_list = _parse_json_list_of_strings(neighbor, "--neighbor") or []

    params = ExpressionDistanceRunParams(
        sigma_param=float(sigma_param),
        r2_threshold=float(r2_threshold),
        maximum_distance=float(maximum_distance),
        leg_threshold=float(leg_threshold),
        cell_count_threshold=int(cell_count_threshold),
        num_iterations=int(num_iterations),
        fdr_thresh=float(fdr_thresh),
        tissue_ids=tissue_list,
        primary_cell_types=primary_list,
        neighbor_cell_types=neighbor_list,
    )

    res = run_expression_as_function_of_distance(
        summary_table_path=summary_table_path,
        genes_names_path=genes_names_path,
        output_path=out,
        params=params,
    )

    typer.echo(f"[DONE] Outputs written to: {out}")
    if isinstance(res, dict):
        if "final_results_csv" in res:
            typer.echo(f"  final results: {res['final_results_csv']}")
        if "n_results" in res:
            typer.echo(f"  significant genes found: {res['n_results']}")
        if "output_directory" in res:
            typer.echo(f"  output directory: {res['output_directory']}")


if __name__ == "__main__":
    app()
