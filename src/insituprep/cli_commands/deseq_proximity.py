from __future__ import annotations

import json
from pathlib import Path
from typing import Optional, Dict, List
import pandas as pd

import typer

from insituprep.analyses.deseq_proximity import run_deseq_pairs, DeseqRunParams

app = typer.Typer(
    help=(
        "Differential expression analysis for proximity-defined cell states.\n"
        "Tests whether genes are differentially expressed in primary-type cells proximal vs distant\n"
        "to neighbor-type cells using DESeq2 (PyDESeq2) with permutation-based validation."
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


@app.command("run")
def run(
    summary_table_path: Path = typer.Option(..., "--summary-table-path",
                                     help="CSV with per-cell metadata including cell centroid coordinates (µm) and gene expression.\n\n"
                                            "Required columns: 'Var1' (cell IDs), 'tissue', 'cell_type'.\n\n"
                                            "The table should also include gene expression columns "
                                            "(one column per gene, raw numeric count values).\n\n"
                                            "Proximity source:\n\n"
                                            "  - If --distance-dir is provided: proximity is computed from the per-tissue distance matrices.\n\n"
                                            "  - If --distance-dir is NOT provided: the summary table MUST include centroid coordinates:\n\n"
                                            "      * Required: 'X_space', 'Y_space'\n\n"
                                            "      * Optional: 'Z_space' (if present, 3D distances are used; otherwise 2D).\n\n"

                                     ),

    distance_matrix_path: Optional[Path] = typer.Option(
    None,
    "--distance-matrix-path",
        help="Optional cell-cell distance matrix CSV for the selected tissue "
            "('distance_matrix_{tissue}.csv'), containing pairwise cell–cell distances (µm). "
            "Rows and columns correspond to cell IDs.\n\n "
            "If omitted, proximity is computed from cell coordinates in the summary table "
            "(X_space, Y_space, and optionally Z_space). "
    ),
    
    genes_names_path: Path = typer.Option(..., "--genes-names-path", help="Text file containing one gene name per line. "
                                                                    "Used to select a subset of genes and extract the corresponding "
                                                                    "gene expression columns from the 'summary_table' CSV."),

    tissue: str = typer.Option(..., "--tissue", help='Single tissue ID as string. Example: "100". One tissue per run.'),
    out: Path = typer.Option(..., "--out", help="Output directory"),
    dist_threshold: float = typer.Option(1.0, "--dist-threshold", help="Distance threshold (µm) for proximity labeling."),
    n_perm: int = typer.Option(1000, "--n-perm", help="Number of permutations/iterations"),
    primary: Optional[str] = typer.Option(
        None,
        "--primary",
        help='JSON list of primary cell types. Example: \'["T_CD3","B"]\'.\n\n '
            'If omitted: use all cell types.',
    ),
    neighbor: Optional[str] = typer.Option(
        None,
        "--neighbor",
        help='JSON list of neighbor cell types. Example: \'["Epithelial","Endothelial"]\'.\n\n '
            'If omitted: use all cell types.',
    ),
    marker_genes_by_tissue_json: Optional[Path] = typer.Option(
        None,
        "--marker-genes-by-tissue-json",
        help="Optional JSON mapping tissue->marker genes csv path.\n\n"
        "Marker filtering removes genes that appear as markers for the neighbor cell type in order to reduce potential segmentation artifacts. "),

    perm_print_every: int = typer.Option(50, "--perm-print-every", help="Print permutation progress every N"),
    quiet_permutations: bool = typer.Option(
        True,
        "--quiet-permutations/--no-quiet-permutations",
        help="Silence PyDESeq2 output during permutations",
    ),
):
    """
    Run DESeq2 proximity analysis for a single tissue (and optionally a subset of primary/neighbor cell types).
    """
    out.mkdir(parents=True, exist_ok=True)

    marker_dict: Optional[Dict[str, str]] = None
    if marker_genes_by_tissue_json is not None:
        marker_dict = json.loads(marker_genes_by_tissue_json.read_text(encoding="utf-8"))

    primary_list = _parse_json_list_of_strings(primary, "--primary")
    neighbor_list = _parse_json_list_of_strings(neighbor, "--neighbor")
    
    if distance_matrix_path is None:
        # Validate that summary table has coordinate columns (run_deseq_pairs will also check, but CLI can fail fast)
        df_head = pd.read_csv(summary_table_path, nrows=5)  # quick check
        if "X_space" not in df_head.columns or "Y_space" not in df_head.columns:
            raise typer.BadParameter(
                "When --distance-matrix-path is not provided, summary table CSV must include X_space and Y_space columns "
                "(and optionally Z_space) for coordinate-based proximity."
            )

    params = DeseqRunParams(
        tissue=str(tissue),
        dist_threshold=float(dist_threshold),
        n_perm=int(n_perm),
    )

    res = run_deseq_pairs(
        summary_table_path=summary_table_path,
        distance_matrix_path=distance_matrix_path,
        genes_names_path=genes_names_path,
        out_dir=out,
        params=params,
        marker_genes_by_tissue=marker_dict,
        primary_cell_types=primary_list,
        neighbor_cell_types=neighbor_list,
        quiet_permutations=quiet_permutations,
        perm_print_every=int(perm_print_every),
    )

    typer.echo(f"[DONE] summary: {res['summary_csv']}")
    typer.echo(f"  pairs done: {res['n_pairs_done']}")
