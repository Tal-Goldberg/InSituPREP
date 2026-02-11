from __future__ import annotations

from pathlib import Path
from typing import Optional, List

import typer
import pandas as pd

from insituprep.analyses.triplet_cell_interaction import (
    TripletParams,
    run_triplet_cell_interaction,
)

app = typer.Typer(add_completion=True, help="Triplet cell interaction analysis (DESeq2 over triplet neighborhoods).")


def _read_genes_list(path: Path) -> List[str]:
    genes = pd.Series(path.read_text().splitlines()).astype(str)
    genes = genes.str.strip()
    genes = genes[genes != ""].tolist()
    return genes


@app.command("run")
def run(
    summary_table_path: Path = typer.Option(
        ...,
        "--summary-table-path",
        exists=True,
        help=(
            "CSV with per-cell metadata including coordinates (X_space, Y_space, Z_space) and gene expression. "
            "Required columns: 'Var1' (cell IDs), 'cell_type'. "
            "Optional column: 'tissue' (string tissue identifier). "
            "The table must also include columns: 'X_space', 'Y_space' (cell centroid coordinates in µm). "
            "Optional column 'Z_space' for 3D coordinates (defaults to 2D if absent). "
            "The table should also include gene expression columns (one column per gene, raw numeric counts values)."
        ),
    ),
    genes_names_path: Path = typer.Option(
        ...,
        "--genes-names-path",
        exists=True,
        help=(
            "Text file containing one gene name per line. "
            "Used to select a subset of genes and extract the corresponding gene-expression columns."
        ),
    ),
    out: Path = typer.Option(..., "--out", help="Output directory."),
    tissue: Optional[str] = typer.Option(
        None,
        "--tissue",
        help=(
            "Tissue identifier (string). "
            "If provided, the analysis is restricted to this tissue only (requires the tissue column to exist). "
            "If not provided and tissue column exists, the analysis is run per tissue and results are concatenated. "
            "If the tissue column does not exist, the analysis runs once on the whole table."
        ),
    ),
    pv_threshold: float = typer.Option(
        0.05,
        "--pv-threshold",
        help=(
            "Adjusted p-value threshold (padj) used for filtering results when a tissue column exists. "
            "Matches the updated script behavior (padj < pv_threshold)."
        ),
    ),
    attach: float = typer.Option(10.0, "--attach", help="Radius (µm) for attaching neighbors."),
    n_cpus: int = typer.Option(8, "--n-cpus", help="Number of CPUs for PyDESeq2."),
    refit_cooks: bool = typer.Option(True, "--refit-cooks/--no-refit-cooks", help="PyDESeq2 refit_cooks flag."),
    tissue_col: str = typer.Option("tissue", "--tissue-col", help="Tissue column name."),
    celltype_col: str = typer.Option("cell_type", "--celltype-col", help="Cell type column name."),
    cell_id_col: str = typer.Option("Var1", "--cell-id-col", help="Cell ID column name."),
):
    out.mkdir(parents=True, exist_ok=True)

    genes = _read_genes_list(genes_names_path)
    if len(genes) == 0:
        raise typer.BadParameter("genes-names file is empty (no gene names found).")


    # Keep dtype for Var1; for tissue dtype, only enforce if column exists after read.
    df = pd.read_csv(summary_table_path, dtype={cell_id_col: str})

    required_cols = {cell_id_col, celltype_col, "X_space", "Y_space"}
    missing_cols = [c for c in required_cols if c not in df.columns]
    if missing_cols:
        raise typer.BadParameter(
            f"Missing required columns in summary table: {missing_cols}"
        )

    if tissue is not None and tissue_col not in df.columns:
        raise typer.BadParameter(
            f"--tissue was provided but tissue column '{tissue_col}' was not found in the summary table."
        )


    if tissue_col in df.columns:
        df[tissue_col] = df[tissue_col].astype(str)

    params = TripletParams(
        attach=attach,
        n_cpus=n_cpus,
        refit_cooks=refit_cooks,
    )

    try:
        res = run_triplet_cell_interaction(
            labels_df=df,
            genes=genes,
            params=params,
            pv_threshold=pv_threshold,
            tissue=tissue,
            tissue_col=tissue_col,
            celltype_col=celltype_col,
            cell_id_col=cell_id_col,
        )
    except ValueError as e:
        raise typer.BadParameter(str(e))

    suffix = f"tissue_{tissue}" if tissue is not None else "all"
    out_path = out / f"triplets_deseq_{suffix}.csv"
    res.to_csv(out_path, index=True)

    typer.echo(f"[OK] Wrote: {out_path}")
