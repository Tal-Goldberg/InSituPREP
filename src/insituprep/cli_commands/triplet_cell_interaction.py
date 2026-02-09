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
        help="CSV with per-cell metadata including coordinates (X_space, Y_space, Z_space) and gene expression. "
            "Required columns: 'Var1' (cell IDs), 'tissue', 'cell_type'. "
            "The table must also include columns: 'X_space', 'Y_space' (cell centroid coordinates in µm). "
            "The table should also include gene expression columns. "
            "(one column per gene, raw numeric counts values). "
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
    
    out: Path = typer.Option(..., "--out", help="Output directory."),
    tissue: Optional[str] = typer.Option(
        None,
        "--tissue",
        help="Tissue identifier (string). "
        "If provided, the analysis is restricted to cells from this tissue only. "
        "If not provided, the analysis is run on the full summary table "
        "(i.e. all tissues together). "
        "Note: running without --tissue may mix cells from different tissues."
    ),
    attach: float = typer.Option(10.0, "--attach", help="Radius (µm) for attaching neighbors."),
    alpha: float = typer.Option(0.05, "--alpha", help="p-value threshold on DESeq2 results."),
    n_cpus: int = typer.Option(8, "--n-cpus", help="Number of CPUs for PyDESeq2."),
    tissue_col: str = typer.Option("tissue", "--tissue-col", help="Tissue column name."),
    celltype_col: str = typer.Option("cell_type", "--celltype-col", help="Cell type column name."),
    cell_id_col: str = typer.Option("Var1", "--cell-id-col", help="Cell ID column name."),
):
    out.mkdir(parents=True, exist_ok=True)

    genes = _read_genes_list(genes_names_path)

    df = pd.read_csv(
        summary_table_path,
        dtype={cell_id_col: str, tissue_col: str},
    )

    params = TripletParams(
        attach=attach,
        alpha=alpha,
        n_cpus=n_cpus,
    )

    res = run_triplet_cell_interaction(
        labels_df=df,
        genes=genes,
        params=params,
        tissue=tissue,
        tissue_col=tissue_col,
        celltype_col=celltype_col,
        cell_id_col=cell_id_col,
    )

    suffix = f"tissue_{tissue}" if tissue is not None else "all_tissues"
    out_path = out / f"triplets_deseq_{suffix}.csv"
    res.to_csv(out_path)

    typer.echo(f"[OK] Wrote: {out_path}")
