# src/insituprep/cli_commands/bacteria_objects.py
from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from insituprep.analyses.bacteria_objects import run_bacteria_objects, BacteriaObjectsParams

app = typer.Typer(
    help=(
        "Detect statistically significant spatial transcript objects for a target gene (e.g., bacterial transcripts) "
        "using DBSCAN + permutation testing, and optionally run DESeq2 per cell type by "
        "comparing cells that are near to or far from significant transcript objects. "
     )
)


@app.command("run", help="Detect statistically significant spatial transcript objects for a target gene (e.g., bacterial transcripts) "
        "using DBSCAN + permutation testing, and optionally run DESeq2 per cell type by "
        "comparing cells that are near to or far from significant transcript objects. ")
def run(
    # ----------------------
    # Required I/O
    # ----------------------
    transcripts_path: Path = typer.Option(
        ...,
        "--transcripts-path",
        exists=True,
        readable=True,
        help=(
            "Path to a transcripts CSV (one row per detected transcript). "
            "Must include: a gene name column, global transcript coordinates (e.g., in microns), "
            "and a cell_id column (string; empty/NA if the transcript is not assigned to a cell)."
        ),
    ),
    out_dir: Path = typer.Option(
        ...,
        "--out-dir",
        exists=False,
        help=("Output directory. All results of the analysis will be written here:\n "
              "- transcripts_with_objects.csv\n "
              "- eps_sweep_summary.csv\n "
              "- deseq2_results/ (per-cell-type DESeq2 outputs, if requested) "
            ),
    ),
    target_gene: str = typer.Option(
        ...,
        "--target-gene",
        help='Target gene to detect objects for (e.g., "BACTERIA").',
    ),

    # ----------------------
    # Transcripts columns
    # ----------------------
    gene_col: str = typer.Option(
        "gene_name",
        "--gene-col",
        help="Column name in transcripts CSV containing gene names.",
    ),
    tx_x_col: str = typer.Option(
        "x_micron",
        "--tx-x-col",
        help="Column name in transcripts CSV for global X coordinate (e.g., microns).",
    ),
    tx_y_col: str = typer.Option(
        "y_micron",
        "--tx-y-col",
        help="Column name in transcripts CSV for global Y coordinate (e.g., microns).",
    ),

    tx_z_col: Optional[str] = typer.Option(
        None,
        "--tx-z-col",
        help=(
            "Column name in transcripts CSV for global Z coordinate (e.g., microns). "
            "If not provided, the analysis runs in 2D (X,Y only)."
        ),
    ),


    tx_cell_id_col: str = typer.Option(
        "cell_id",
        "--tx-cell-id-col",
        help="Column name in transcripts CSV containing the global cell_id (string). Empty/NA means not inside a cell.",
    ),

    # ----------------------
    # DBSCAN / permutation
    # ----------------------
    eps_min: int = typer.Option(1, "--eps-min", min=1, help="Minimum eps value in the sweep."),
    eps_max: int = typer.Option(150, "--eps-max", min=1, help="Maximum eps value in the sweep."),
    min_samples: int = typer.Option(3, "--min-samples", min=1, help="DBSCAN min_samples."),
    n_perms: int = typer.Option(10_000, "--n-perms", min=1, help="Number of permutations per cluster for p-values."),
    alpha: float = typer.Option(0.05, "--alpha", min=0.0, max=1.0, help="Cluster significance threshold (p-value)."),
    eps_final: str = typer.Option(
        "auto",
        "--eps-final",
        help=(
            'Final eps used for the final assignment. Use an integer (e.g., "29") '
            'or "auto" to choose the eps that maximizes global F1. '
            '"auto" selects the eps value that maximizes the global F1 score across the eps sweep.'
        ),
    ),

    seed: Optional[int] = typer.Option(
        None,
        "--seed",
        help=(
            "Random seed for reproducibility (in permutations). "
            "If not provided, the analysis runs with a random seed."
        ),
    ),

    # ----------------------
    # resume / reuse existing object-detection outputs
    # ----------------------
    reuse_existing_objects: bool = typer.Option(
        False,
        "--reuse-existing-objects",
        help=(
            "If transcripts_with_objects.csv and eps_sweep_summary.csv already exist in --out-dir, "
            "load them and skip DBSCAN+permutations. Useful for resuming after DESeq2 errors."
        ),
    ),

    # ----------------------
    # Optional DESeq2 step
    # ----------------------
    summary_table_path: Optional[Path] = typer.Option(
        None,
        "--summary-table-path",
        exists=True,
        readable=True,
        help="Optional CSV with per-cell metadata including coordinates (X_space, Y_space, Z_space) and gene expression. "
            "Required columns: 'Var1' (cell IDs), 'tissue', 'cell_type'. "
            "The table must also include columns: 'X_space', 'Y_space' (cell centroid coordinates in µm). "
            "The table should also include gene expression columns. "
            "(one column per gene, raw numeric counts values). "
            "Optional column 'Z_space' for 3D coordinates (defaults to 2D if absent). "

    ),
    genes_names_path: Optional[Path] = typer.Option(
        None,
        "--genes-names-path",
        exists=True,
        readable=True,
        help=(
            "Optional text file listing gene names (one per line). "
            "When running DESeq2, ONLY these genes will be used as expression columns, "
            "and the code verifies they all exist in the 'summary_table' CSV."
        ),
    ),

    celltype_col: str = typer.Option(
        "cell_type",
        "--celltype-col",
        help="Cell type column in cells summary.",
    ),
    cell_x_col: str = typer.Option(
        "X_space",
        "--cell-x-col",
        help="Cell X coordinate column (same coordinate frame/units as transcripts).",
    ),
    cell_y_col: str = typer.Option(
        "Y_space",
        "--cell-y-col",
        help="Cell Y coordinate column (same coordinate frame/units as transcripts).",
    ),
    cell_z_col: Optional[str] = typer.Option(
        None,
        "--cell-z-col",
        help=(
            "Cell Z coordinate column (same coordinate frame/units as transcripts). "
            "If not provided, DESeq2 grouping uses 2D distances (X,Y only)."
        ),
    ),

    dist_threshold: float = typer.Option(
        10.0,
        "--dist-threshold",
        min=0.0,
        help=(
            "Distance threshold used to assign cells to near_sig vs far_sig. "
            "near_sig: min distance to ANY significant target transcript <= threshold. "
            "far_sig:  min distance > threshold."
        ),
    ),
    min_cells_per_group: int = typer.Option(
        20,
        "--min-cells-per-group",
        min=1,
        help="Minimum cells required in EACH group (near_sig & far_sig) per cell type to run DESeq2.",
    ),

):
    # If user requested DESeq2, enforce genes_names_path
    if summary_table_path is not None and genes_names_path is None:
        raise typer.BadParameter("--genes-names-path is required when --summary-table-path is provided.")

    eps_final_str = eps_final.strip().lower()
    eps_final_val = None if eps_final_str == "auto" else int(eps_final_str)


    params = BacteriaObjectsParams(
        input_csv=transcripts_path,
        out_dir=out_dir,

        target_gene=target_gene,
        gene_col=gene_col,
        tx_x_col=tx_x_col,
        tx_y_col=tx_y_col,
        tx_z_col=tx_z_col,
        tx_cell_id_col=tx_cell_id_col,

        eps_min=eps_min,
        eps_max=eps_max,
        min_samples=min_samples,
        n_perms=n_perms,
        alpha=alpha,
        eps_final=eps_final_val,
        seed=seed,

        reuse_existing_objects=reuse_existing_objects,

        summary_table_path=summary_table_path,
        celltype_col=celltype_col,
        cell_x_col=cell_x_col,
        cell_y_col=cell_y_col,
        cell_z_col=cell_z_col,
        dist_threshold=dist_threshold,
        min_cells_per_group=min_cells_per_group,

        genes_names_path=genes_names_path,
    )

    run_bacteria_objects(params)
