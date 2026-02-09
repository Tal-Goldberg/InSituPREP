from __future__ import annotations

import json
from pathlib import Path
from typing import List, Optional

import typer

from insituprep.analyses.comparison_across_platforms_and_patients import (
    run_comparison_across_platforms_and_patients,
    ComparisonAcrossPlatformsRunParams,
)

app = typer.Typer(
    help=(
        "Compare proximity-based DE signatures across platforms and patients/tissues.\n"
        "Aggregates per-tissue-and-cell types pair DESeq outputs and evaluates cross-platform and patients status concordance\n"
        "using permutation tests."
    )
)


def _parse_json_list_of_strings(s: str, opt_name: str) -> List[str]:
    try:
        x = json.loads(s)
    except Exception as e:
        raise typer.BadParameter(f"{opt_name} must be a valid JSON list of strings.") from e
    if not isinstance(x, list) or not all(isinstance(v, str) for v in x):
        raise typer.BadParameter(f"{opt_name} must be a JSON list of strings.")
    return x


@app.command("run", help ="Compare proximity-based DE signatures across platforms and patients/tissues.\n"
        "Aggregates per-tissue-and-cell types pair DESeq outputs and evaluates cross-platform and patients status concordance\n"
        "using permutation tests.")
def run(
    results_dir: Path = typer.Option(
        ...,
        "--deseq-results-dir",
        exists=True,
        file_okay=False,
        dir_okay=True,
        readable=True,
        help=(
            "Folder with DESeq results CSVs named like: "
            "final_genes_{tissue}_{primary_cell_type}_vs_{neighbor_cell_type}_LFC_*.csv"
        ),
    ),
    out: Path = typer.Option(
        ...,
        "--out",
        file_okay=False,
        dir_okay=True,
        help="Output directory.",
    ),
    tissues: str = typer.Option(
        ...,
        "--tissues",
        help='JSON list of tissues (strings). Example: \'["364","330","MERFISH_330"]\'',
    ),
    non_tumor_cell_types: str = typer.Option(
        ...,
        "--non-tumor-cell-types",
        help='JSON list of non-tumor cell types. Example: \'["B","Endothelial","T_CD3"]\'',
    ),
    tumor_cell_type: str = typer.Option(
        ...,
        "--tumor-cell-type",
        help="Tumor cell type name (as used in the final_genes files).",
    ),
    hr_positive_tissues: str = typer.Option(
        '["982","880","330","783","MERFISH_982","MERFISH_880"]',
        "--hr-positive-tissues",
        help="JSON list of tissues considered HR+/HER2-.",
    ),
    n_perm_pca: int = typer.Option(
        10000,
        "--n-perm-pca",
        help="Permutations for PC explained variance p-values.",
    ),
    n_perm_dist: int = typer.Option(
        10000,
        "--n-perm-dist",
        help="Permutations for platform distance and receptor proximity p-values.",
    ),
    rng_seed: Optional[int] = typer.Option(
        None,
        "--rng-seed",
        help="Random seed for reproducibility. If not provided, the run is stochastic.",
    ),
    platform_prefix: str = typer.Option(
        "MERFISH_",
        "--platform-prefix",
        help="Prefix marking the second platform tissue IDs.",
    ),
    verbose: bool = typer.Option(
        True,
        "--verbose/--no-verbose",
        help="Print progress messages during the run.",
    ),
    perm_progress_every: int = typer.Option(
        1000,
        "--perm-progress-every",
        help="When verbose: print permutation progress every N permutations (0 disables).",
    ),
):
    """
    Run cross-platform / cross-patient comparison of proximity-based DE signatures.
    """
    tissues_list = _parse_json_list_of_strings(tissues, "--tissues")
    non_tumor_ct_list = _parse_json_list_of_strings(non_tumor_cell_types, "--non-tumor-cell-types")
    hr_pos_list = _parse_json_list_of_strings(hr_positive_tissues, "--hr-positive-tissues")

    out.mkdir(parents=True, exist_ok=True)

    params = ComparisonAcrossPlatformsRunParams(
        tissues=tissues_list,
        non_tumor_cell_types=non_tumor_ct_list,
        tumor_cell_type=tumor_cell_type,
        platform_prefix=platform_prefix,
        hr_positive_tissues=hr_pos_list,
        n_perm_pca=int(n_perm_pca),
        n_perm_dist=int(n_perm_dist),
        verbose=bool(verbose),
        perm_progress_every=int(perm_progress_every),
        final_genes_prefix_template=(
            "final_genes_{tissue}_{primary_cell_type}_vs_{neighbor_cell_type}_LFC"
        ),
        gene_column="gene",
        pval_column="pvalue",
    )

    res = run_comparison_across_platforms_and_patients(
        results_dir=results_dir,
        out_dir=out,
        params=params,
        rng_seed=rng_seed,
    )

    typer.echo(f"[DONE] Summary CSV: {res['summary_csv']}")
    typer.echo(f"Rows: {res['n_rows']}")
