from __future__ import annotations

import typer


from insituprep.cli_commands.deseq_proximity import app as deseq_proximity_app
from insituprep.cli_commands.dispersion_pca import app as dispersion_pca_app
from insituprep.cli_commands.comparison_across_platforms import app as comparison_across_platforms_app
from insituprep.cli_commands.rna_velocity_app import app as rna_velocity_app
from insituprep.cli_commands.bacteria_objects import app as bacteria_objects_app
from insituprep.cli_commands.expression_as_function_of_distance import app as exprdist_app
from insituprep.cli_commands.nbrs_count_regression_cli import app as nbrs_app
from insituprep.cli_commands.triplet_cell_interaction import app as triplet_app

app = typer.Typer(add_completion=True, help="InSituPREP command-line interface")
app.add_typer(deseq_proximity_app, name="deseq-proximity")
app.add_typer(dispersion_pca_app, name="dispersion-pca")
app.add_typer(comparison_across_platforms_app, name="platform-comparison")
app.add_typer(rna_velocity_app, name="rna-velocity")
app.add_typer(bacteria_objects_app, name="bacteria-objects")
app.add_typer(exprdist_app, name="expression-distance")
app.add_typer(nbrs_app, name="nbrs-count-regression")
app.add_typer(triplet_app, name="triplet-cell-interaction")

@app.callback()
def _root() -> None:
    """InSituPREP CLI."""
    pass


def main() -> None:
    app()


if __name__ == "__main__":
    main()
