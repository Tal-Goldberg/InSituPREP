# InSituPREP

**InSituPREP** is a modular Python package for the analysis of **spatial transcriptomics data at single-cell resolution**.  
The package provides a collection of complementary, proximity-aware analyses designed to decompose and interpret how **spatial organization shapes gene expression, cellular states, and cell–cell interactions**.

Rather than relying on a single analytical paradigm, InSituPREP is built as a **multi-layered analytical toolkit**, where each analysis captures a different aspect of spatial structure in tissue.

---

## Conceptual overview

Spatial transcriptomics data contains multiple intertwined signals:

- Continuous spatial gradients
- Binary proximity effects (near vs far)
- Neighborhood composition
- Multivariate expression structure
- Dynamic transcriptional states

InSituPREP addresses these through **separate but interoperable analyses**, each answering a distinct biological question:

- Does gene expression change smoothly with distance to another cell type?
- Does spatial proximity contributes to within cell type gene expression variability? 
- Are transcriptional dynamics altered by spatial context?
- Do local cellular neighborhoods predict expression changes?
- Are spatial effects consistent across tissues, platforms, or patients?

Each analysis is implemented as:
- A standalone Python module
- A reproducible CLI command
- A documented method with clearly defined inputs and outputs

---

## Package structure

The framework is organized into three main layers:

### 1. Analyses (`insituprep.analyses`)

Each analysis lives in its own module and has a **dedicated README** under `docs/` describing:
- Conceptual overview
- Required inputs
- Command-line parameters
- Output formats
- Example CLI usage

Implemented analyses include:

- Expression as a function of distance
- Expression as a function of #neighbors
- DESeq2-based proximity analysis
- Dispersion of cells in PCA space
- Spatial RNA velocity analysis (stage 1 + stage 2)
- Spatial bacterial objects detection + DESeq2-based proximity analysis
- Triplet cells interactions
- Cross-platform and cross-patient comparisons

---

### 2. Command-line interface (`insituprep`)

All analyses are exposed through a unified CLI entry point:

```bash
insituprep <analysis-name> run [options]
```

To list available analyses and global options:

```bash
insituprep --help
```

---

### 3. Documentation (`docs/`)

Each analysis has its own detailed documentation file:

| Analysis | Documentation |
|--------|---------------|
| Expression as a function of distance | `docs/README_expression_as_function_of_distance.md` |
| Expression as a function of neighbors counts | `docs/README_nbrs_count_regression.md` |
| Proximity-related DESeq2 | `docs/README_deseq_proximity.md` |
| Dispersion in PCA space | `docs/README_dispersion_in_pca_space.md` |
| Spatial RNA velocity | `docs/README_rna_velocity.md` |
| Bacteria objects | `docs/README_bacteria_objects.md` |
| Expression affected by triplet cells neighborhoods | `docs/README_triplet_cell_interaction.md` |
| Platform / patient comparison | `docs/README_comparison_across_platforms_and_patients.md` |

Users are encouraged to consult the **analysis-specific README** for methodological and usage details.

---

## Installation

### Environment setup

InSituPREP must be installed inside a dedicated Python environment.  
Installing the package into the system Python is strongly discouraged.

Any standard environment manager (e.g. `venv`, `conda`, `mamba`) can be used.

The framework requires **Python ≥ 3.11**.

Example using `conda`:

```bash
conda create -n insituprep python=3.11
conda activate insituprep
```

### Installation

Clone the repository and install in editable mode:

```bash
git clone <https://github.com/Tal-Goldberg/InSituPREP.git>
cd InSituPREP_framework
pip install -e .
```
---

## Input data assumptions

Across analyses, InSituPREP assumes a **cell-level metadata table** with:

- Unique cell IDs (column named "Var1")
- Tissue identifier (column named "tissue")
- Cell type annotations (columns named "cell_type")
- Gene expression counts (raw counts, one column per gene)

Depending on the analysis, spatial information can be provided as:

- Cell centroid coordinates: `X_space`, `Y_space`, and optionally `Z_space` (in µm).
- Precomputed distance matrices (separate CSV files; analysis-specific naming and format)

Example (wide format; metadata + gene columns):

| Var1       | cell_type     | tissue | X_space        | Y_space       | Z_space      | ACTA2 | ACTG2 | ACTR3B | ADGRL4 |
|------------|---------------|--------|----------------|---------------|--------------|-------|-------|--------|--------|
| 100.26.007 | B             | 100    | 182.5507909    | 48.5805818    | 4.4420121    | 0     | 10    | 0      | 0      |
| 982.14.016 | Endothelial   | 982    | 149.6374364    | 70.9100727    | 7.3111273    | 4     | 0     | 0      | 3      |
| 313.76.019 | Epithelial    | 313    | 140.4967091    | 76.4971091    | 6.6119515    | 0     | 0     | 1      | 0      |
| 100.1.056  | T_CD3         | 100    | 121.1400273    | 36.8719009    | 4.2549212    | 1     | 0     | 0      | 0      |
| 330.5.001  | T_CD8         | 330    | 121.2483273    | 31.0298673    | 6.3241939    | 2     | 0     | 0      | 5      |

Each analysis README specifies its exact requirements.

---

## Getting started

If you are new to the framework, a recommended entry point is:

1. Read `docs/README_expression_as_function_of_distance.md`
2. Run the corresponding CLI on a single tissue
3. Explore additional analyses as complementary layers

---

## Citation

If you use InSituPREP in your work, please cite the corresponding manuscript or repository (details will be added):

*InSituPREP: InSituPREP enables 3D single-cell mapping of interaction-associated gene programs in the breast cancer tumor microenvironment*,  
Repository URL: 

A link to the preprint or published manuscript will be added here when available.

---

## License

InSituPREP is released under the MIT License.
See the LICENSE file for details.
