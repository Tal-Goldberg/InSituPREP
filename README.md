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
- Does spatial proximity contributes to within cell type variability in expression? 
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
- Biological motivation
- Statistical methodology
- Required inputs
- Output formats
- Example CLI usage

Implemented analyses include:

- Expression as a function of distance
- Expression as a function of #neighbors
- DESeq2-based proximity analysis
- Dispersion in PCA space
- Spatial bacterial objects + DESeq2-based proximity analysis
- RNA velocity (stage 1 + stage 2)
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
| Bacteria objects | `docs/README_bacteria_objects.md` |
| Expression vs distance | `docs/README_expression_as_function_of_distance.md` |
| DESeq2 proximity | `docs/README_deseq_proximity.md` |
| dispersion in PCA  | `docs/README_dispersion_in_pca_space.md` |
| Neighbor count regression | `docs/README_nbrs_count_regression.md` |
| RNA velocity | `docs/README_rna_velocity.md` |
| Platform / patient comparison | `docs/README_comparison_across_platforms_and_patients.md` |

Users are encouraged to consult the **analysis-specific README** for methodological and usage details.

---

## Installation

### Environment setup

InSituPREP must be installed inside a dedicated Python environment.  
Installing the package into the system Python is strongly discouraged.

The framework requires **Python ≥ 3.11**.

### Installation

Clone the repository and install in editable mode:

```bash
git clone <repository-url>
cd InSituPREP_framework
pip install -e .
```
---

## Input data assumptions

Across analyses, InSituPREP assumes a **cell-level table** with:

- Unique cell IDs (column named "Var1")
- Tissue identifier (column named "tissue")
- Cell type annotations (columns named "cell_type")
- Gene expression counts (raw counts, one column per gene)

Depending on the analysis, spatial information can be provided as:

- Cell centroid coordinates (`X_space`, `Y_space`, optional `Z_space`)
- Precomputed distance matrices

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

*InSituPREP: A multi-layered framework for proximity-aware analysis of spatial transcriptomics data*,  
Repository URL: <repository-url>

A link to the preprint or published manuscript will be added here when available.

---

## License

InSituPREP is released under the MIT License.
See the LICENSE file for details.
