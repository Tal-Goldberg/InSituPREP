# Triplet Cell Interaction (`insituprep triplet-cell-interaction`)

This document describes the **Triplet Cell Interaction** analysis
implemented in the `insituprep` package.

The goal of this analysis is to detect **gene-expression changes in a
primary (primary) cell type** that are associated with the **presence of
specific multi--cell-type neighborhoods ("triplets")** within a
user-defined spatial radius.

Concretely, the method: - uses cell centroid coordinates (`X_space`,
`Y_space`, optional `Z_space`) to define local neighborhoods, -
enumerates **all cell-type combinations automatically** (no need to
provide lists of cell types), - builds DESeq2 comparisons between
*different neighborhood compositions* for the **same primary cell
type**, - runs **PyDESeq2** to identify genes that differ between the
compared neighborhood contexts.

The pipeline is implemented as an `insituprep` CLI command.

To see all command-line options, run:

``` bash
insituprep triplet-cell-interaction --help
```

------------------------------------------------------------------------

## Conceptual overview: what this analysis does and why

### What is a "triplet" in this implementation?

For each **primary cell**, we define its neighborhood as:

> the set of **unique cell types** observed among all cells whose
> centroids fall within a radius `attach` (µm) around the primary cell.

This neighborhood type-set is stored as a **comma-separated string**
called `friends`, for example: - `"Endothelial,Smooth muscle"` -
`"Endothelial,Epithelial,T_CD3"`

Only primary cells whose neighborhood contains **≤ 3 unique cell types**
are retained (i.e. neighborhood strings with at most two commas).

### The main question being tested

For each primary cell type **i**, and for each baseline neighborhood
composition **k** containing **at least two cell types**, the analysis
compares:

-   **Group G (baseline context):** primary cells of type **i** whose
    neighborhood is exactly `k`
-   **Group A (augmented context):** primary cells of the same type **i**
    whose neighborhood is exactly `k + "," + j`, where **j** is an
    additional cell type

Each comparison corresponds to a tuple:

(**primary type i**, **baseline neighborhood k**, **added cell type j**)

**Important:** neighborhood matching is **order-sensitive**. Only exact
matches of the form `k + "," + j` are tested; the reverse order is not
considered.

------------------------------------------------------------------------

## Methodology (step-by-step)

1.  **Load the summary table**
    -   Read a per-cell CSV containing metadata, spatial coordinates,
        and gene expression counts.
2.  **Optional: restrict to one tissue**
    -   If `--tissue` is provided, the table is filtered to cells
        belonging to that tissue.
    -   The tissue identifier is always treated as a **string**.
    -   If `--tissue` is not provided and a tissue column exists, the
        analysis is run **separately for each tissue** and results are
        concatenated.
    -   If no tissue column exists, the analysis is run once on the full
        table.
3.  **Choose coordinate dimensionality (2D vs 3D) automatically**
    -   If `X_space`, `Y_space`, and `Z_space` exist → use **3D**
        coordinates.
    -   Otherwise, if `X_space` and `Y_space` exist → use **2D**
        coordinates.
    -   If neither condition is met, the analysis fails.
4.  **Collect neighborhood type-sets**
    -   For each cell *i*, compute distances to all other cells using
        `scipy.spatial.distance.cdist`.
    -   Cells with more than one neighbor (including the primary cell
        itself) within radius `attach` are retained.
    -   For each retained cell, extract the unique cell types within the
        radius, sort them, and concatenate into a comma-separated
        `friends` string.
    -   Only neighborhoods with **≤ 3 unique cell types** are kept.

    The output of this step is a table with:
    -   `cell`: primary cell ID
    -   `type`: primary cell type
    -   `friends`: neighborhood type-set string
5.  **Build DESeq2 inputs**
    -   For each primary cell type **i**:
        -   For each baseline neighborhood `k` containing at least two
            cell types:
            -   **Group G:** cells with `friends == k`
            -   For each candidate added type **j**:
                -   **Group A:** cells with `friends == k + "," + j`
            -   If both groups contain at least one cell, a DESeq2
                dataset is built.
6.  **Run DESeq2 (PyDESeq2)**
    -   Design formula: `~ condition`
    -   Contrast: `["condition", "G", "A"]`
    -   Results are returned **without internal p-value filtering**.
7.  **Multiple-testing correction and filtering**
    -   Filtering behavior depends on the presence of a tissue column:
        -   **No tissue column:** results are returned without p-value
            or FDR filtering.
        -   **With tissue column:** for each tissue separately:
            -   non-finite p-values are set to 1
            -   adjusted p-values (`padj`) are computed using
                `scipy.stats.false_discovery_control`
            -   only genes with `padj < pv_threshold` are retained.

------------------------------------------------------------------------

## Input files

### 1) Summary table file path (`--summary-table-path`) (required)

CSV with per-cell metadata, coordinates, and gene expression.

**Required columns (default names):**
- `Var1` — unique cell IDs
- `tissue` — tissue identifier (**treated as string**)
- `cell_type` — cell type annotation
- `X_space`, `Y_space` — cell centroid coordinates in µm

**Optional:**
- `Z_space` — if present, analysis runs in 3D

**Gene expression columns:**
- One column per gene (raw numeric counts per cell)

Example:

| Var1       | tissue | cell_type       | X_space | Y_space | Z_space | ACTA2 | EPCAM |
|------------|--------|-----------------|--------:|--------:|--------:|------:|------:|
| 100.1.001  | 100    | Endothelial     | 182.55  | 48.58   | 4.44    | 0     | 1     |
| 100.1.002  | 100    | Smooth muscle   | 149.64  | 70.91   | 7.31    | 4     | 0     |

Notes:
- Units for spatial coordinates should be **micrometers (µm)**.
- The neighborhood radius `attach` is interpreted in the same units.

------------------------------------------------------------------------

### 2) Gene list (`--genes-names-path`) — required

Plain text file with **one gene name per line**.

Notes (current behavior):
- The analysis **expects all listed genes to exist** in the summary table.
  If a listed gene is missing, the run stops with an error.

Example:

```text
ACTA2
EPCAM
VIM
CD3D
```

------------------------------------------------------------------------

## Command-line parameters

-   `--summary-table-path PATH` (required)\
    Per-cell CSV with metadata, coordinates, and gene counts.
    
-   `--genes-names-path PATH` (required)\
    Text file with one gene name per line.
    
-   `--out PATH` (required)\
    Output directory.
    
-   `--tissue TEXT` (optional)\
    Tissue identifier (**string**).  
    If provided, analysis is restricted to this tissue only.  
    If not provided, analysis is run once on the entire table (all tissues together).

-   `--attach FLOAT` (default: 10.0)\
    Neighborhood radius in **µm** used to query neighboring cells for each focal cell.

-   `--pv-threshold FLOAT` (default: 0.05)\
    Adjusted p-value (FDR) threshold used to filter DESeq2 results **when a tissue column exists**.  
    For each tissue separately, genes with `padj < pv-threshold` are retained.  
    If the summary table does not contain a tissue column, no p-value or FDR filtering is applied.
    
-   `--n-cpus INT` (default: 8)\
    CPU count passed to PyDESeq2 inference.
    
-   `--tissue-col TEXT` (default: `tissue`)\
    Column name for tissue identifier.
    
-   `--celltype-col TEXT` (default: `cell_type`)\
    Column name for cell type annotation.
    
-   `--cell-id-col TEXT` (default: `Var1`)
    Column name for cell IDs.
    
------------------------------------------------------------------------

## Output

A single CSV file written under the output directory:

-   `triplets_deseq_tissue_<TISSUE>.csv` if `--tissue` is provided
-   `triplets_deseq_all.csv` otherwise

Each row corresponds to one gene tested in one triplet comparison.

### Triplet annotation columns

-   `primary type` --- primary cell type (**i**)
-   `neighbor type` --- baseline neighborhood composition
-   `neighbor added type` --- added cell type (**j**)

### DESeq2 result columns

All remaining columns are taken directly from PyDESeq2 output
(e.g. `baseMean`, `log2FoldChange`, `stat`, `pvalue`, `padj`).

------------------------------------------------------------------------

## Example usage

``` bash
insituprep triplet-cell-interaction run \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --tissue "100" \
  --attach 10 \
  --pv-threshold 0.05 \
  --out results/triplet_cell_interaction_tissue100
```

------------------------------------------------------------------------

## Interpretation notes

-   All comparisons are performed **within the same primary cell type**.
-   Neighborhood composition matching is **order-sensitive**.
