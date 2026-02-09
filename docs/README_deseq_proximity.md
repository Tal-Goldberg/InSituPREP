# DESeq Proximity Analysis (`insituprep deseq-proximity`)

This document describes the **DESeq proximity analysis** implemented in the `insituprep` package.
The goal of this analysis is to identify genes whose expression in a given cell type
(**primary cell type**) depends on the spatial proximity of another cell type
(**neighbor cell type**) within a tissue.

The analysis is designed for spatial transcriptomics data and combines:
- PyDESeq2-based differential expression testing
- Spatial proximity labeling
- Permutation-based statistical validation
- filtering genes by absolute log2 fold-change thresholds
- multiple-testing correction (FDR / BH adjustment).
- Optional marker-gene filtering

The pipeline is implemented as an insituprep CLI command.

To view the full list of command-line options and defaults, run:

```bash
insituprep deseq-proximity --help
```
---

## Conceptual overview: what this analysis does and why

For a selected tissue, the analysis proceeds as follows:

1. **Select a tissue**
   - Only cells belonging to the requested tissue are analyzed.

2. **Define primary and neighbor cell types**
   - Primary cells are the cells whose gene expression will be tested.
   - Neighbor cells define the spatial context.

3. **Assign spatial proximity labels**
   - Each primary cell is classified as:
     - **proximal**: at least one neighbor cell is found within a specified distance threshold (µm).
     - **non-proximal**: no neighbor cell is found within that distance (µm).
   - This converts spatial information into a binary condition suitable for differential expression testing.

4. **Differential expression analysis**
   - Gene expression in proximal vs. non-proximal primary cells is compared using PyDESeq2
     (a Python implementation of DESeq2).
   - This tests whether proximity to the neighbor cell type is associated with changes in gene expression.

5. **Permutation-based significance assessment**
   - proximity labels are randomly shuffled multiple times.
   - An empirical (permutation-based) p-value is computed for each gene,
     reflecting how often the observed signal could arise by chance.
   Conceptually, proximity labels are shuffled across primary cells multiple times and the differential expression analysis (PyDESeq2) is re-run. The permutation-based p-value reflects how often a randomized spatial configuration produces an effect as strong as, or stronger than, the observed one.

6. **Filtering and refinement**
   - Genes are filtered using an automatically selected, data-driven log2 fold-change threshold.
   - multiple-testing correction (FDR / BH adjustment) is applied
   - Genes are filtered using nominal p-value thresholds and false discovery rate (FDR) thresholds.

   - Optionally, genes that are known markers of the neighbor cell type are removed
     to reduce potential segmentation artifacts.

7. **Results are written to disk**
   - A summary table describing all analyzed cell-type pairs
   - A final gene list for each primary–neighbor pair
   
**Important note on log2 fold-change thresholds**

In this analysis, the log2 fold-change (log2FC) threshold is **not predefined by the user**.
Instead, it is **automatically selected during the analysis**, separately for each primary–neighbor cell-type pair.
The threshold is chosen from a predefined set of candidate values, based on the observed data,
in order to balance effect-size stringency while avoiding overly restrictive filtering.

---

## Filtering logic

After differential expression and permutation testing, genes are filtered in several steps to retain only robust and biologically meaningful signals.

For a given absolute log2 fold-change (log2FC) threshold, the following filtering steps are applied:

1. **Effect size filtering (optional)**  
   Genes are retained only if their absolute log2 fold change exceeds the selected threshold.  
   This step ensures that retained genes show a minimal magnitude of expression change between proximal and non-proximal primary cells.

2. **Multiple-testing correction**  
   False discovery rate (FDR) correction using the Benjamini–Hochberg procedure is applied separately to:
   - the standard (analytic) p-values obtained from differential expression testing
   - the permutation-based p-values

3. **Statistical significance filtering**  
   A gene is retained only if it satisfies *all* of the following criteria:
   - nominal analytic p-value ≤ 0.05  
   - nominal permutation-based p-value ≤ 0.05  
   - FDR-adjusted analytic p-value ≤ 0.1  
   - FDR-adjusted permutation-based p-value ≤ 0.1  

   This combined criterion ensures that the gene is significant both analytically and empirically.

4. **Marker-gene removal (optional)**  
   If marker genes are provided for the analyzed tissue, genes that are known markers of the **neighbor cell type** are removed.  
   This step is intended to reduce false positives arising from segmentation artifacts or signal contamination from neighboring cells.

---

## Log2 fold-change threshold selection

Rather than using a fixed log2FC cutoff, the analysis automatically selects an appropriate threshold for each primary–neighbor pair.

The selection procedure is as follows:

1. **Define a baseline reference**  
   First, the number of genes that are nominally significant based on both analytic and permutation p-values (without applying fold-change, FDR, or marker filtering) is recorded.  
   This serves as a reference representing the initial signal strength before stringent filtering.

2. **Evaluate candidate thresholds**  
   For each candidate log2FC threshold in the predefined set  
   `(0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)`,  
   the full filtering procedure is applied and the number of retained genes is recorded.

3. **Select the optimal threshold**  
   The log2FC threshold is chosen such that the number of final genes is closest to the baseline reference.  
   This approach balances effect-size stringency with sensitivity, avoiding overly permissive or overly restrictive thresholds.

The selected value is reported as `chosen_lfc` in the summary file and as `chosen_lfc_threshold` in the final gene list.

---

## Input files

### 1) Summary table file path (`--summary-table-path`) (required)

CSV with per-cell metadata and gene expression counts from all tissues.

**Requirements:**
- Required columns:
  - `Var1` (cell IDs) 
  - `tissue` (string tissue identifier. Even if the CSV stores tissue as numeric, the tool loads and treats it as string))
  - `cell_type` (string)
- Should include gene expression columns (one column per gene, numeric raw counts)
- Optional columns (used if distance matrix is not provided):
  - `X_space`, `Y_space` (cell centroid coordinates in µm)
  - `Z_space` (optional, for 3D analysis)
  
Example:

| Var1       | cell_type        | tissue | X_space        | Y_space        | Z_space        | ACTA2 | ACTG2 | ACTR3B | ADGRL4 |
|------------|------------------|--------|----------------|----------------|----------------|-------|-------|--------|--------|
| 100.26.007 | B                | 100    | 182.5507909    | 48.5805818     | 4.4420121      | 0     | 0     | 0      | 0      |
| 982.14.016 | Endothelial      | 982    | 149.6374364    | 70.9100727     | 7.3111273      | 4     | 0     | 0      | 0      |
| 313.76.019 | Epithelial       | 313    | 140.4967091    | 76.4971091     | 6.6119515      | 0     | 0     | 1      | 0      |
| 100.1.056  | T_CD3            | 100    | 121.1400273    | 36.8719009     | 4.2549212      | 1     | 0     | 0      | 0      |
| 330.5.001  | T_CD8            | 330    | 121.2483273    | 31.0298673     | 6.3241939      | 2     | 0     | 0      | 0      |

Notes:
- Gene expression values represent raw transcript counts per cell.
- Only cells belonging to the selected tissue are used during the analysis.
- Only genes listed in `--genes-names-path` are considered.
- Coordinates must be numeric and in micrometers (µm).

---

### 2) Distance matrix path (`--distance-matrix-path`) (optional)

CSV file containing pairwise spatial distances between cells:

distance_matrix_<TISSUE>.csv

Example (`distance_matrix_313.csv`):

|             | 313.1.001 | 313.1.021 | 313.1.022 | 313.1.023 |
|-------------|-----------|-----------|-----------|-----------|
| 313.1.001   | inf       | 425.33    | 399.38    | 537.90    |
| 313.1.021   | 425.33    | inf       | 322.86    | 658.58    |
| 313.1.022   | 399.38    | 322.86    | inf       | 302.78    |
| 313.1.023   | 537.90    | 658.58    | 302.78    | inf       |

Rows and columns correspond to cell IDs within the same tissue.
Each entry represents the pairwise spatial distance between two cells (in micrometers).
Diagonal entries are set to `inf`.

Usage:
For each primary cell, distances to all neighbor cells are queried from this matrix.
A primary cell is labeled *proximal* if at least one neighbor cell has distance ≤ the specified threshold.

If `--distance-matrix-path` is **not provided**, proximity is computed from `X_space`, `Y_space` (and optionally `Z_space` if present) in the summary table file.

---

### 3) Gene list path (`--genes-names-path`) (required)

A plain text file with **one gene name per line**.

Notes:
- Only genes present both in this list and in the summary table file are analyzed.
- Genes missing from the summary table file are ignored with a warning.
- If no overlap exists, the analysis stops with an error.


Example (`genes_names.txt`):

ACTA2
AGR2
EPCAM
VIM
CD3D
CD3E

---

### 4) Marker genes mapping (optional)  
`--marker-genes-by-tissue-json`

A JSON file mapping tissue identifiers to marker-gene CSV file path:

```json
{
  "100": "markers_100.csv",
  "313": "markers_313.csv"
}
```

Example (markers_313.csv):
| CellType     | Marker   |
|--------------|----------|
| Endothelial  | PECAM1   |
| Endothelial  | VWF      |
| Fibroblast   | COL1A1   |
| Fibroblast   | COL1A2   |
| T_CD3        | CD3D     |


Behavior:
- If this argument is **not provided**, marker filtering is skipped.
- If the tissue key is missing from the JSON mapping, or the mapped value is empty/`"None"`, marker filtering is skipped.
- If a path exists but the file is missing, a warning is printed and marker filtering is skipped.

Marker genes CSV requirements:
- Must contain columns: `CellType`, `Marker`

Implementation details:
- Marker filtering removes genes that appear as markers for the **neighbor cell type** in order to reduce potential segmentation artifacts.


### Other parameters
- `--tissue TEXT` (required)  
  Tissue identifier to analyze.  
  This parameter is treated strictly as a **string** and must match exactly one of the values (strings) in the `tissue` column of the summary table.  
  This analysis is designed to run on **one tissue per execution**.

- `--out PATH` (required)
  Output directory
  
- `--primary TEXT`  (optional)
  JSON list of primary cell types (string). Example:
  ```bash
  --primary '["Endothelial", "T_CD3"]'
  ```
  If omitted: uses **all** cell types found in the tissue.

- `--neighbor TEXT`  (optional)
  JSON list of neighbor cell types (string). Example:
  ```bash
  --neighbor '["Epithelial", "Endothelial"]'
  ```
  If omitted: uses **all** cell types found in the tissue.

- `--dist-threshold FLOAT` (optional) (default: `1`)
  Spatial distance threshold (in µm) used to classify primary cells as proximal or non-proximal relative to the specified neighbor cell type.

- `--n-perm INTEGER` (optional) (default: `1000`)
  Number of permutations
  
- `--perm-print-every INTEGER` (optional) (default: `50`)
  Print permutation progress every N permutations. 
  
- `--marker-genes-by-tissue-json PATH` (optional)
  Optional JSON mapping tissue->marker genes csv path.    
  
- `--quiet-permutations / --no-quiet-permutations` (optional) (default: quiet)
  silence PyDESeq2 output during permutations. 

---


## Examples

### run all cell types combinations

```bash
insituprep deseq-proximity run \
  --summary-table-path data/summary_table.csv \
  --distance-matrix-path data/distance_matrix_100.csv \
  --genes-names-path data/genes_names.txt \
  --tissue "100" \
  --dist-threshold 3.3 \
  --out results/deseq_out
```

### Restrict to specific primary/neighbor lists - All combinations between the specified primary and neighbor cell types are evaluated

```bash
insituprep deseq-proximity run \
  --summary-table-path data/summary_table.csv \
  --distance-matrix-path data/distance_matrix_100.csv \
  --genes-names-path data/genes_names.txt \
  --tissue "100" \
  --dist-threshold 3.3 \
  --primary '["Endothelial", "T_CD3"]' \
  --neighbor '["Epithelial", "B"]' \
  --out results/deseq_out
```

### With marker genes mapping + without distance natrix input

```bash
insituprep deseq-proximity run \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --tissue "100" \
  --dist-threshold 10 \
  --primary '["Endothelial"]' \
  --neighbor '["Epithelial"]' \
  --out results/deseq_out \
  --marker-genes-by-tissue-json data/marker_genes_by_tissue.json
  
```

---

## Output files

All output files are written to the directory provided in `--out`.

This analysis produces:
1. A **summary table** (one row per analyzed primary–neighbor pair)
2. A **final gene table** per pair (the filtered gene list that passed all criteria)

---

### 1) Summary file

**Filename**
- `deseq_summary_tissue_<TISSUE>.csv`

**Meaning**
A compact overview of the run results for each analyzed **primary–neighbor** pair (one row per pair).

**Columns (exact)**

- `tissue`  
  Tissue ID used for the run (the value passed to `--tissue`).

- `primary`  
  The primary cell type (the cells whose expression is tested).

- `neighbor`  
  The neighbor cell type (cells used to define whether a primary cell is “proximal” or “nonproximal”).

- `baseline_n`  
  A *pre-filter* count of genes that are nominally significant by both:
  - the standard DESeq2-style p-value (`pvalue` ≤ `pv_thresh`), and
  - the permutation-based p-value (`pvalue_perm` ≤ `pv_thresh`)
  
  This value is computed **before** applying:
  - the log2 fold-change cutoff filtering
  - marker-gene removal (if enabled).
  - multiple-testing correction (FDR / BH adjustment),

  
  In other words, `baseline_n` is an initial “how many genes look interesting before stricter filtering” count.

- `chosen_lfc`  
  The automatically selected absolute log2 fold-change threshold used for the final gene list.  
  The threshold is chosen from:
  `(0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45)`.

- `n_final_genes`  
  Number of genes in the final filtered output for this pair (i.e., number of rows in the final gene table).

- `out_csv`  
  Full path to the final gene table CSV written for this pair.

- `n_perm_done`  
  Number of permutations that completed successfully.  
  This can be smaller than `--n-perm` if some permutation runs failed (e.g., numerical issues).

- `marker_genes_used`  
  `TRUE` if marker genes were loaded for this tissue and marker filtering was applied; otherwise `FALSE`.

---

### 2) Final genes file

**Filename pattern**
- `final_genes_<TISSUE>_<PRIMARY>_vs_<NEIGHBOR>_LFC_<THR>.csv`

**Meaning**
This file contains the **final filtered list of genes** for a specific primary–neighbor pair.
A gene appears in this table only if it passes all filtering steps (statistical + permutation + effect-size + log2 fold change cutoff + optional marker removal + FDR).

**Columns (exact)**

- `gene`  
  Gene name.

- `baseMean` + `log2FoldChange` + `lfcSE` + `stat`  
  as reported by PyDESeq2.

- `pvalue`  
  PyDESeq2’s standard (analytic) p-value as reported by PyDESeq2.

- `padj`  
  PyDESeq2’s adjusted p-value as reported by PyDESeq2.  
  Note: the pipeline also computes its own BH-adjusted columns (see below) as part of genes filtering.

- `pvalue_perm`  
  Permutation-based p-value that quantifies how likely the observed signal is under a null model where proximity labels are random.  
  Conceptually:
  - Proximity labels are shuffled across primary cells many times,
  - PyDESeq2 is re-run,
  - The permutation p-value reflects how often a randomized shuffle produces a signal as strong as (or stronger than) the observed one.

- `p_value_adj_new`  
  Benjamini–Hochberg (BH/FDR) adjusted p-values computed from the standard `pvalue` column after log2 fold change cutoff filtering.

- `p_value_adj_new_perm`  
  Benjamini–Hochberg (BH/FDR) adjusted p-values computed from the permutation p-values (`pvalue_perm`) after log2 fold change cutoff filtering.

- `chosen_lfc_threshold`  
  The absolute log2FC threshold selected automatically for this specific pair.

- `primary_cell_type`  
  Primary cell type label for this pair.

- `neighbor_cell_type`  
  Neighbor cell type label for this pair.

- `tissue`  
  Tissue ID used for this pair.
