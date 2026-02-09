# Triplet Cell Interaction (`insituprep triplet-cell-interaction`)

This document describes the **Triplet Cell Interaction** analysis implemented in the `insituprep` package.

The goal of this analysis is to detect **gene-expression changes in a focal(primary) cell type** that are associated with the
**presence of specific multi–cell-type neighborhoods ("triplets")** within a user-defined spatial radius.
Concretely, the method:
- uses cell centroid coordinates (`X_space`, `Y_space`, optional `Z_space`) to define local neighborhoods,
- enumerates **all cell-type combinations automatically** (no need to provide lists of cell types),
- builds DESeq2 comparisons between *different neighborhood compositions* for the **same focal cell type**,
- runs **PyDESeq2** to identify genes that differ between the compared neighborhood contexts.

The pipeline is implemented as an `insituprep` CLI command.

To see all command-line options, run:

```bash
insituprep triplet-cell-interaction --help
```

---

## Conceptual overview: what this analysis does and why

### What is a "triplet" in this implementation?

For each **focal cell** (a single cell in the tissue), we define its neighborhood as:

> the set of **unique cell types** observed among all cells whose centroids fall within a radius `attach` (µm)
> around the focal cell.

This neighborhood type-set is stored as a **comma-separated string** called `friends`, e.g.:
- `"Endothelial,Smooth muscle"`  (2 unique types in the radius)
- `"Endothelial,Epithelial,T_CD3"`  (3 unique types in the radius)

The analysis keeps only focal cells whose neighborhood contains **≤ 3 unique cell types**.
(Internally: neighborhood strings with at most 2 commas.)

### The main question being tested

For each focal cell type **i**, and for each neighborhood pattern **k** containing **at least two cell types**, the analysis compares:

- **Group G (baseline context):** focal cells of type **i** whose neighborhood type-set is exactly `k`
- **Group A (augmented context):** focal cells of the same type **i** whose neighborhood type-set is exactly `k` **plus one additional cell type j**

This produces a DESeq2 test per tuple:
(**focal type i**, **pair-types k**, **added type j**)

---

## Methodology (step-by-step)

1. **Load the summary table**
   - Read a per-cell CSV containing metadata, coordinates, and gene expression counts.

2. **Optional: restrict to one tissue**
   - If `--tissue` is provided, the table is filtered to:
     ```text
     df[tissue_col] == tissue
     ```
   - **Important:** `tissue` is treated as a **string** (even if the CSV stores it as numbers).

   If `--tissue` is not provided, the analysis runs once on the full table (all tissues together),
   which may mix cells from different tissues.

3. **Choose coordinate dimensionality (2D vs 3D) automatically**
   - If columns `X_space`, `Y_space`, `Z_space` exist → use **3D** coordinates.
   - Otherwise, if `X_space`, `Y_space` exist → use **2D** coordinates.
   - If neither condition is met → the analysis fails with an error.

   (Implementation: `_get_coords()` chooses `[[X,Y,Z]]` if present else `[[X,Y]]`.)  

4. **Build a spatial index and collect neighborhood types**
   - Build a `scipy.spatial.cKDTree` from the coordinates.
   - For each cell `i`, query all cells within radius `attach`:
     ```python
     idxs = tree.query_ball_point(coords[i], r=attach)
     ```
   - Extract the **unique** cell types among these indices.
   - Create a sorted comma-separated string `friends = ",".join(sorted(unique_types))`.
   - Keep only cells where `friends` contains **≤ 3 unique cell types**.

   Output of this step is a table with:
   - `cell`: focal cell ID
   - `type`: focal cell type
   - `friends`: neighborhood type-set string

5. **Build DESeq2 inputs (counts + metadata) for each triplet comparison**
   - For each focal type `i_type`:
     - For each neighborhood string `k` that contains at least one comma (≥ 2 types):
       - Define **Group G** as focal cells of type `i_type` whose `friends == k`.
       - For every other cell type `j_type != i_type`:
         - Define **Group A** as focal cells of the same `i_type` whose neighborhood is exactly:
           - `k + "," + j_type` **or** `j_type + "," + k`
           - (both orders are checked)
         - If both groups contain at least one cell, build a DESeq2 dataset:
           - `counts_df`: stacked gene-count matrix for Group G and Group A
           - `meta_df`: one column `condition` with values `"G"` or `"A"`

6. **Run DESeq2 (PyDESeq2) per triplet comparison**
   - Design: `~ condition`
   - Contrast: `["condition", "G", "A"]`
     - Interpreted as **G vs A** (log2 fold-change for G relative to A).

7. **Filter DESeq2 results**
   - Keep only genes with finite p-values and:
     - `pvalue < alpha` (where `alpha` is the CLI `--alpha` parameter)

8. **Combine results across all comparisons**
   - Results across all triplet comparisons are concatenated into one output table.
   - The group key is parsed into 3 explicit columns:
     - `triplet_focal_type`
     - `triplet_pair_types`
     - `triplet_neighbor_type`

---

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

---

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

---

## Command-line parameters

- `--summary-table-path PATH` (required)  
  Per-cell CSV with metadata, coordinates, and gene counts.

- `--genes-names-path PATH` (required)  
  Text file with one gene name per line.

- `--out PATH` (required)  
  Output directory.

- `--tissue TEXT` (optional)  
  Tissue identifier (**string**).  
  If provided, analysis is restricted to this tissue only.  
  If not provided, analysis is run once on the entire table (all tissues together).

- `--attach FLOAT` (default: `10.0`)  
  Neighborhood radius in **µm** used to query neighboring cells for each focal cell.

- `--alpha FLOAT` (default: `0.05`)  
  P-value threshold used to filter the DESeq2 results (`pvalue < alpha`).

- `--n-cpus INT` (default: `8`)  
  CPU count passed to PyDESeq2 inference.

- `--tissue-col TEXT` (default: `tissue`)  
  Column name for tissue identifier.

- `--celltype-col TEXT` (default: `cell_type`)  
  Column name for cell type annotation.

- `--cell-id-col TEXT` (default: `Var1`)  
  Column name for cell IDs.

---

## Output files

All outputs are written under the directory specified by `--out`.

### 1) Main results table  
`triplets_deseq_tissue_<TISSUE>.csv`  (if `--tissue` is provided)  
or  
`triplets_deseq_all_tissues.csv`  (if `--tissue` is not provided)


Each **row** corresponds to:
- a **single gene** that passed the p-value filter, for
- a specific triplet comparison (**focal type i**, **pair-types k**, **added type j**).

#### Triplet identifier columns (added by this analysis)

- `triplet_focal_type`  
  The focal cell type **i**. DESeq2 is run using only cells of this type.

- `triplet_pair_types`  
  The **baseline neighborhood type-set** `k` (comma-separated cell types, sorted).
  This is the neighborhood composition for **Group G**.

- `triplet_neighbor_type`  
  The **additional cell type** `j` that is present in the **augmented neighborhood**.
  Group A neighborhoods contain the types in `triplet_pair_types` **plus** this type.

#### DESeq2 result columns (from PyDESeq2)

The remaining columns come directly from `PyDESeq2`'s `results_df`.
Typical columns include (names may depend on PyDESeq2 version, but these are the standard ones):

- `baseMean` — mean of normalized counts across all samples
- `log2FoldChange` — log2 fold-change for the specified contrast (**G vs A**)
- `lfcSE` — standard error of the log2 fold-change
- `stat` — Wald test statistic
- `pvalue` — Wald test p-value
- `padj` — adjusted p-value (if present in the results object)

**Filtering applied by this analysis:**
- keep only finite `pvalue`
- keep only `pvalue < alpha`

---

## Example usage

### Run on a single tissue
```bash
insituprep triplet-cell-interaction run \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --tissue "100" \
  --attach 10 \
  --alpha 0.05 \
  --n-cpus 8 \
  --out results/triplet_cell_interaction_tissue100
```

### Run on all tissues together (single pooled run)
```bash
insituprep triplet-cell-interaction run \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --attach 10 \
  --alpha 0.05 \
  --n-cpus 8 \
  --out results/triplet_cell_interaction_all_tissues
```

---

## Interpretation notes

- The analysis is **neighborhood-composition aware**: it tests whether adding a specific cell type `j`
  to a baseline neighborhood `k` is associated with differential expression in the focal type `i`.
- The DESeq2 comparison is **within the same focal cell type** (reduces confounding by cell-type identity).
- Running without `--tissue` may mix tissues and introduce batch/tissue effects into DESeq2,
  because the run is performed as one pooled dataset.

---

## Recommended usage

- Prefer running **one tissue at a time** (`--tissue`) unless you intentionally want a pooled analysis.
- Start with an `attach` radius that matches your biological scale (e.g., immediate neighborhood vs broader microenvironment).

---