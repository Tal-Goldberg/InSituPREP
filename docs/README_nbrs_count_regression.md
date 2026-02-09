# Neighbor Count Regression Analysis (`insituprep nbrs-count-regression`)

This document describes the **Neighbor Count Regression Analysis** implemented in the `insituprep` package.
The goal of this analysis is to identify genes whose expression in a given cell type
(**primary cell type**) correlates with the **number of neighboring cells** of another cell type
(**neighbor cell type**) within a spatial distance threshold.

Unlike proximity-based differential expression approaches, this analysis tests for a **continuous relationship**
between gene expression and neighbor count, rather than a binary proximal/non-proximal comparison.

The method is designed for spatial transcriptomics data and integrates:
- Three regression methods (linear, sampling-based, weighted least squares)
- Spatial neighbor counting using cKDTree or pre-computed distance matrices
- Gene filtering by expression threshold and variance quantile
- Multiple-testing correction (FDR / BH adjustment)
- Visualization (4-panel gene plots, cross-tissue heatmaps)
- Optional marker-gene filtering

To see all command-line options, run:

```bash
insituprep nbrs-count-regression --help
```
The pipeline is implemented as an insituprep CLI command.

---

## Conceptual overview: what this analysis does and why

For each tissue and each ordered pair of cell types (primary → neighbor), the analysis proceeds as follows:

1. **Select tissue(s)**
   - All cells belonging to each selected tissue are analyzed independently.
   - If no tissue is specified, all tissues in the data are analyzed.

2. **Define primary and neighbor cell types**
   - **Primary cell type**: cells whose gene expression is analyzed.
   - **Neighbor cell type**: cells used only to compute neighbor counts.
   - All ordered pairs (primary ≠ neighbor) are evaluated unless restricted by the user.

3. **Count neighbors within distance threshold**
   - For each primary cell, count how many neighbor cells lie within the distance threshold (µm).
   - Distance is computed from coordinates (`X_space`, `Y_space`, optionally `Z_space`) using cKDTree.
   - Alternatively, pre-computed distance matrices can be provided.

4. **Filter genes**
   - **Expression filter**: Keep genes with max expression > threshold (default: 10).
   - **Variance filter**: Keep top X% most variable genes (default: top 20%).
   - Optionally, marker genes of the neighbor cell type can be removed.

5. **Regression analysis**
   - For each gene, regress expression ~ neighbor_count using one of three methods:
   
   **(a) Method 1: Linear regression**
   - Standard OLS regression of expression on neighbor count.
   - Fast and suitable when neighbor count distribution is balanced.
   
   **(b) Method 2: Sampling-based regression**
   - Addresses imbalanced neighbor count distributions.
   - **2a**: Average slopes from bootstrap samples, compute optimal intercept.
   - **2b**: Fit regression line to all sampled lines ("line-on-line").
   
   **(c) Method 3: Weighted least squares**
   - Weights cells inversely proportional to their neighbor count frequency.
   - **3.1**: Unnormalized weights.
   - **3.2**: Normalized weights (sum to 1).

6. **Multiple-testing correction**
   - Benjamini–Hochberg (FDR) correction is applied across all tested genes.

7. **Results are written to disk**
   - CSV files with all results and significant-only results
   - Heatmaps showing -log10(q-value) across tissues
   - 4-panel gene plots for significant genes

---

## Input files


### 1) Summary table file path (`--summary-table-path`) (required)

CSV with per-cell metadata and gene expression counts from all tissues.

**Requirements:**
- Required columns:
  - `Var1` (cell IDs) 
  - `tissue` (string tissue identifier. Even if the CSV stores tissue as numeric, the tool loads and treats it as string)
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
- `Z_space` is optional — if absent or all NaN, 2D distances are used automatically.

---

### 2) Gene list (`--genes-names-path`) — required

Plain text file with **one gene name per line**.

Notes:
- Only genes present in both this list and the summary table file are used.
- Missing genes are ignored with a warning.
- If no overlap exists, the analysis stops.


Example (`genes_names.txt`):

ACTA2
AGR2
EPCAM
VIM
CD3D
CD3E

---

### 3) Distance matrices directory (`--distance-dir`) — optional

Directory containing CSV files with pre-computed pairwise distances:

```
distance_matrix_<TISSUE>.csv
```

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
Neighbor count is the number of neighbor cells with distance ≤ the specified threshold.

If `--distance-dir` is **not provided**, neighbor counts are computed from `X_space`, `Y_space` (and optionally `Z_space`) in the summary table file using cKDTree.

---

### 4) Marker genes mapping (optional)
`--marker-genes-by-tissue-json`

JSON mapping tissue IDs to marker-gene CSVs or gene lists:

```json
{
  "100": "markers_100.csv",
  "313": "markers_313.csv"
}
```

Example (markers_100.csv):
| CellType     | Marker   |
|--------------|----------|
| Endothelial  | PECAM1   |
| Endothelial  | VWF      |
| Fibroblast   | COL1A1   |
| T_CD3        | CD3D     |


Behavior:
- If this argument is **not provided**, marker filtering is skipped.
- If the tissue key is missing from the JSON mapping, or the mapped value is empty/`"None"`, marker filtering is skipped.
- If a path exists but the file is missing, a warning is printed and marker filtering is skipped.

Implementation details:
- Marker genes are removed before regression to reduce potential segmentation artifacts.

---

## Command-line parameters

- `--summary-table-path`, `-s` PATH (required)  
  CSV with cell data (cell_id, tissue, cell_type, coordinates, genes).

- `--genes-names-path`, `-g` PATH (required)  
  Text file with gene names to analyze (one per line).

- `--out`, `-o` PATH (required)  
  Output directory for results and plots.

- `--tissue`, `-t` TEXT (optional)  
  JSON list of tissue IDs as strings. Example:
  ```bash
  --tissue '["100","313"]'
  ```
  If omitted: uses **all** tissues found in the data.

- `--distance-dir` PATH (optional)  
  Directory containing pre-computed distance matrices.
  If not provided, neighbor counts are computed from coordinates using cKDTree.

- `--dist-threshold` FLOAT (optional) (default: `15.0`)  
  Spatial distance threshold for counting neighbors (in the same units as the coordinate columns).

- `--method`, `-m` TEXT (optional) (default: `1`)  
  Regression method: `1`, `2a`, `2b`, `3.1`, or `3.2`.

- `--min-expression` INTEGER (optional) (default: `10`)  
  Keep genes with max expression > threshold.

- `--variance-quantile` FLOAT (optional) (default: `0.20`)  
  Keep top X% most variable genes (0.20 = top 20%).

- `--min-neighbor-bins` INTEGER (optional) (default: `4`)  
  Minimum neighbor count bins required (Method 2/3).

- `--min-cells-per-bin` INTEGER (optional) (default: `10`)  
  Minimum cells per bin (Method 2/3).

- `--n-iterations` INTEGER (optional) (default: `1000`)  
  Number of sampling iterations (Method 2).

- `--sample-size` INTEGER (optional) (default: `10`)  
  Cells sampled per bin per iteration (Method 2).

- `--fdr-alpha` FLOAT (optional) (default: `0.05`)  
  FDR significance threshold.

- `--primary-cell-types` TEXT (optional)  
  JSON list of primary cell types. Example:
  ```bash
  --primary-cell-types '["Endothelial","Smooth muscle"]'
  ```
  If omitted: uses **all** cell types found in the tissue.

- `--neighbor-cell-types` TEXT (optional)  
  JSON list of neighbor cell types.
  If omitted: uses **all** cell types found in the tissue.

- `--marker-genes-by-tissue-json` PATH (optional)  
  JSON mapping tissue ID to marker genes CSV path.

- `--rng-seed` INTEGER (optional)  
  Random seed for reproducibility (Method 2 sampling).

- `--workers`, `-w` INTEGER (optional) (default: `1`)  
  Number of parallel workers.
  Positive = exact count, -1 = all CPUs, -2 = all CPUs minus 1.

---

## Output files

All outputs are written under the directory specified by `--out`.

This analysis produces:
1. **Results CSV files** (all results + significant-only)
2. **Heatmaps** showing -log10(q-value) across tissues
3. **Gene plots** (4-panel visualization for each significant gene)

---

### 1) Results CSV files

**Filename pattern**
- `method{M}_{primary}_vs_{neighbor}_dist{D}_all.csv` — All results
- `method{M}_{primary}_vs_{neighbor}_dist{D}_significant.csv` — FDR-significant only

**Meaning**
Each row corresponds to one gene analyzed for a specific primary–neighbor pair.

**Columns (exact)**

- `tissue`  
  Tissue ID used for the run.

- `target_cell`  
  The primary cell type (cells whose expression is analyzed).

- `neighbor_cell`  
  The neighbor cell type (cells used to compute neighbor counts).

- `gene`  
  Gene name.

- `method`  
  Regression method used (`1`, `2a`, `2b`, `3.1`, or `3.2`).

- `p_value`  
  Raw p-value from the regression analysis.

- `q_value`  
  Benjamini–Hochberg FDR-adjusted p-value.

- `r_value`  
  Correlation coefficient (Pearson r for Method 1).

- `r_squared`  
  Coefficient of determination (R²).

- `slope`  
  Regression slope.
  Positive slope: expression increases with more neighbors.
  Negative slope: expression decreases with more neighbors.

- `intercept`  
  Regression intercept.

- `n_cells`  
  Number of primary cells included in the analysis.

---

### 2) Heatmaps

**Filename pattern**
- `plots/heatmap_method{M}_{primary}_vs_{neighbor}_dist{D}.png` — Per cell pair
- `plots/heatmap_method{M}_{primary}_vs_{neighbor}_dist{D}_ALL_TISSUES.png` — Across tissues

**Meaning**
Heatmaps show -log10(q-value) for significant genes.
Rows are genes, columns are tissues.
Higher values (brighter colors) indicate stronger significance.

---

### 3) Gene plots

**Location:** `plots/genes/`

**Filename pattern**
- `{gene}_T{tissue}_method{M}_{primary}_vs_{neighbor}_dist{D}.png`

**Meaning**
4-panel visualization for each significant gene:

- **Top-left**: Expression vs. neighbor count scatter + regression line
- **Top-right**: Neighbor count histogram
- **Bottom-left**: Mean expression ± std vs. neighbor count
- **Bottom-right**: Fit on raw data (Method 2/3 only)

---

## Example usage

### Run all tissues with Method 1 (linear regression)
```bash
insituprep nbrs-count-regression run \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --out results/nbrs_method1 \
  --method 1 \
  --dist-threshold 15
```

### Run specific tissues with Method 2a (sampling-based)
```bash
insituprep nbrs-count-regression run \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --out results/nbrs_method2a \
  --tissue '["100", "982"]' \
  --method 2a \
  --dist-threshold 15 \
  --rng-seed 42
```

### Restrict to specific cell types
```bash
insituprep nbrs-count-regression run \
  --summary-table-path data/summary_10tissues_micron.csv \
  --genes-names-path data/genes_names.txt \
  --out results/nbrs_subset \
  --primary-cell-types '["Endothelial","T_CD3"]' \
  --neighbor-cell-types '["Epithelial","B"]' \
  --method 1
```

### With parallel processing
```bash
insituprep nbrs-count-regression run \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --out results/nbrs_parallel \
  --method 1 \
  --workers -1
```

### With marker gene filtering
```bash
insituprep nbrs-count-regression \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --out results/nbrs_markers \
  --method 1 \
  --marker-genes-by-tissue-json data/marker_genes.json
```

### With marker gene filtering and distance matrix 
```bash
insituprep nbrs-count-regression \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --out results/nbrs_markers \
  --method 1 \
  --marker-genes-by-tissue-json data/marker_genes.json
  --distance-dir data/distance_matrix \
```

---

## Interpretation notes

- A significant **positive slope** indicates that gene expression **increases** with more neighbors.
- A significant **negative slope** indicates that gene expression **decreases** with more neighbors.
- The heatmap provides a quick overview of which genes are significant across multiple tissues.
- Gene plots help visually validate whether the regression fit is meaningful.

---

## Recommended usage

- Use Method 1 for initial exploration (fastest).
- Use Method 2 or 3 when neighbor count distributions are highly imbalanced.
- Use `--rng-seed` for reproducible results with Method 2.
- Use `--workers -1` to leverage all CPU cores for large datasets.
