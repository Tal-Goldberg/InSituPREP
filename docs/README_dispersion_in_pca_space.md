# Dispersion in PCA Space (`insituprep dispersion-pca`)

This document describes the **Dispersion in PCA space analysis** implemented in the `insituprep` package.
The goal of this analysis is to test whether **cells of a given primary cell type occupy distinct regions
in a low-dimensional expression space depending on their spatial proximity to another cell type**.

Unlike differential expression–based approaches, this analysis focuses on **multivariate dispersion
and spatial structure in expression space**, rather than gene-by-gene changes.

The method is designed for spatial transcriptomics data and integrates:
- PCA or t-SNE embedding of primary-cell expression profiles
- Spatial proximity labeling using distance matrices
- Permutation-based null models
- Max-type test statistics combining multiple dispersion metrics
- Multiple-testing correction (FDR / BH)
- Several sensitivity analyses (marker removal, shuffled proximity, t-SNE)

The pipeline is implemented as an insituprep CLI command.

To see all command-line options, run:

```bash
insituprep dispersion-pca --help
```

---

## Conceptual overview: what this analysis does and why

For each tissue and each ordered pair of cell types (primary → neighbor), the analysis proceeds as follows:

1. **Select tissue(s)**
   - All cells belonging to each selected tissue are analyzed independently.

2. **Define primary and neighbor cell types**
   - **Primary cell type**: cells whose expression profiles are embedded and analyzed.
   - **Neighbor cell type**: cells used only to define spatial proximity.
   - All ordered pairs (primary ≠ neighbor) are evaluated unless restricted by the user.

3. **Compute expression embedding**
   - For each primary cell type:
     - Gene expression vectors are standardized.
     - A 2D embedding is computed:
       - **PCA** for the main analysis
       - **t-SNE** for sensitivity analyses

4. **Assign spatial proximity labels**
   - For each primary cell:
     - **Proximal** if at least one neighbor cell lies within a distance threshold (µm).
     - **Distant** otherwise.
   - Proximity is a **binary attribute** derived from continuous spatial distances.

5. **Quantify dispersion in embedding space**
   Two complementary statistics are computed:

   **(a) Single-centroid statistic**
   - Sum of distances from proximal cells to the overall centroid of all primary cells.

   **(b) Two-centroid statistic**
   - Euclidean distance between the centroids of proximal and distant primary cells.

6. **Permutation-based significance testing**
   - Proximity labels are randomly shuffled across primary cells.
   - The two statistics are recomputed for each permutation.
   - Empirical null distributions are generated.

7. **Max-type test**
   - For each pair, a max statistic is computed:
     ```
     max(z_first_statistic, z_second_statistic)
     ```
   - This captures the strongest dispersion signal across metrics.
   - A Gaussian-approximated p-value is computed for the max statistic.

8. **Multiple-testing correction**
   - Benjamini–Hochberg (FDR) correction is applied across all tested pairs.

9. **Sensitivity analyses**
   - Shuffled proximity (sanity check)
   - Marker-gene removal to reduce potential segmentation artifacts.
   - t-SNE instead of PCA
   - t-SNE with shuffled proximity labels (all pairs)

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
| 313.76.019 | Epithelial.      | 313    | 140.4967091    | 76.4971091     | 6.6119515      | 0     | 0     | 1      | 0      |
| 100.1.056  | T_CD3            | 100    | 121.1400273    | 36.8719009     | 4.2549212      | 1     | 0     | 0      | 0      |
| 330.5.001  | T_CD8            | 330    | 121.2483273    | 31.0298673     | 6.3241939      | 2     | 0     | 0      | 0      |

Notes:
- Gene expression values represent raw transcript counts per cell.
- Only cells belonging to the selected tissue are used during the analysis.
- Only genes listed in `--genes-names-path` are considered.
- Coordinates must be numeric and in micrometers (µm).

---

### 2) Distance matrices folder (`--distance-dir`) - optional

Directory containing CSV files with pre-computed pairwise distances.

Important:
You must provide the directory path only (do not provide a specific distance matrix CSV file).
The directory must contain one CSV file per tissue, and each file name must follow exactly this naming convention:

```
distance_matrix_<TISSUE>.csv
```
Where <TISSUE> matches the tissue ID used in the summary table.

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

For each tissue <TISSUE>, the tool automatically loads:
```
<distance_dir>/distance_matrix_<TISSUE>.csv
```

Usage:
For each primary cell, distances to all neighbor cells are queried from this matrix.
A primary cell is labeled *proximal* if at least one neighbor cell has distance ≤ the specified threshold.

If `--distance-dir` is **not provided**, proximity is computed from `X_space`, `Y_space` (and optionally `Z_space` if present) in the summary table file.

---

### 3) Gene list (`--genes-names-path`) — required

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

### 4) Marker genes mapping - optional
`--marker-genes-by-tissue-json`

JSON mapping tissue IDs to marker-gene CSV file path:

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

---

## Command-line parameters

- `--out PATH` (required)  
  Output directory

- `--dist-threshold FLOAT` (default: `1.0`)  
  Distance threshold (µm) for proximity labeling between primary-type cell and neighbor-type cell

- `--n-perm INTEGER` (default: `10000`)  
  Number of permutations

- `--fdr-alpha FLOAT` (default: `0.05`)  
  FDR cutoff for significant pairs

- `--tissue-ids TEXT` (optional)  
  JSON list of tissues to analyze  
  Example:
  ```bash
  --tissue '["100","313"]'
  ```

- `--primary TEXT` (optional)  
  JSON list of primary cell types
  ```bash
  --primary '["Endothelial", "T_CD3"]'
  ```
  If omitted: uses **all** cell types found in the tissue, evaluating all ordered (primary → neighbor) pairs.
  
- `--neighbor TEXT` (optional)  
  JSON list of neighbor cell types
  ```bash
  --neighbor '["Epithelial", "Endothelial"]'
  ```
  If omitted: uses **all** cell types found in the tissue, evaluating all ordered (primary → neighbor) pairs.
    
- `--tsne-perplexity INTEGER` (default: `30`)  
  Perplexity for t-SNE sensitivity analysis

- `--tsne-random-state INTEGER` (default: `None`)  
  Random seed for t-SNE

- `--rng-seed INTEGER` (default: `None`)  
  Random seed for reproducibility
---

## Output files

All outputs are written under the directory specified by `--out`.

---
### 1) All-pairs PCA results  
`pca_pairs_all_results.csv`  
(and its filtered version: `pca_pairs_significant.csv`)

Each row corresponds to one ordered pair:

(**tissue**, **primary cell type → neighbor cell type**)

Expression profiles of the *primary cell type* are embedded in PCA space, and spatial proximity is defined relative to the *neighbor cell type*.

---

#### Identifier columns
- `tissue`  
  Tissue ID, as defined in the summary table file.

- `cell_type_1` *(primary cell type)*  
  The cell type whose gene expression profiles are embedded in PCA space and analyzed.

- `cell_type_2` *(neighbor cell type)*  
  The cell type used only to define spatial proximity of the primary cells
  (proximal vs distant).

---

#### Basic counts
- `num_cells_ct1_used`  
  Number of primary cell type cells included in the analysis for this pair,
  after restricting to the tissue, selecting the gene list, and removing cells
  with missing expression values.

- `num_proximal`  
  Number of primary cell type cells that have at least one neighbor cell of the
  specified neighbor cell type within the distance threshold.

---

#### Raw distance statistics in PCA space

- `original_dist_single_centroid`  
  Sum of Euclidean distances from *proximal* primary cells to the overall centroid
  of all primary cell type cells in PCA space.

- `original_dist_all_to_centroid`  
  Sum of Euclidean distances from *all* primary cell type cells to the overall
  centroid in PCA space.

- `original_dist_two_centroids`  
  Euclidean distance between the centroid of proximal primary cells and the centroid
  of distant primary cells in PCA space.

---

#### Relative single-centroid statistic
- `relative_original_dist_single_centroid`  
  The single-centroid distance statistic normalized by the total dispersion of all
  primary cells. This measures how dispersed proximal primary cells are relative to
  the overall spread of the primary cell population.

- `mean_relative_perm_dist_single_centroid`  
  Mean value of the relative single-centroid statistic across permutations, where
  proximity labels are randomly reassigned while preserving the number of proximal
  cells.

- `std_relative_perm_dist_single_centroid`  
  Standard deviation of the relative single-centroid statistic across permutations.

---

#### Z-scores from permutation-based null distributions
- `z_single`  
  Z-score of the observed single-centroid distance statistic relative to its
  permutation-based null distribution.

- `z_two`  
  Z-score of the observed two-centroid distance statistic relative to its
  permutation-based null distribution.

- `max_z_original`  
  The maximum of `z_single` and `z_two`. This is the main effect-size statistic
  summarizing dispersion differences for each (primary → neighbor) pair.

---

#### Gaussian-approximated max-Z statistic and p-values
- `max_z_gaus`  
  Z-score of `max_z_original` relative to the permutation-based distribution of
  max-Z values.

- `max_z_pval_gaus`  
  One-sided p-value computed from the standard normal distribution based on
  `max_z_gaus`.

- `pval_gaus_fdr_bh`  
  Benjamini–Hochberg FDR-adjusted p-value computed across all tested pairs.
  This value is used to define significant pairs (e.g. FDR ≤ 0.05).

---

#### Related output files
- `pca_pairs_significant.csv`  
  Subset of `pca_pairs_all_results.csv` containing only pairs that pass the
  FDR threshold.

- `pca_significant_pairs_panels.tif`  
  Multi-panel PCA visualization for significant pairs, showing the PCA embedding
  of primary cell type cells, proximal vs distant labeling, centroids, and summary
  statistics per pair.

---

### 2) Shuffled-proximity PCA results (sanity check)  
`pca_pairs_all_results_shuffled_proximity.csv`

This file contains results from a **sanity-check analysis** in which spatial proximity labels
are randomly shuffled across primary cell type cells.

#### Purpose
- To verify that significant dispersion signals observed in the main PCA analysis
  are **not driven by chance or by the geometry of the embedding alone**.
- Under shuffled proximity labels, true biological spatial structure should be destroyed.

#### Implementation
- The PCA embedding is identical to the main analysis.
- The number of proximal cells is preserved.
- Proximity labels are randomly reassigned across primary cells.
- All dispersion statistics and permutation tests are recomputed.

#### Interpretation
- The distribution of p-values in this file should be approximately uniform.
- Few or no pairs should pass FDR correction.

---

### 3) Marker-removal sensitivity analysis  
`pca_pairs_sig_marker_removal_sensitivity.csv`

This file contains results for **significant PCA pairs** after removing marker genes
associated with the **neighbor cell type**.

#### Purpose
- To reduce potential segmentation artifacts.

#### Implementation
- Applied only to pairs that were significant in the main PCA analysis.
- For each (primary → neighbor) pair:
  - Genes listed as markers of the neighbor cell type are removed.
  - PCA embedding and dispersion statistics are recomputed using the remaining genes.

---

### 4) t-SNE sensitivity analysis (significant pairs)  
`tsne_pairs_sig_sensitivity.csv`

This file contains results from repeating the dispersion analysis using
**t-SNE instead of PCA** for significant pairs.

#### Purpose
- To assess robustness of the results to the choice of dimensionality-reduction method.
- To test whether dispersion signals persist in a non-linear embedding.

#### Implementation
- Applied only to pairs significant in the PCA-based analysis.
- Gene expression is embedded in 2D using t-SNE.
- Proximity labels and dispersion statistics are computed as in the PCA analysis.

#### Interpretation
- Consistent results between PCA and t-SNE suggest robust spatial structure.
- Differences may reflect sensitivity to embedding geometry or non-linear effects.

---

### 5) t-SNE with shuffled proximity (all pairs)  
`tsne_pairs_all_sensitivity_shuffled_proximity.csv`

This file contains a **combined sanity-check and sensitivity analysis**:
t-SNE embedding with **shuffled proximity labels**, applied to all pairs (not only on previous significant pairs).

#### Purpose
- To test whether t-SNE alone can generate spurious dispersion signals
  when proximity information is randomized.
- To provide a conservative null model for non-linear embeddings.

#### Implementation
- t-SNE embedding is computed once per primary cell type.
- Proximity labels are randomly shuffled across primary cells.
- Dispersion statistics and permutation tests are computed for all pairs.

#### Interpretation
- Significant results are not expected in this file.

---

### 6) Summary counts  
`dispersion_pca_summary_counts.csv`

This file provides a compact overview of how many pairs are significant
at each stage of the analysis.

#### Contents
- Number of all tested PCA pairs
- Number of significant PCA pairs (after FDR)
- Number of significant pairs after shuffled proximity
- Number of significant pairs after marker removal
- Number of significant pairs in t-SNE analyses with or without shuffled proximity.

#### Purpose
- To facilitate quick comparison across analysis stages.
- To assess overall robustness and specificity of detected spatial dispersion signals.

---

## Example usage

### Run on all tissues and all cell types
```bash
insituprep dispersion-pca run \
  --summary-table-path data/summary_table \
  --distance-dir data/distance_matrix \
  --genes-names-path data/genes_names.txt \
  --dist-threshold 3.3 \
  --marker-genes-by-tissue-json data/marker_genes_by_tissue.json \
  --out results/dispersion_pca
```

### Restrict to selected cell types and tissue
```bash
insituprep dispersion-pca run \
  --summary-table-path data/summary_table \
  --distance-dir data/distance_matrix \
  --genes-names-path data/genes_names.txt \
  --primary '["T_CD3","Endothelial"]' \
  --neighbor '["Epithelial","B"]' \
  --tissue-ids '["100"]' \
  --marker-genes-by-tissue-json data/marker_genes_by_tissue.json \
  --out results/dispersion_pca_subset
```

---

## Interpretation notes

- This analysis tests **spatial structuring of expression programs**, not individual genes.
- A significant pair indicates that proximity to the neighbor cell type is associated
  with a systematic shift or dispersion in expression space.
- The max-type statistic provides robustness against metric-specific effects.

---

## Recommended usage

- Use PCA-based results for primary conclusions.
- Use shuffled-proximity, t-SNE and neighbor's marker removal analyses as **sanity checks**.
