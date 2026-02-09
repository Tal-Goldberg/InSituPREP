# Comparison Across Platforms and Patients (`insituprep platform-comparison`)

This document describes the **cross-platform / cross-patient comparison analysis** implemented in the `insituprep` package
(via `insituprep platform-comparison run`).

The goal of this analysis is to **compare proximity-based differential-expression (DE) signatures across platforms**
(e.g. ExSeq vs MERFISH) and across patient/tissue samples, using the **DESeq-proximity** outputs (`final_genes_*.csv`)
as input.

Conceptually, for each cell-type direction, the analysis builds a **tissue-by-gene matrix of p-values** derived from
DESeq-proximity results, embeds tissues in a 2D PCA space, and then tests whether:
1) **Matched platform pairs** (e.g. `330` vs `MERFISH_330`) are **closer than expected by chance** in PCA space, and
2) **specific pathological subgroup (e.g. HR+/HER2-) tissues** cluster more tightly than expected by chance.

  
The method integrates:
- Reading `final_genes_*` tables produced by `insituprep deseq-proximity`
- Building p-value matrices (genes × tissues) for each cell-type direction
- PCA embedding of tissues based on `-log10(p-value)` patterns across genes
- Permutation-based p-values for:
  - PCA explained-variance ratios (PC1 and PC2)
  - Cross-platform tissue-pair distances in PCA space
  - Specific pathological subgroup (e.g. Receptor-status (HR+/HER2-)) proximity in PCA space
- Publication-ready PCA scatter plots per tested cell-type direction

The pipeline is implemented as an insituprep CLI command.

To see all command-line options, run:

```bash
insituprep platform-comparison --help
```

---

## Conceptual overview: what this analysis does and why

The analysis runs in **two directions** for each non-tumor cell type:

1) **NotTumor → Tumor**  
   Primary = each non-tumor cell type (one at a time)  
   Neighbor = tumor cell type (fixed)

2) **Tumor → NotTumor**  
   Primary = tumor cell type (fixed)  
   Neighbor = each non-tumor cell type (one at a time)

For each direction (i.e., for each `(primary, neighbor)` pair), the pipeline does:

### Step 1 — Load DESeq p-values across tissues
For the requested `(primary, neighbor)` cell-type pair, the code looks under `--deseq-results-dir` for files matching:

```
final_genes_{tissue}_{primary_cell_type}_vs_{neighbor_cell_type}_LFC*.csv
```

From each matching file, it loads:
- gene names (default column: `gene`)
- raw p-values (default column: `pvalue`)

The result is a **p-value matrix**:
- rows = genes
- columns = tissues
- join strategy = **inner join across tissues** (only genes present in all loaded tissues remain)

> Important: this analysis compares *patterns across genes*, so it requires overlap of genes across tissues for each pair.

### Step 2 — Transform p-values and build the tissue embedding
- Replace p-values of 0 (if any) with a small positive value (to avoid `-log10(0)`).
- Convert to `-log10(pvalue)` so that stronger signals correspond to larger values.
- Standardize the tissue-by-gene matrix with `StandardScaler`.
- Fit PCA with 2 components and record:
  - PC coordinates for each tissue (PC1, PC2)
  - explained variance ratios for PC1 and PC2

### Step 3 — Validate that PCA structure is non-random (PC permutation test)
To avoid interpreting PCA structure that could arise by chance, the analysis computes permutation p-values for PC1/PC2:

- It permutes the standardized matrix **by shuffling values within each gene column** across tissues.
- For each permutation, PCA is refit and explained variance ratios are recorded.
- PC p-values are computed as:

```
p(PCk) = P( perm_explained_variance_ratio_k >= real_explained_variance_ratio_k )
```

**Filtering rule:** a pair is kept only if **both** PCs pass:

- `PC1 p-value <= pc_pvalue_alpha` AND `PC2 p-value <= pc_pvalue_alpha`  
(default `pc_pvalue_alpha = 0.05`)

If this fails, the pair is skipped (no plot and no row in the final summary).

### Step 4 — Cross-platform similarity (ExSeq vs MERFISH) for matched tissue IDs
Tissues are treated as **paired** when they share the same numeric suffix:

- `330` (ExSeq) pairs with `MERFISH_330`

For each paired tissue ID, the analysis computes:
- the observed Euclidean distance between the two platform points in PCA space
- a permutation null distribution of distances, built by permuting the tissue-by-gene matrix and refitting PCA each time

The distance p-value is computed as:

```
p_dist = P( perm_distance <= observed_distance )
```

**Interpretation:**  
- **Small p_dist** ⇒ the matched ExSeq vs MERFISH tissues are **unusually close** in PCA space (higher concordance)  
- **Large p_dist** ⇒ no evidence of unusually high concordance (distance is not smaller than expected)

> Note: internally, PCA coordinates from each permutation are rescaled back to the original PC1/PC2 ranges
> so distances are comparable across permutations.

### Step 5 — Specific pathological subgroup similarity (Receptor-status (HR+/HER2-)) among tissues
Optionally, the analysis also tests whether ExSeq tissues labeled HR+/HER2- cluster more tightly than expected.

- Tissues are first filtered to “signal-present” tissues: tissues with **min p-value ≤ `min_pv_for_tissue_filter`**
  (default `min_pv_for_tissue_filter = 0.01`) in the loaded `final_genes` table.
- Among these, it selects **ExSeq (non-prefixed) tissues** whose tissue name is listed in `--hr-positive-tissues`.
- It computes the **mean pairwise distance** among the selected ExSeq HR+/HER2- tissues in PCA space.
- A permutation p-value is computed as:

```
p_receptor = P( perm_mean_distance <= observed_mean_distance )
```

**Interpretation:**
- **Small p_receptor** ⇒ HR+/HER2- ExSeq tissues are **unusually close** (more similar) in this pair’s PCA embedding
- **Large p_receptor** ⇒ no evidence of unusually tight clustering

If there are ≤ 1 HR+/HER2- ExSeq tissues after filtering, this test is skipped and `p_receptor` is reported as `NaN`.

### Step 6 — Plotting
For each kept `(primary, neighbor)` pair, the analysis generates a PCA scatter plot:
- Points = tissues with min p-value ≤ `min_pv_for_tissue_filter`
- Colors = matched across platforms by numeric tissue ID (e.g. `330` and `MERFISH_330` share the same color)
- Shapes = receptor status
  - `*` for HR+/HER2-
  - `o` for all other tissues

The plot title is:

```
PCA of DeSeq P-Values: <Primary> proximal to <Neighbor>
```

---

## Input files

### 1) DESeq results directory (`--deseq-results-dir`) — required

A directory containing per-tissue-and-cell types pair CSVs produced by `insituprep deseq-proximity`, named like:

```
final_genes_{tissue}_{primary_cell_type}_vs_{neighbor_cell_type}_LFC_*.csv
```

Example filenames:
- `final_genes_330_Endothelial_vs_Epithelial_LFC_0.30.csv`
- `final_genes_MERFISH_330_Endothelial_vs_Epithelial_LFC_0.30.csv`

| gene    | baseMean | log2FoldChange | lfcSE | stat.   | pvalue     | padj    | ......  | primary_cell_type | neighbor_cell_type | tissue |
|---------|----------|----------------|-------|---------|------------|---------|---------|-------------------|--------------------|--------|
| CCND1   | 2.2667   | -0.4339        | 0.2002| -2.1667 | 0.03026    | 0.07158 | ......  | Endothelial       | Epithelial         | 330    |
| CD24    | 52.0225  | -0.5355        | 0.2349| -2.2796 | 0.02263    | 0.06099 | ......  | Endothelial       | Epithelial         | 330    |
| CD69    | 3.3015   | 0.6759         | 0.2176| 3.1063  | 0.00189    | 0.00838 | ......  | Endothelial       | Epithelial         | 330    |
| CDK6    | 2.0265   | 0.6762         | 0.1735| 3.8967  | 9.75e-05   | 0.00073 | ......  | Endothelial       | Epithelial         | 330    |
| CDKN2A  | 3.0544   | -1.2763        | 0.2923| -4.3666 | 1.26e-05   | 0.00012 | ......  | Endothelial       | Epithelial         | 330    |



**File requirements:**
- Must include a gene column (default: `gene`)
- Must include a p-value column (default: `pvalue`)
- May include additional DESeq columns (ignored by this analysis)

**Important behavior:**
- For each tissue, the analysis picks the **first** matching file for that prefix (deterministic due to sorting).
- If no matching file is found for a tissue, that tissue is skipped for the current pair.
- Genes are combined across tissues using an **inner join** (only shared genes remain).

---

## Command-line parameters

### Required
- `--deseq-results-dir PATH`  
  Directory with `final_genes_*.csv` outputs from `insituprep deseq-proximity`.

- `--out PATH`  
  Output directory for this analysis results.

- `--tissues TEXT`  
  JSON list of tissue IDs (strings). Must include both platforms if you want cross-platform pairing (note that for each platform its own `final_genes_*.csv` is used).  
  Example:
  ```bash
  --tissues '["364","330","MERFISH_330","MERFISH_364"]'
  ```

- `--non-tumor-cell-types TEXT`  
  JSON list of non-tumor cell types to test against the tumor cell type (in both directions).  
  Example:
  ```bash
  --non-tumor-cell-types '["B","Endothelial","T_CD3"]'
  ```

- `--tumor-cell-type TEXT`  
  Tumor cell type name **exactly as used in the DESeq `final_genes` filenames**.  
  Example:
  ```bash
  --tumor-cell-type "Epithelial"
  ```

### Optional (defaults shown)
- `--hr-positive-tissues TEXT`   
  JSON list of tissue IDs corresponding to a specific pathological subgroup (e.g. HR+/HER2-).
  These tissues are used to test whether samples sharing this pathological status cluster more tightly than expected by chance in PCA space.  
  These are used for the receptor-status proximity test and for plot markers (`*`).  
  Example:
  ```bash
  --hr-positive-tissues '["982","880","330","783","MERFISH_982","MERFISH_880"]'
  ```

- `--platform-prefix TEXT` (default: `"MERFISH_"`)  
  Prefix used to identify the second platform tissues.  
  Tissues are paired by numeric suffix:
  - `330` ↔ `MERFISH_330`

- `--n-perm-pca INT` (default: `10000`)  
  Number of permutations for PC explained-variance p-values.

- `--n-perm-dist INT` (default: `10000`)  
  Number of permutations for platform-distance and receptor-proximity p-values.

- `--rng-seed INT` (default: `None`)  
  Random seed for reproducibility. If omitted, permutations are stochastic."

- `--verbose / --no-verbose` (default: verbose)  
  Print progress messages during the run.

- `--perm-progress-every INT` (default: `1000`)  
  When verbose: print permutation progress every N permutations. `0` disables progress printing.

---

## Output files

All outputs are written under the directory specified by `--out`.

### 1) Summary table
**Filename**
- `comparison_across_platforms_and_patients_summary.csv`

**Meaning**  
One row per tested cell-type direction that passed the PC permutation filter.

**Columns**
- `Primary cell type`  
  The primary cell type for this row.

- `Neighbor cell type`  
  The neighbor cell type for this row.

- `Direction`  
  Either:
  - `NotTumor→Tumor`
  - `Tumor→NotTumor`

- `PC1 p-value`, `PC2 p-value`  
  Permutation p-values testing whether each PC’s explained variance ratio is larger than expected by chance.
  (Rows are included only if both pass `pc_pvalue_alpha`.)

- `ExSeq-MERFISH similarity p-value`  
  A per-tissue pairing p-values:
  - computed for matched numeric tissue IDs that exist in both platforms
  - **smaller is more similar** (unusually small PCA distance)

  Example format:
  ```
  Tissue 330: p-value=0.003 ; Tissue 364: p-value=0.12
  ```

- `Receptor status similarity p-value`  
  Permutation p-value for whether ExSeq HR+/HER2- tissues are unusually close (mean distance) in PCA space.
  `NaN` if not enough tissues were available for this test.

- `Plot path`  
  Path to the saved PCA plot (`.tif`) for this pair.

### 2) PCA plots
Plots are saved under:

- `not_tumor_to_tumor/plots/`
- `tumor_to_not_tumor/plots/`

**Filename pattern**
- `pca_<direction>_primary_<Primary>_neighbor_<Neighbor>.tif`

Each plot shows:
- PC1/PC2 coordinates for tissues passing the min-pvalue tissue filter
- colors consistent between paired platforms (same numeric tissue ID)
- shape indicating receptor status (HR+/HER2- vs Other)
- a legend showing each tissue and its `(min_pvalue, max_pvalue)` within the loaded `final_genes` table for that tissue

> Optional dependency: If `adjustText` is installed (`pip install adjustText`), plot labels are de-overlapped.
> If not installed, a warning is printed and labels may overlap.

---

## Example usage

### Run on a mixed ExSeq + MERFISH tissue set

```bash
insituprep platform-comparison run \
  --deseq-results-dir /path/to/deseq_proximity_results \
  --out /path/to/output/platform_comparison \
  --tissues '["330","364","783","880","982","MERFISH_330","MERFISH_364","MERFISH_880","MERFISH_982"]' \
  --non-tumor-cell-types '["B","Endothelial","T_CD3","Fibroblast"]' \
  --tumor-cell-type "Epithelial" \
  --hr-positive-tissues '["982","880","330","783","MERFISH_982","MERFISH_880"]' \
  --n-perm-pca 10000 \
  --n-perm-dist 10000 \
  --rng-seed 0
```
---

## Interpretation notes

- This analysis **does not** re-run DESeq; it **reuses DESeq-proximity outputs**.
- A “signature” here means the **pattern of DESeq p-values across genes** for a given (primary, neighbor) definition.
- **Small cross-platform p-values** (distance test) suggest that ExSeq vs MERFISH tissues form a **consistent signature** in PCA space.
- The receptor-status p-value is an additional patient-level consistency test (optional) focused on HR+/HER2- ExSeq tissues.

---
