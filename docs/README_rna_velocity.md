# RNA Velocity (`insituprep rna-velocity`)

This document describes the **RNA velocity analysis** implemented in the `insituprep` package.
The workflow consists of **Stage 1 (a–d)** and **Stage 2**.

The pipeline is implemented as an insituprep CLI command.

To see all command-line options, run:

```bash
insituprep rna-velocity --help
```

---

## Conceptual overview: what this analysis does and why

RNA velocity uses **spliced (cytoplasmic)** and **unspliced (nuclear)** RNA counts to infer a **directional trend** in transcriptional state for each cell.
The analysis models transcriptional kinetics to compute per-gene splicing rates (γ), calculate RNA velocity (v = ds/dt), and predict future transcriptomic states s(t). The resulting dynamics are visualized as trajectories in PCA space, providing insight into how cells transition between states within the tissue.

--------------------------------------------------
BIOLOGICAL MODEL:

Let:
  u(t) = unspliced (nuclear) RNA as a function of time  
  s(t) = spliced (cytoplasmic) RNA as a function of time

For each gene in each cell:
  ds/dt = v = u(t) − γ·s(t)
  
where:
  γ = splicing rate constant, estimated using k-nearest neighbor (KNN) pooling.

Each cell’s initial state is defined as s(t=0) = s₀, and its predicted future state is:
  s(t) = s₀ + v·t

--------------------------------------------------

In this implementation:

- **Stage 1** learns per-gene kinetic slopes (**gamma**) from KNN-pooled spliced/unspliced counts, generates **phase portraits** for manual gene QC, learns an expression cutoff (`best_l`) for filtering genes, and computes velocity-derived **future spliced states** across a grid of abstract time points **t** (unitless).
- **Stage 2** takes the Stage 1 outputs (especially `spliced_values_in_multiple_t.pkl`) and:
  1) computes per-cell physical **proximity** to each cell type (continuous min-distance and binary thresholded proximity),
  2) builds simple **velocity geometry features** in PCA space (delta vectors, magnitudes, phases),
  3) tests **primary → neighbor** cell types associations between physical proximity and those velocity features (with optional permutation p-values + BH-FDR),
  4) optionally repeats association tests **per tissue**,
  5) optionally runs a **gene-level** analysis correlating per-gene |velocity| with proximity for significant primary→neighbor pairs.

> Important note about **t**  
> `t` is an **abstract, unitless time parameter** used to project cells forward along their inferred velocity direction.
> It is *not* a direct biological time scale. You choose a `t` by visual convergence of arrow plots (Stage 1-d),
> and Stage 2 then analyzes a selected `t` value.

---

## Input files

### Stage 1 inputs

#### 1) PCA origin file (`--pc-origin-path`) — required (Stage 1 + Stage 2)
CSV with PCA coordinates of each cell.

**Requirements:**
- Index must be **cell IDs**
- Required columns: `PC1`, `PC2`, `PC3`
- Optional column: `cell_type` (used only for coloring plots)

Typical columns:

| cell_id     | PC1        | PC2        | PC3        | cell_type (optional) |
|-------------|------------|------------|------------|----------------------|
| 100.26.007  | 0.0042664  | 2.8630893  | 1.7721467  | B                    |
| 982.14.016  | -2.8451111 | 5.3620301  | 0.0659634  | Endothelial          |
| 313.76.019  | -1.7862769 | 5.4585959  | -0.7030395 | Epithelial           |
| 100.1.056   | -2.0181643 | 3.8052376  | -0.0557922 | T_CD3                |


#### 2) Spliced counts (`--spliced-counts-path`) — required (Stage 1)
CSV with spliced (cytoplasmic) RNA expression matrix.

**Requirements:**
- Rows correspond to **cell IDs**
- The first column must be `tissue` (tissue ID)
- Followed by gene columns (cells × genes)

| cell_id      | tissue | ACTA2 | ACTG2 | ACTR3B | ADGRL4 | AGR2 | AHR |
|--------------|--------|-------|-------|--------|--------|------|-----|
| 100.26.007   | 100    | 2     | 0     | 0      | 0      | 0    | 0   |
| 982.14.016   | 982    | 0     | 0     | 0      | 0      | 0    | 0   |
| 313.76.019   | 313    | 0     | 0     | 0      | 0      | 0    | 0   |


#### 3) Unspliced counts (`--unspliced-counts-path`) — required (Stage 1)
CSV with unspliced (nuclear) RNA expression matrix.


| cell_id     | tissue | ACTA2 | ACTG2 | ACTR3B | ADGRL4 | AGR2 | AHR |
|-------------|--------|-------|-------|--------|--------|------|-----|
| 100.26.007   | 100    | 5     | 0     | 0      | 0      | 0    | 0   |
| 982.14.016   | 982    | 0     | 4     | 0      | 14     | 0    | 0   |
| 313.76.019   | 313    | 0     | 0     | 0      | 0      | 0    | 7   |

**Requirements:**
- Same structure as spliced CSV:
  - rows = cells
  - first column = `tissue`
  - same gene columns (cells × genes)

---

### Stage 1 curated input (for Stage 1-c and Stage 1-d)

#### 4) Ground truth velocity genes (`--ground-truth-path`) — required (Stage 1-c and Stage 1-d)
Plain text file: **`ground_truth_genes.txt`**, one gene per line.

These genes are the subset the user decided to **keep** for downstream velocity calculations after inspecting the phase portraits
(see Stage 1-b). In practice:

- You inspect gamma/phase portrait plots
- You keep genes whose spliced–unspliced geometry matches expected RNA velocity behavior
  (e.g., characteristic elliptical structure / clear linear trend around the fitted line)

Example (`ground_truth_genes.txt`):

```
ACTA2
VIM
EPCAM
CD3D
...
```

---

### Stage 2 additional inputs

### 5) Summary table file path (`--summary-table-path`) (required)

CSV with per-cell metadata and gene expression counts from all tissues.

**Requirements:**
- Required columns:
  - `Var1` (cell IDs. The tool sets 'Var1' as the DataFrame index) 
  - `tissue` (string tissue identifier. Even if the CSV stores tissue as numeric, the tool loads and treats it as string))
  - `cell_type` (string)
- May include gene expression columns (one column per gene, numeric raw counts) - these columns are ignored by RNA velocity analysis (since we use the 'spliced' and 'unspliced' counts matrices instead).
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
- Coordinates must be numeric and in micrometers (µm).


#### 6) Distance matrices folder (`--distance-dir`) — optional (Stage 2)

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

If `--distance-dir` is **not provided**, proximity is computed from `X_space`, `Y_space` (and optionally `Z_space`) in the summary table file.

#### 7) Stage 1 pickle dictionary (`--gene-dict-pkl`) — required (Stage 2)
Pickle produced by Stage 1, typically:

- `spliced_values_in_multiple_t.pkl`

Stage 2 expects this pickle to contain (for the selected `t`) at least:
- `st_norm_filtered` — normalized predicted spliced expression for the future state s(t)
- `filtered_v` — filtered per-gene velocities (used for the gene-level analysis)

---

## Stage-by-stage: what happens in each step

## Stage 1-a: Gamma grid computation (`insituprep rna-velocity stage1-a`)

**Goal:** scan a grid of `(k_neighbors, quantile_outliers)` values and compute per-gene **gamma** (slope) estimates.

What happens:
- Compute (or reuse) KNN relationships (in PCA space) between cells (saved as `neighbors.csv`)
- For each `k` in `--k-values` and each `q` in `--quantiles`:
  - Pool spliced/unspliced counts across KNN
  - Fit a regression (with outlier handling via `quantile_outliers`)
  - Store per-gene gamma values

Key parameters (from CLI):
- `--k-values`: JSON list of K values (ints), e.g. `"[5,10,20,30]"`
- `--quantiles`: JSON list of outlier quantiles, e.g. `"[0.05,0.075,0.1]"`

---

## Stage 1-b: Gamma portrait plots (`insituprep rna-velocity stage1-b`)

**Goal:** generate per-gene **phase portraits** (gamma plots) for a **chosen** `(k, quantile)`.

What happens:
- Load/reuse `neighbors.csv`
- For each gene:
  - Plot spliced vs unspliced (KNN-pooled)
  - Overlay fitted gamma line
  - Optionally color points by `cell_type` if provided in `pc_origin`

These plots are used for **manual QC** and to build `ground_truth_genes.txt`.

Key parameters:
- `--k`: the chosen K (neighbors) based on the previous step
- `--quantile`: the chosen outlier quantile based on the previous step.

---

## Stage 1-c: Learn expression threshold `best_l` (`insituprep rna-velocity stage1-c`)

**Goal:** choose a global normalized-expression cutoff (`best_l`) by maximizing an F-score based on your curated ground truth genes.

What happens:
- Load ground truth gene list (`ground_truth_genes.txt`)
- Scan candidate thresholds `l`
- For each threshold:
  - mark genes as "kept" vs "filtered out" based on expression normalization rules
  - compute precision/recall against the ground truth
- Select `best_l` = threshold that maximizes F-score

Output includes a precision/recall table and a summary plot.

---

## Stage 1-d: Compute s(t) across multiple t (`insituprep rna-velocity stage1-d`)

**Goal:** compute velocity-derived future spliced states **s(t)** across a grid of `t` values, and write the main pickle output.

What happens:
- Uses chosen `k`, `quantile`, your curated gene list, and `best_l`
- For a grid of `t` values:
  - compute filtered velocities
  - compute predicted future spliced expression `s(t)`
  - sample `--n-cells` cells and generate PCA arrow plots for each `t`

**How to choose `t_final`:**
- Inspect arrow plots across increasing `t`
- Choose the point where arrows **converge** (do not change drastically from one `t` to the next)

Stage 2 will then analyze one selected `t` (via `--t`).

---

## Stage 2: Proximity ↔ velocity-geometry association analysis (`insituprep rna-velocity stage2`)

Stage 2 uses Stage 1 outputs (for a selected `t`) and tests whether **physical proximity to a neighbor cell type**
is associated with **velocity geometry features** in PCA space.

### 2.1 Proximity tables
Stage 2 builds two proximity tables:

- **Continuous**: for each cell and each cell type, the **minimum distance** to any cell of that type
- **Binary**: thresholded proximity (1 if min-distance ≤ `--distance-threshold`, else 0)

Proximity is computed either from:
- distance matrices (`--distance-dir`) *or*
- coordinates (`X_space`, `Y_space`, optional `Z_space`) if `--distance-dir` is not provided

### 2.2 End-state PCA and delta vectors
For the selected `t`, Stage 2:
- loads `st_norm_filtered` (normalized future-state spliced expression matrix)
- performs PCA on that matrix to obtain end-state coordinates (`PC1_end`, `PC2_end`, …)
- merges these with the origin PCs (`PC1`, `PC2`, …)
- computes delta vectors in PCA space:
  - `delta_x = PC1_end - PC1`
  - `delta_y = PC2_end - PC2`

### 2.3 Velocity geometry features (y-features)
Stage 2 derives four per-cell features:

- `magnitude_delta` : length of the delta vector in PCA space  
- `magnitude_s0`    : distance of the origin point from (0,0) in PC1/PC2  
- `phase_deg_delta` : angle (degrees) of the delta vector  
- `phase_deg_s0`    : angle (degrees) of the origin PC1/PC2 point  

### 2.4 Pairwise primary→neighbor association tests
For each **primary cell type** and each **neighbor cell type**:

- Use only cells of the primary type
- Take x = proximity(primary cell → neighbor type)  
- Take y = chosen velocity feature (one of the four y-features)
- Compute correlation:
  - continuous proximity → Pearson correlation
  - binary proximity → point-biserial correlation

Optional permutation testing (enabled by default):
- permute x among primary cells `--n-perm` times
- compute an empirical permutation p-value
- apply BH-FDR

### 2.5 Tissue-stratified analysis (optional)
If enabled (`--compute-tissue`), Stage 2 repeats the primary→neighbor tests **within each tissue** separately.

### 2.6 Gene-level analysis for significant pairs (optional)
If enabled (`--compute-gene-level`), Stage 2:
- selects significant primary→neighbor pairs from one chosen feature table (default `magnitude_delta`)
- for each selected pair:
  - uses `filtered_v` (per-cell, per-gene velocity matrix from Stage 1)
  - correlates |velocity_gene| with proximity across primary cells
  - optionally computes gene-level permutation p-values (can be very slow)

---

## Command-line parameters

Run `--help` to see defaults and all options:

```bash
insituprep rna-velocity stage1-a --help
insituprep rna-velocity stage1-b --help
insituprep rna-velocity stage1-c --help
insituprep rna-velocity stage1-d --help
insituprep rna-velocity stage2 --help
```

Below is a readable summary of the most important parameters.

### Common
- `--pc-origin-path`  
  PCA embedding CSV. Index = cell IDs. Required columns: `PC1`, `PC2`, `PC3` (for Stage 1).

- `--out`  
  Output directory.

### Stage 1-a
- `--k-values`  
  JSON list of KNN sizes to scan (ints).
- `--quantiles`  
  JSON list of quantile outliers to scan (floats).

### Stage 1-b / 1-c / 1-d
- `--k`  
  Chosen number of PCA nearest neighbors for KNN pooling.
- `--quantile`  
  Chosen outlier quantile used for gamma fitting.
- `--ground-truth-path`  
  Path to `ground_truth_genes.txt` (one gene per line).
- `--best-l` *(Stage 1-d)*  
  Normalized-expression cutoff chosen in Stage 1-c.
- `--n-cells` *(Stage 1-d)*  
  Number of cells to sample for arrow plots.

### Stage 2
- `--summary-table-path`  
  Per-cell metadata + expression.
- `--distance-dir`  
  Folder with `distance_matrix_<tissue>.csv` files. If omitted, use X/Y(/Z) columns in summary table.
- `--gene-dict-pkl`  
  Stage 1 pickle (typically `spliced_values_in_multiple_t.pkl`) includes the calculated future spliced counts and genes velocities.
- `--t`  
  Which timepoint to analyze (must exist as a key in the pickle dict).
- `--proximity-mode`  
  `continuous` or `binary` (which proximity table is used in association tests).
- `--distance-threshold`  
  Threshold (µm) for building the binary proximity table.
- `--maximum-distance-threshold`  
  For coordinate-based proximity only: distances above this are set to NaN and excluded (reducing calculation time when generating proximity tables)
- `--min-cells`  
  Minimum number of primary cells required to test a pair.
- `--n-perm`  
  Number of permutations for pairwise permutation testing.
- `--compute-pair-permutations/--no-compute-pair-permutations`  
  Enable/disable pairwise permutations.
- `--compute-tissue/--no-compute-tissue`  
  Enable/disable tissue-stratified analysis.
- `--compute-gene-level/--no-compute-gene-level`  
  Enable/disable gene-level analysis.
- `--run-gene-permutations/--no-run-gene-permutations`  
  Enable/disable gene-level permutation p-values (very slow).

---

## Output files

All outputs are written under the directory specified by `--out`.

### Stage 1 outputs

#### `neighbors.csv`
KNN neighbor indices saved to avoid recomputing neighbors across Stage 1 steps.

#### `gamma_per_gene_per_quantile_per_k_dict.pkl`
Pickle containing per-gene gamma estimates across the scanned `(k, quantile)` grid.

#### `gamma_plots/`
Directory of per-gene phase portraits generated in Stage 1-b.
File naming pattern includes the gene name and K (e.g. `{GENE}_K{k}.tif`).

#### `precision_recall_unspliced.csv`
Table produced in Stage 1-c containing precision/recall/F-score across scanned `l` thresholds.

#### `Precision_Recall_F-score_plot.tif`
Plot summarizing Stage 1-c threshold selection.

#### `stage1_selected_params.txt`
Text summary of selected parameters (e.g., chosen k/quantile/best_l).

#### `spliced_values_in_multiple_t.pkl`  **(critical / required for Stage 2)**
Main Stage 1 output used by Stage 2.

This pickle stores a dictionary keyed by `t` values. For each `t`, it contains (at minimum):
- `st_norm_filtered`: predicted future spliced expression matrix for s(t) (cells × genes)
- `filtered_v`: filtered velocity matrix (cells × genes)

#### `t_{t}-{cell_type}.tif`
Arrow plots generated in Stage 1-d showing origin PCA locations and projected endpoints for each `t`.

---

### Stage 2 outputs

Stage 2 writes outputs under:

- `stage2_inputs/`  (inputs computed by Stage 2, mainly proximity tables)
- `stage2_outputs/` (association results)

#### Stage 2 input tables (`stage2_inputs/`)

- `proximity_to_cell_types_continuous.csv`  
  Continuous min-distance to each cell type (cells × cell_types), plus a `cell_type` column.

- `proximity_to_cell_types_binary_thr{distance_threshold}.csv`  
  Binary proximity (0/1) using `--distance-threshold`.

#### Stage 2 merged per-cell table (`stage2_outputs/`)

- `cells_info_{proximity_mode}.csv`  
  Per-cell merged table including:
  - origin PCs (from `pc_origin`)
  - end PCs (`PC1_end`, `PC2_end`, …) from PCA of `st_norm_filtered`
  - delta vectors (`delta_x`, `delta_y`)
  - magnitudes + phases
  - proximity columns to each neighbor cell type

This file is the easiest place to debug and understand what Stage 2 is using per cell.

#### Stage 2 pairwise association tables (`stage2_outputs/`)

For each y-feature:

- `pairwise_magnitude_delta_{proximity_mode}.csv`
- `pairwise_magnitude_s0_{proximity_mode}.csv`
- `pairwise_phase_deg_delta_{proximity_mode}.csv`
- `pairwise_phase_deg_s0_{proximity_mode}.csv`

Each table contains rows like:

- `primary_cell_type`
- `neighbor_cell_type`
- `correlation`
- `p_value`
- `p_value_fdr`
- (if permutations enabled) permutation-based p-values and FDR columns

#### Stage 2 tissue-stratified table (`stage2_outputs/`)

- `pairwise_by_tissue_{proximity_mode}.csv`  
  Same idea as pairwise tables, but computed within each tissue.

#### Stage 2 gene-level outputs (`stage2_outputs/gene_level/`)

For each significant primary→neighbor pair:

- `gene_velocity_vs_proximity__{primary}__{neighbor}__{proximity_mode}.csv`

Contains, for each gene:
- association between |velocity_gene| and proximity
- (optional) gene-level permutation p-values and BH-FDR

---

## Example usage

### Stage 1-a: gamma grid scan
```bash
insituprep rna-velocity stage1-a \
  --pc-origin-path data/pc_origin.csv \
  --spliced-counts-path data/spliced_counts.csv \
  --unspliced-counts-path data/unspliced_counts.csv \
  --k-values "[10,20,30,40]" \
  --quantiles "[0.05,0.075,0.1]" \
  --out results/rna_velocity
```

### Stage 1-b: gamma plots for chosen parameters
```bash
insituprep rna-velocity stage1-b \
  --pc-origin-path data/pc_origin.csv \
  --spliced-counts-path data/spliced_counts.csv \
  --unspliced-counts-path data/unspliced_counts.csv \
  --k 30 \
  --quantile 0.075 \
  --out results/rna_velocity
```

### Stage 1-c: choose best_l using ground truth genes
```bash
insituprep rna-velocity stage1-c \
  --pc-origin-path data/pc_origin.csv \
  --spliced-counts-path data/spliced_counts.csv \
  --unspliced-counts-path data/unspliced_counts.csv \
  --k 30 \
  --quantile 0.075 \
  --ground-truth-path data/ground_truth_genes.txt \
  --out results/rna_velocity
```

### Stage 1-d: compute s(t) for multiple t and generate arrow plots
```bash
insituprep rna-velocity stage1-d \
  --pc-origin-path data/pc_origin.csv \
  --spliced-counts-path data/spliced_counts.csv \
  --unspliced-counts-path data/unspliced_counts.csv \
  --k 30 \
  --quantile 0.075 \
  --ground-truth-path data/ground_truth_genes.txt \
  --best-l 0.038 \
  --n-cells 500 \
  --out results/rna_velocity
```

### Stage 2: run proximity↔velocity association analysis for one t
```bash
insituprep rna-velocity stage2 \
  --pc-origin-path data/pc_origin.csv \
  --summary-table-path data/sammary_table.csv \
  --distance-dir data/distance_matrix \
  --gene-dict-pkl results/rna_velocity/spliced_values_in_multiple_t.pkl \
  --t 3.0 \
  --proximity-mode continuous \
  --distance-threshold 1.0 \
  --n-perm 10000 \
  --out results/rna_velocity_stage2
```

---

## Interpretation notes

- Stage 1 is intentionally interactive: gene selection and `t` selection are based on diagnostic plots.
- Stage 2 does **not** recompute RNA velocity; it **analyzes** Stage 1’s inferred future-state geometry and its association with spatial proximity.
- Significant primary→neighbor results in Stage 2 mean:
  > among cells of the primary type, proximity to the neighbor type is associated with velocity-derived geometry features
  (e.g., larger delta magnitudes or consistent delta directions).

---

## Recommended usage

1. **Stage 1-a**: scan a reasonable grid of K and quantiles.
2. **Stage 1-b**: inspect phase portraits; curate `ground_truth_genes.txt`.
3. **Stage 1-c**: compute `best_l`.
4. **Stage 1-d**: inspect arrow plots across t; decide a stable `t` (convergence).
5. **Stage 2**: pick that `t`, choose proximity mode (continuous first is often more informative), and run pairwise + (optional) gene-level analysis.
