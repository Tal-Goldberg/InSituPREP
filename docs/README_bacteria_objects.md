# Bacteria Objects Analysis

## Overview

This analysis identifies **spatial bacterial objects** based on clustered bacterial transcripts and tests whether **gene expression in nearby host cells is altered by proximity to these objects**.

Bacterial objects are defined as **spatial clusters of transcripts of a bacterial marker gene **. 
Each object does **not necessarily represent a single bacterium**, but rather indicates **local bacterial presence**.

After detecting significant bacterial objects, the analysis compares gene expression between:
- **Cells proximal to bacterial objects**
- **Cells distal from bacterial objects**

Differential expression is performed using **DESeq2 (via PyDESeq2)**.

Note: Although this analysis is described in the context of bacterial detection, the pipeline is fully generic.
The term bacteria can be replaced by any user-defined target gene, and the same workflow applies to spatial clustering of transcripts from any gene of interest.

The pipeline is implemented as an insituprep CLI command.

---

## Conceptual workflow

1. Load transcript-level spatial data.
2. Select transcripts belonging to a user-defined *target gene*.
3. Cluster target-gene transcripts in space to define bacterial objects.
4. Filter objects based on size and spatial properties.
5. Assign host cells as *near* or *far* from bacterial objects.
7. Run DESeq2 to identify genes associated with bacterial proximity to the host cells.
8. Save intermediate and final results.

---

## Input files

### 1. Transcript file (`--transcripts-path`) — required

Per-tissue CSV file with **single-transcript spatial coordinates**.

Required columns:

| Column name     | Description |
|-----------------|-------------|
| transcript_id   | Unique transcript identifier |
| gene_name       | Gene name |
| tissue          | Tissue ID (not necessary since it includes transcripts from only single tissue) |
| cell_id         | Cell ID the transcript assigned to (can be empty for unassigned transcripts) |
| x_micron        | X coordinate (µm). Global transcript coordinates aligned with cell centroids in the summary table (Input 2) |
| y_micron        | Y coordinate (µm). Global transcript coordinates aligned with cell centroids in the summary table (Input 2) |
| z_micron        | Z coordinate (µm). Global transcript coordinates aligned with cell centroids in the summary table (Input 2) |

Example:

| transcript_id | gene_name | tissue | cell_id   | x_micron | y_micron | z_micron |
|---------------|-----------|--------|-----------|----------|----------|----------|
| 1             | MYB       | 58     | 58.14.001 | 1591.85  | 179.39   | 4.24     |

---

### 2) Summary table file path (`--summary-table-path`) (required)

CSV with per-cell metadata and gene expression counts from all tissues.

**Requirements:**
- Required columns:
  - `Var1` (cell IDs)
  - `tissue` (string tissue identifier. Even if the CSV stores tissue as numeric, the tool loads and treats it as string)
  - `cell_type` (string)
- Should include gene expression columns (one column per gene, numeric raw counts)
- Required coordinates:
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
- This file is used as input to DESeq2.

---

### 3. Gene list (`--genes-names-path`) — required

A plain text file with **one gene name per line**.

This parameter is required. If --genes-names-path is not provided, the analysis cannot extract gene expression counts from the summary table file and will terminate with an error.
If `--genes-names-path` is not provided, the analysis runs object detection only and skips the DESeq2 step.

Notes:
- Only genes present both in this list and in the summary table file are analyzed.
- Genes missing from the summary table file (columns) are ignored with a warning.
- If no overlap exists, the analysis stops with an error.


Example (`genes_names.txt`):

ACTA2
AGR2
EPCAM
VIM
CD3D
CD3E

---

### Main parameters

| Parameter | Description |
|---------|-------------|
| `--transcripts-path` | Transcript-level spatial data |
| `--summary-table-path` | Cell-level counts + metadata |
| `--genes-names-path` | Gene names corresponding to expression columns in the summary table file, used to extract raw counts for DESeq2 |
| `--target-gene` | Gene name used to detect bacterial objects |
| `--out` | Output directory |

### Transcripts file column mapping

| Parameter | Description |
|----------|-------------|
| `--gene-col` | Column name in the transcripts CSV containing gene names |
| `--tx-x-col` | Column name for transcript X coordinate (µm). Only cells whose IDs appear in the transcripts table (column cell_id, via --tx-cell-id-col) are included in the downstream per-cell analyses (distance labeling and optional DESeq2). Cells present in the summary table but absent from the transcripts file are ignored. |
| `--tx-y-col` | Column name for transcript Y coordinate (µm) |
| `--tx-z-col` | Column name for transcript Z coordinate (µm). If not provided, clustering is performed in 2D |
| `--tx-cell-id-col` | Column name containing the cell_id associated with each transcript (empty/NA if not assigned to a cell) |

### Object detection parameters

| Parameter        | Description |
|------------------|-------------|
| `--min-samples`  | Minimum number of target-gene transcripts required to form a valid object |
| `--eps-max`      | Maximum DBSCAN radius (µm) evaluated during the parameter sweep |
| `--eps-min`      | Minimum DBSCAN radius (µm) evaluated during the parameter sweep |
| `--eps-final`    | Final DBSCAN radius (µm) used for object assignment. Use an integer value or "auto" to automatically select the eps that maximizes the global F1 score across the eps sweep |
| `--n-perms`      | Number of permutations used to compute p-values for transcript object significance |
| `--alpha`        | Significance threshold (p-value) used to define statistically significant transcript objects |
| `--seed` | Random seed for permutation testing, enabling reproducible object significance results |
| `--reuse-existing-objects` | Reuse previously computed object detection results (transcripts_with_objects.csv and eps_sweep_summary.csv) from the output directory and skip DBSCAN and permutation testing |


To guide the choice of the DBSCAN radius, the analysis scans a range of eps values between eps-min and eps-max.
This allows assessing the sensitivity of object detection to the clustering radius, before selecting a final value (eps-final) for the main analysis.

### Cells summary column mapping (summary table file)

| Parameter | Description |
|----------|-------------|
| `--celltype-col` | Column containing cell type annotations |
| `--cell-x-col` | Cell centroid X coordinate (µm) |
| `--cell-y-col` | Cell centroid Y coordinate (µm) |
| `--cell-z-col` | Cell centroid Z coordinate (µm). If not provided, distances are computed in 2D |

### Cell proximity parameters

| Parameter | Description |
|-----------|-------------|
| `--dist-threshold`      | Distance threshold (µm) used to assign cells to DESeq2 groups. Cells with minimum distance to any significant target-gene object ≤ threshold are labeled `near_sig`, otherwise `far_sig` |
| `--min-cells-per-group` | Minimum number of cells required in EACH group (near_sig and far_sig) per cell type in order to run DESeq2 |


---

## Step-by-step analysis details

### Step 1: Target gene selection
Only transcripts with:
```
gene_name == target_gene
```
are used to define bacterial objects.

---

### Step 2: Spatial clustering (DBSCAN)

Target-gene transcripts are clustered using **DBSCAN** in 3D/2D space.

Each cluster corresponds to a **candidate bacterial object**.

Computed per object:
- Number of transcripts
- Centroid (x, y, z)
- Max radius

Objects failing size or radius thresholds are discarded.

---

### Step 3: Cell–object distance calculation

For each host cell:
- Compute the minimum Euclidean distance to any transcript belonging to a
  statistically significant target-gene object.
- Assign label:
  - `near_sig` if the minimum distance ≤ `dist-threshold`
  - `far_sig` otherwise
- Cells containing at least one significant target-gene transcript are always
  assigned to `near_sig`.

---

### Step 4: Build DESeq2 count matrix

Cells are grouped into:
- `near_sig`
- `far_sig`

---

### Step 5: Differential expression (PyDESeq2)

DESeq2 model:
```
~ condition
```

Where `condition ∈ {near_sig, far_sig}`.

Statistics computed:
- Log2 fold change
- Wald test p-value
- Benjamini–Hochberg FDR

---
## Output files

All outputs are written to the directory specified by `--out-dir`.

---

### 1. DBSCAN eps sweep summary  
`eps_sweep_summary.csv`

This file summarizes the DBSCAN radius (`eps`) parameter sweep performed to detect
spatial transcript objects for the target gene.

Each row corresponds to a single `eps` value evaluated during the sweep.

| Column | Description |
|------|-------------|
| eps | DBSCAN radius value (µm) |
| n_target_transcripts | Total number of target-gene transcripts in the dataset |
| n_target_in_clusters | Number of target-gene transcripts assigned to clusters at this eps |
| recall | Fraction of target-gene transcripts assigned to clusters |
| n_clusters_total | Total number of DBSCAN clusters detected |
| n_clusters_significant | Number of clusters passing the statistical significance threshold |
| precision | Precision of cluster assignment (fraction of clustered transcripts belonging to significant clusters) |
| f1 | Global F1 score combining precision and recall |
| eps_selected | Final eps value selected for the analysis (used when `--eps-final auto`) |

This file is used to:
- Evaluate sensitivity of object detection to the DBSCAN radius
- Select the final eps value (`--eps-final`)
- Provide transparency and reproducibility for parameter selection

---

### 2. Transcript-level object assignments  
`transcripts_with_objects.csv`

This file contains **transcript-level annotations** resulting from the object
detection and significance testing.

Each row corresponds to a **single transcript** from the input transcripts file.
There is **no separate object-level table**; detected objects are represented
implicitly through transcript group assignments.

The file contains the original transcript columns, with additional object-related
annotations:

| Column | Description |
|------|-------------|
| transcript_id | Unique transcript identifier |
| gene_name | Gene name |
| group_id | Identifier of the detected transcript object (NA if the transcript is not assigned to any object) |
| p_value_group | Permutation-based p-value assessing the statistical significance of the transcript object (NA if not part of an object) |

Only transcripts belonging to statistically significant objects (based on the
specified `alpha` threshold) receive a valid `group_id` and `p_value_group`.

---

### 3. DESeq2 results (optional)  
`deseq2_results/`

This directory is created **only if** all DESeq2-related inputs are provided:
- `--summary-table-path`
- `--genes-names-path`
- Sufficient cells per group (`--min-cells-per-group`)

The directory contains one CSV file per cell type.

Each file corresponds to a **DESeq2 comparison for a single cell type**, comparing
cells that are near to vs far from statistically significant target-gene objects.

The results are generated using **PyDESeq2**, and therefore follow the standard
DESeq2 output format. For details on the statistical model and result columns,
please refer to the PyDESeq2 / DESeq2 documentation.

---


## Command-line arguments

```bash
insituprep bacteria-objects run \     
  --transcripts-path data/transcripts.csv \
  --out-dir output_dir \
  --summary-table-path data/summary_table.csv \
  --genes-names-path data/genes_names.txt \
  --target-gene BACTERIA \
  --gene-col gene_name \
  --tx-x-col x_micron \
  --tx-y-col y_micron \
  --tx-z-col z_micron \
  --eps-min 10 \
  --eps-max 150 \
  --min-samples 3 \
  --n-perms 10000 \
  --alpha 0.05 \
  --eps-final auto \
  --celltype-col cell_type \
  --cell-x-col X_space \
  --cell-y-col Y_space \
  --cell-z-col Z_space \
  --dist-threshold 10.0 \
  --min-cells-per-group 20
```
---

## Interpretation notes

- Significant genes indicate **host transcriptional responses associated with bacterial proximity**.
- Objects represent **spatial bacterial presence**, not single bacteria.
- Results depend on clustering and distance thresholds; sensitivity analysis is recommended.

---

## Recommended usage

- Start with conservative clustering parameters.
- Inspect `bacteria_objects.csv` visually.
- Validate DESeq2 assumptions before downstream interpretation.
