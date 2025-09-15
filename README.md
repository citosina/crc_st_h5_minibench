# crc_st_h5_minibench

This mini-benchmark reproduces a simplified **Visium CRC** analysis using a precomputed 10x/Space Ranger count matrix (`filtered_feature_bc_matrix.h5`) derived from:

> *Spatial transcriptomics reveals tumor–immune architecture and signaling in colorectal cancer.*  
> Processed Visium outputs hosted on Zenodo (see Data Sources).

The task focuses on running lightweight, fully local scripts that compute spot-level QC and naïve marker summaries **starting from the H5 matrix**. No alignment or counting is performed here.

---

## 1) Data sources

- **Visium processed output (H5)**: `filtered_feature_bc_matrix.h5` from Zenodo record DOI `10.5281/zenodo.7551712`.  
  We used the CRC section **SN84_A120838_Rep1** in testing.  
  *(Only the H5 is required for this mini-benchmark.)*

- **This repository ships empty folders and expects you to place the H5 at:**  
  `sample_submission_boilerplate/data/filtered_feature_bc_matrix.h5`

  The verification step will create a symlink under `sample_submission_boilerplate/results/filtered_feature_bc_matrix.h5` (and under `outs/` for convenience).

> **Note:** The H5 is *not* committed to the repo. You must source it externally (e.g., from Zenodo or your internal data share) and drop it into `sample_submission_boilerplate/data/` before running the workflow.

---

## 2) How to place the data

Place the single file:

```
sample_submission_boilerplate/data/
  filtered_feature_bc_matrix.h5
```

Then run the verify step (below).

---

## 3) Workflow steps

All scripts live in `sample_submission_boilerplate/workflow/` and write results under `sample_submission_boilerplate/results/`.

### Step 0 – Verify H5 and prepare links
Checks that the file is readable and symlinks it into `results/` (and `outs/`).

```bash
cd sample_submission_boilerplate
bash workflow/00_verify_h5.sh
```

**Creates/links:**

```
results/filtered_feature_bc_matrix.h5   -> data/filtered_feature_bc_matrix.h5
outs/filtered_feature_bc_matrix.h5      -> data/filtered_feature_bc_matrix.h5
```

### Step 1 – Metrics
Computes basic spot metrics (counts per spot, detected genes per spot, medians/means) from the H5.

```bash
bash workflow/03_metrics.sh
```

**Outputs:**

```
results/metrics/spot_metrics.tsv   # “metric<TAB>value” table
results/metrics/top_genes.txt      # top 5 genes by total counts (one per line)
```

### Step 2 – Markers
Ranks genes by total counts across spots and exports a simple top list.
Also includes an epithelial CMS marker panel used for a toy enrichment call.

```bash
bash workflow/04_markers.sh
```

**Outputs:**

```
results/markers/marker_scores.tsv  # “gene<TAB>score” (top 500 by default)
```

### Step 3 – Evaluate
Parses the metrics and markers, assembles the answers (JSON).

```bash
python workflow/05_evaluate.py
```

**Output:**

```
results/results.json
```

Example content (values will reflect your H5):

```json
{
  "q1": 99.3,
  "q2": 328,
  "q3": 3958,
  "q4": ["MT-CO2", "MT-CO3", "RPL13", "RPS21", "MT-ATP6"],
  "q5": "CMS2_like"
}
```

---

## 4) Questions & Answers

- `sample_submission_boilerplate/questions.yaml` – defines the five benchmark questions that are answerable **from the H5**.
- `sample_submission_boilerplate/answers.yaml` – canonical answers used for validation when testing this template.

These mirror the analysis performed by the workflow scripts above.

---

## 5) Environment

A minimal, reproducible conda/mamba environment is defined in `sample_submission_boilerplate/workflow/env.yml`:

```bash
mamba env create -f sample_submission_boilerplate/workflow/env.yml
mamba activate crc-st-h5
```

Key dependencies:
- Python 3.11
- `scanpy`, `anndata`, `h5py`, `numpy`, `scipy`, `pandas`, `matplotlib`

> You do **not** need STAR, kallisto, Space Ranger, or other aligners/counters for this H5-only mini-benchmark.

---

## 6) Directory structure

```
sample_submission_boilerplate/
├── answers.yaml
├── data/
│   └── filtered_feature_bc_matrix.h5      # (you provide)
├── metadata.yaml
├── outs/
│   └── filtered_feature_bc_matrix.h5 -> data/filtered_feature_bc_matrix.h5
├── questions.yaml
├── results/
│   ├── filtered_feature_bc_matrix.h5 -> data/filtered_feature_bc_matrix.h5
│   ├── markers/
│   │   └── marker_scores.tsv
│   ├── metrics/
│   │   ├── spot_metrics.tsv
│   │   └── top_genes.txt
│   └── results.json
└── workflow/
    ├── 00_verify_h5.sh
    ├── 03_metrics.sh
    ├── 04_markers.sh
    ├── 05_evaluate.py
    ├── markers_crc_cms.tsv
    └── env.yml
```

---

## 7) Runtime expectations

- Data size: ~2 MB (single H5) in this mini-benchmark layout.
- Runtime: seconds to a few minutes on a laptop (depends on SciPy sparse operations).
- Human time to reproduce: ~10 minutes once the H5 is in place.

