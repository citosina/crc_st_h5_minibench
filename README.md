# crc_st_h5_minibench

This mini-benchmark recreates a **condensed Visium spatial transcriptomics analysis**
starting **from a processed 10x h5 matrix** (`filtered_feature_bc_matrix.h5`) and asks a
set of questions that span distinct stages of the analysis (QC, PCA,
clustering, DE and epithelial CMS enrichment).

> Paper data source: Visium CRC dataset; processed Space Ranger outputs are
> available on Zenodo (DOI: **10.5281/zenodo.7551712**). We use the Rep1
> `filtered_feature_bc_matrix.h5` as the starting point for both problems.

---

## 1) Data sources

- **Primary matrix**: `filtered_feature_bc_matrix.h5` (10x Genomics format), obtained from the
  Zenodo record above (Rep1).  
- This repository **does not** commit raw FASTQs or Space Ranger outputs.
- Two zips are typically prepared outside this repo:
  - **crc_st_h5_minibench_data_only.zip** – contains only the H5 under
    `sample_submission_boilerplate/data/`.
  - **crc_st_h5_minibench_with_workflow.zip** – includes the H5 **plus**
    the workflow scripts in `sample_submission_boilerplate/workflow/` (no
    derived results).

Directory layout expected **after** unzipping/copying the H5:

```
sample_submission_boilerplate/
└── data/
    └── filtered_feature_bc_matrix.h5
```

---

## 2) How to get the data into place

If you have the prepared “data-only” zip:

```bash
# from repo root
unzip -o crc_st_h5_minibench_data_only.zip -d ./
# you should now have: sample_submission_boilerplate/data/filtered_feature_bc_matrix.h5
```

Or copy the H5 you downloaded from Zenodo into the path below:

```
sample_submission_boilerplate/data/filtered_feature_bc_matrix.h5
```

---

## 3) Environment

A small, reproducible conda/mamba environment is defined at
`sample_submission_boilerplate/workflow/env.yml`.

```bash
# create & activate
mamba env create -f sample_submission_boilerplate/workflow/env.yml -n crc-st-h5
mamba activate crc-st-h5

# if needed (depending on platform), install Leiden:
mamba install -c conda-forge leidenalg
```

> Tested with Python 3.11 on macOS (arm64). Key deps: `scanpy`, `anndata`,
> `numpy`, `pandas`, `h5py`, `leidenalg`.

---

## 4) Workflow steps

All scripts live in `sample_submission_boilerplate/workflow/`.
They assume the H5 is in `sample_submission_boilerplate/data/`.
Outputs are written under `sample_submission_boilerplate/results/`.

### Step 0 – Verify H5 and prepare links
Verifies the H5 is readable and symlinks it into `results/` (and `outs/` for compatibility).

```bash
bash workflow/00_verify_h5.sh
```

### Step 1 – Gene summary (feature count)
Writes total number of features/genes detected in the matrix.

```bash
bash workflow/02_gene_summary.sh
# outputs: results/qc/features_total.txt
```

### Step 2 – Normalize & PCA
Normalizes counts, selects HVGs, runs PCA and writes how many PCs explain ≥80% variance.

```bash
bash workflow/02b_normalize_pca.sh
# outputs: results/norm/pcs_for_80pct.txt
# env var: export CRCST_HVG_FLAVOR=seurat    # (default) avoids scikit-misc
```

### Step 3 – Clustering (Leiden)
Builds neighbors (using the PCs above) and runs Leiden clustering at resolution 0.4.

```bash
bash workflow/04a_cluster_leiden.sh
# outputs: results/clustering/leiden_res0.4.tsv
```

### Step 4 – Differential expression (top 3 genes for one cluster)
Ranks genes for the epithelial-enriched cluster and writes the **top 3** marker genes.

```bash
bash workflow/04b_de_top3.sh
# outputs: results/de/top3_marker_genes.txt
```

### Step 5 – Basic metrics & expression ranking
Computes per-spot QC and global expression ranking; writes top 5 genes by total counts.

```bash
bash workflow/03_metrics.sh
bash workflow/04_markers.sh
# outputs: results/metrics/spot_metrics.tsv, results/metrics/top_genes.txt,
#          results/markers/marker_scores.tsv
```

### Step 6 – Epithelial CMS enrichment
Scores epithelial (EPCAM-high) cluster against CMS-like epithelial signatures.

```bash
bash workflow/04c_epcam_cms_enrich.sh
# outputs: results/enrichment/epcam_cluster_cms.txt
```

### Step 7 – Evaluation
Collects the requested answers into `results/results.json`.

```bash
python workflow/05_evaluate.py
# outputs: results/results.json
```

---

## 5) Questions & answers

- **Questions**: `sample_submission_boilerplate/questions.yaml`  
  The five questions span **different stages** of the analysis:
  - Q1: QC/feature counting
  - Q2: Normalization/PCA (PCs to explain 80% variance)
  - Q3: Clustering (number of clusters at res=0.4)
  - Q4: Differential expression (top 3 genes for epithelial-enriched cluster)
  - Q5: Epithelial CMS enrichment (best-matching CMS-like signature)
- **Answers**: `sample_submission_boilerplate/answers.yaml` (with tolerances), which
  match the outputs produced by the scripts above.

---

## 6) Directory structure

```
sample_submission_boilerplate/
├── answers.yaml
├── data/
│   └── filtered_feature_bc_matrix.h5
├── metadata.yaml
├── outs/                              # symlinked convenience path to the H5
├── questions.yaml
├── results/                           # generated outputs (not committed)
│   ├── qc/
│   ├── norm/
│   ├── clustering/
│   ├── de/
│   ├── enrichment/
│   ├── markers/
│   └── metrics/
└── workflow/
    ├── 00_verify_h5.sh
    ├── 02_gene_summary.sh
    ├── 02b_normalize_pca.sh
    ├── 03_metrics.sh
    ├── 04a_cluster_leiden.sh
    ├── 04b_de_top3.sh
    ├── 04_markers.sh
    ├── 04c_epcam_cms_enrich.sh
    ├── 05_evaluate.py
    ├── env.yml
    └── markers_crc_cms.tsv
```

---

## 7) Runtime expectations

- **Data size**: the H5 is ~2 MB (Rep1 filtered matrix).
- **Runtime**: typically **< 5–10 minutes** on a laptop (8 cores) for all steps.
- **Reproducibility**: setting the same Scanpy seeds and using the provided env.yml
  should yield identical or near-identical outputs within the specified tolerances.

---

## 8) Notes for “no-workflow” problem

The **no-workflow** variant uses the same input H5 but expects answers to be derived
from the data with standard bioinformatics tooling. The questions remain stage-based,
but can be answered using popular libraries (e.g., Scanpy/Seurat) without relying on
the provided shell scripts.

---

## 9) Acknowledgments & citation

If you use this mini-benchmark, please cite the originating CRC spatial transcriptomics
data (Zenodo DOI **10.5281/zenodo.7551712**). See `metadata.yaml` for a full paper
reference and tool versions.
