#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
H5="${ROOT}/results/filtered_feature_bc_matrix.h5"
OUTD="${ROOT}/results/clustering"
mkdir -p "$OUTD"

if [[ ! -e "$H5" ]]; then
  echo "ERROR: Missing $H5 (run workflow/00_verify_h5.sh first)"; exit 1
fi

# keep consistent with normalization step (avoids scikit-misc if seurat/cell_ranger)
export CRCST_HVG_FLAVOR="${CRCST_HVG_FLAVOR:-seurat}"

python - "$H5" "$OUTD" <<'PY'
import os, sys, numpy as np, pandas as pd, scanpy as sc
from pathlib import Path

h5   = Path(sys.argv[1])            # passed from bash; absolute or repo-local path
outd = Path(sys.argv[2]); outd.mkdir(parents=True, exist_ok=True)

flavor = os.environ.get("CRCST_HVG_FLAVOR","seurat")
if flavor not in {"seurat","cell_ranger","seurat_v3"}:
    flavor = "seurat"

# reuse PCs-for-80% variance if present (cap 10..50)
root = outd.parents[1]
pcs_file = root/"results/norm/pcs_for_80pct.txt"
n_pcs = 30
if pcs_file.exists():
    try:
        n_pcs = int(pcs_file.read_text().strip())
        n_pcs = max(10, min(n_pcs, 50))
    except Exception:
        n_pcs = 30

adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=2000, subset=True)
sc.pp.scale(adata, zero_center=True, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=max(n_pcs, 30))
sc.pp.neighbors(adata, n_pcs=n_pcs)
sc.tl.leiden(adata, resolution=0.4, key_added="leiden_r0_4")

pd.DataFrame({"barcode": adata.obs_names,
              "cluster": adata.obs["leiden_r0_4"].astype(str)}) \
  .to_csv(outd/"leiden_res0.4.tsv", sep="\t", index=False)

print(f"[ok] wrote {outd/'leiden_res0.4.tsv'} using n_pcs={n_pcs}, flavor={flavor}")
PY
