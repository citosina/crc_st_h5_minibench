#!/usr/bin/env bash
set -euo pipefail

# Paths
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
H5="${ROOT}/results/filtered_feature_bc_matrix.h5"
OUTD="${ROOT}/results/norm"

# Checks
if [[ ! -e "$H5" ]]; then
  echo "ERROR: Missing $H5 (run workflow/00_verify_h5.sh first)"; exit 1
fi
mkdir -p "$OUTD"

echo "[info] PCA on normalized data → $OUTD"

# Use CRCST_HVG_FLAVOR=seurat (default) unless overridden (valid: seurat, cell_ranger, seurat_v3)
H5="$H5" OUTD="$OUTD" CRCST_HVG_FLAVOR="${CRCST_HVG_FLAVOR:-seurat}" python - <<'PY'
import os, numpy as np, pandas as pd, scanpy as sc
from pathlib import Path

h5   = Path(os.environ["H5"])
outd = Path(os.environ["OUTD"])
flavor = os.environ.get("CRCST_HVG_FLAVOR","seurat")
if flavor not in {"seurat","cell_ranger","seurat_v3"}:
    flavor = "seurat"

adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()

# Normalize & log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# HVGs (no skmisc dependency with 'seurat'/'cell_ranger')
sc.pp.highly_variable_genes(adata, flavor=flavor, n_top_genes=2000, subset=True)

# PCA
sc.pp.scale(adata, zero_center=True, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

expl = adata.uns["pca"]["variance_ratio"]
cum  = np.cumsum(expl)
k = int(np.argmax(cum >= 0.80) + 1) if np.any(cum >= 0.80) else len(cum)

(outd/"pcs_for_80pct.txt").write_text(f"{k}\n")
pd.DataFrame({"pc": np.arange(1, len(expl)+1),
              "var_ratio": expl,
              "cum": cum}).to_csv(outd/"pca_explained_variance.tsv", sep="\t", index=False)

print(f"[ok] wrote {outd/'pcs_for_80pct.txt'} = {k}")
PY

echo "[done] PCA summary"
