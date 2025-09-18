







#!/usr/bin/env bash
set -euo pipefail

# Force HVG flavor to Seurat (avoids scikit-misc requirement on macOS ARM)
export CRCST_HVG_FLAVOR=seurat

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

H5="${CRCST_H5:-$ROOT/results/filtered_feature_bc_matrix.h5}"
[[ -f "$H5" ]] || { echo "[error] missing H5: $H5"; exit 1; }

METDIR="$ROOT/results/metrics"
mkdir -p "$METDIR"
OUT_TSV="$METDIR/spot_metrics.tsv"

python - "$H5" "$OUT_TSV" <<'PY'
import scanpy as sc, pandas as pd, numpy as np, sys
from pathlib import Path

h5  = Path(sys.argv[1])
out = Path(sys.argv[2])

adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()

X = adata.X.tocsr() if hasattr(adata.X, "tocsr") else adata.X
umis  = np.asarray(X.sum(axis=1)).ravel()
genes = np.asarray((X>0).sum(axis=1)).ravel()

metrics = {
    "spots_total":            int(adata.n_obs),
    "spots_nonzero":          int((umis > 0).sum()),
    "median_umis":            float(np.median(umis)),
    "median_genes":           float(np.median(genes)),
    "median_umis_per_spot":   float(np.median(umis)),
    "median_genes_per_spot":  float(np.median(genes)),
    "mean_umis":              float(np.mean(umis)),
    "mean_genes":             float(np.mean(genes)),
    "mean_umis_per_spot":     float(np.mean(umis)),
    "mean_genes_per_spot":    float(np.mean(genes)),
}

pd.DataFrame({"metric": list(metrics.keys()), "value": list(metrics.values())}) \
  .to_csv(out, sep="\t", index=False)
print(f"[ok] wrote {out}")
PY
