#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.."; pwd)"
H5="${ROOT}/results/filtered_feature_bc_matrix.h5"

mkdir -p "${ROOT}/results/metrics"

python - <<'PY'
import scanpy as sc, numpy as np, pandas as pd, sys, pathlib
root = pathlib.Path(__file__).resolve().parents[1]
h5 = root / "results" / "filtered_feature_bc_matrix.h5"
assert h5.exists(), f"Missing {h5}"

adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()

X = adata.X.tocsr() if hasattr(adata.X, "tocsr") else adata.X
genes_nonzero = int((np.asarray(X.sum(axis=0)).ravel() > 0).sum())  # any expression across spots
features_total = int(adata.n_vars)

df = pd.DataFrame([
    {"metric": "features_total", "value": features_total},
    {"metric": "genes_nonzero",  "value": genes_nonzero},
])
(df).to_csv(root / "results" / "metrics" / "gene_summary.tsv", sep="\t", index=False)
print("[ok] wrote results/metrics/gene_summary.tsv")
PY
