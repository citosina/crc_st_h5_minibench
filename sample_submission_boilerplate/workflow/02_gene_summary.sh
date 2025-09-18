#!/usr/bin/env bash
set -euo pipefail

# Always work from the boilerplate root
cd "$(dirname "$0")/.."
ROOT="$(pwd)"
H5="${ROOT}/results/filtered_feature_bc_matrix.h5"

# Make sure the H5 is linked under results/
if [ ! -e "$H5" ]; then
  echo "[info] results H5 missing; running 00_verify_h5.sh to create it…"
  bash "$ROOT/workflow/00_verify_h5.sh"
fi

if [ ! -e "$H5" ]; then
  echo "[error] still missing: $H5"
  exit 1
fi

mkdir -p "$ROOT/results/qc"

# Write features_total.txt (number of genes/features in the matrix)
export ROOT
python - <<'PY'
import os, pathlib, scanpy as sc
root = pathlib.Path(os.environ["ROOT"])
h5   = root / "results" / "filtered_feature_bc_matrix.h5"
qc   = root / "results" / "qc"
adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()
n_features = int(adata.n_vars)
(qc / "features_total.txt").write_text(str(n_features) + "\n")
print(f"[ok] wrote {qc/'features_total.txt'} = {n_features}")
PY

echo "[done] gene summary"
