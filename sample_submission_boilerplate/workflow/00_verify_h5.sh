#!/usr/bin/env bash
set -euo pipefail

# project root (one level up from this script)
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# default H5 location; override with:
#   CRCST_H5=/abs/path/to/filtered_feature_bc_matrix.h5 bash workflow/00_verify_h5.sh
H5_DEFAULT="$ROOT/data/filtered_feature_bc_matrix.h5"
H5="${CRCST_H5:-$H5_DEFAULT}"

if [[ ! -f "$H5" ]]; then
  echo "[error] expected H5 at: $H5"
  echo "        Set CRCST_H5 to override, e.g.:"
  echo "        CRCST_H5=/absolute/path/to/filtered_feature_bc_matrix.h5 bash workflow/00_verify_h5.sh"
  exit 1
fi

# Resolve to absolute path (portable on macOS)
H5_ABS="$(cd "$(dirname "$H5")" && pwd)/$(basename "$H5")"

# Sanity: try opening the file and confirm the 'matrix' group exists
python - "$H5_ABS" <<'PY'
import sys, h5py
p = sys.argv[1]
with h5py.File(p, "r") as f:
    assert "matrix" in f.keys(), "H5 missing 'matrix' group"
print("[ok] H5 is readable:", p)
PY

# Primary outcome: link into results/
mkdir -p "$ROOT/results"
ln -sfn "$H5_ABS" "$ROOT/results/filtered_feature_bc_matrix.h5"
echo "[ok] linked results → $ROOT/results/filtered_feature_bc_matrix.h5"

# Optional compatibility: also link into outs/
mkdir -p "$ROOT/outs"
ln -sfn "$H5_ABS" "$ROOT/outs/filtered_feature_bc_matrix.h5"

# Show both links
ls -l "$ROOT/results/filtered_feature_bc_matrix.h5" || true
ls -l "$ROOT/outs/filtered_feature_bc_matrix.h5" || true
