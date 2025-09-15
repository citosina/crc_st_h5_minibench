#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

H5="${CRCST_H5:-$ROOT/results/filtered_feature_bc_matrix.h5}"
[[ -f "$H5" ]] || { echo "[error] missing H5: $H5"; exit 1; }

MKDIR="$ROOT/results/markers"
METDIR="$ROOT/results/metrics"
mkdir -p "$MKDIR" "$METDIR"

MK_TSV="$MKDIR/marker_scores.tsv"
TOP5_TXT="$METDIR/top_genes.txt"

python - "$H5" "$MK_TSV" "$TOP5_TXT" <<'PY'
import scanpy as sc, pandas as pd, numpy as np, sys
from pathlib import Path

h5       = Path(sys.argv[1])
mk_tsv   = Path(sys.argv[2])
top5_txt = Path(sys.argv[3])

adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()

X = adata.X.tocsr() if hasattr(adata.X, "tocsr") else adata.X
gene_totals = np.asarray(X.sum(axis=0)).ravel()

mk = (pd.DataFrame({"gene": adata.var_names, "score": gene_totals})
        .sort_values("score", ascending=False))

mk.to_csv(mk_tsv, sep="\t", index=False)
top5 = mk["gene"].head(5).tolist()
top5_txt.write_text("\n".join(top5) + "\n")

print(f"[ok] wrote {mk_tsv}")
print(f"[ok] wrote {top5_txt} (top 5 genes)")
PY
