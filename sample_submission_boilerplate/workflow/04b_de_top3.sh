#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
H5="${ROOT}/results/filtered_feature_bc_matrix.h5"
CL="${ROOT}/results/clustering/leiden_res0.4.tsv"
OUTD="${ROOT}/results/de"
mkdir -p "${OUTD}"

[[ -e "${H5}" ]] || { echo "ERROR: Missing ${H5} (run workflow/00_verify_h5.sh first)"; exit 1; }
[[ -e "${CL}" ]] || { echo "ERROR: Missing ${CL} (run workflow/04a_cluster_leiden.sh first)"; exit 1; }

python - "${H5}" "${CL}" "${OUTD}" <<'PY'
import sys, numpy as np, pandas as pd, scanpy as sc
from pathlib import Path

h5, lab, outd = map(Path, sys.argv[1:4])

adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()

# attach cluster labels
labdf = pd.read_csv(lab, sep="\t")  # columns: barcode, cluster
labdf.index = labdf["barcode"].astype(str)

# align
keep = adata.obs_names.intersection(labdf.index)
adata = adata[keep].copy()
adata.obs["cluster"] = labdf.loc[adata.obs_names, "cluster"].astype(str).values

# pick top-expressed cluster by median UMIs
X = adata.X.tocsr() if hasattr(adata.X, "tocsr") else adata.X
umis = np.asarray(X.sum(axis=1)).ravel()
adata.obs["umis"] = umis
topc = adata.obs.groupby("cluster")["umis"].median().sort_values(ascending=False).index[0]

# normalize & log
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# rank genes: top cluster vs rest (NOTE: must pass groupby)
sc.tl.rank_genes_groups(adata, groupby="cluster", groups=[topc], reference="rest", method="wilcoxon")
dfres = sc.get.rank_genes_groups_df(adata, group=topc).sort_values("scores", ascending=False)
top3 = dfres["names"].head(3).tolist()

outp = Path(outd) / "top3_marker_genes.txt"
outp.write_text("\n".join(top3) + "\n")
print(f"[ok] wrote {outp} for cluster={topc} → {top3}")
PY
