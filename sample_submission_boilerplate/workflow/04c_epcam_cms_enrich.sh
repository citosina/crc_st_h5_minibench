#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
H5="${ROOT}/results/filtered_feature_bc_matrix.h5"
CL="${ROOT}/results/clustering/leiden_res0.4.tsv"
MARK="${ROOT}/workflow/markers_crc_cms.tsv"
OUTD="${ROOT}/results/enrichment"
mkdir -p "${OUTD}"

[[ -e "${H5}"  ]] || { echo "ERROR: Missing ${H5} (run workflow/00_verify_h5.sh first)"; exit 1; }
[[ -e "${CL}"  ]] || { echo "ERROR: Missing ${CL} (run workflow/04a_cluster_leiden.sh first)"; exit 1; }
[[ -e "${MARK}" ]] || { echo "ERROR: Missing ${MARK}"; exit 1; }

python - "${H5}" "${CL}" "${MARK}" "${OUTD}" <<'PY'
import sys, numpy as np, pandas as pd, scanpy as sc
from pathlib import Path

h5, labp, markp, outd = map(Path, sys.argv[1:5])

adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()

# attach Leiden labels
lab = pd.read_csv(labp, sep="\t")  # columns: barcode, cluster
lab.index = lab["barcode"].astype(str)
keep = adata.obs_names.intersection(lab.index)
adata = adata[keep].copy()
adata.obs["cluster"] = lab.loc[adata.obs_names, "cluster"].astype(str).values

# normalize for scoring
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# --- 1) choose EPCAM-high cluster ---
def mean_gene_by_cluster(adata, gene):
    g = gene if gene in adata.var_names else None
    if g is None:
        return None
    vals = adata[:, g].X
    if hasattr(vals, "A"): vals = vals.A
    vals = np.asarray(vals).ravel()
    return adata.obs.assign(v=vals).groupby("cluster")["v"].mean().sort_values(ascending=False)

# primary marker
rank = mean_gene_by_cluster(adata, "EPCAM")
if rank is None:  # fallback epithelial markers
    epi_fallback = [g for g in ["KRT8","KRT18","KRT19","KRT7"] if g in adata.var_names]
    if not epi_fallback:
        # final fallback: pick cluster with highest median total UMIs
        X = adata.X.tocsr() if hasattr(adata.X,"tocsr") else adata.X
        umis = np.asarray(X.sum(axis=1)).ravel()
        adata.obs["umis"] = umis
        top_cluster = adata.obs.groupby("cluster")["umis"].median().idxmax()
    else:
        sub = adata[:, epi_fallback].X
        if hasattr(sub, "A"): sub = sub.A
        sub = np.asarray(sub)
        score = sub.mean(axis=1)
        adata.obs["epi_score"] = score
        top_cluster = adata.obs.groupby("cluster")["epi_score"].mean().idxmax()
else:
    top_cluster = rank.index[0]

# --- 2) CMS set scoring on the EPCAM-high cluster ---
mk = pd.read_csv(markp, sep="\t")
# try to auto-detect column names
def pick(colnames, candidates):
    for c in colnames:
        if c.lower() in candidates: return c
    return None
gcol = pick(mk.columns, {"gene","symbol","gene_name","gene_symbol"})
scol = pick(mk.columns, {"cms","set","class","label","subtype"})
if gcol is None or scol is None:
    raise SystemExit(f"Cannot infer columns in {markp}. Need gene + cms/set columns.")
mk[gcol] = mk[gcol].astype(str)
mk[scol] = mk[scol].astype(str)

# restrict to genes present
mk = mk[ mk[gcol].isin(adata.var_names) ].copy()
cms_sets = {s: mk.loc[mk[scol]==s, gcol].tolist() for s in sorted(mk[scol].unique()) if len(mk.loc[mk[scol]==s])>0}

# cluster mask
cl_mask = (adata.obs["cluster"] == top_cluster).values

def avg_expr(adata, genes, mask):
    if not genes: return -np.inf
    sub = adata[:, genes].X
    if hasattr(sub, "A"): sub = sub.A
    sub = np.asarray(sub)
    return float(sub[mask].mean())

scores = []
for s, genes in cms_sets.items():
    in_mean  = avg_expr(adata, genes, cl_mask)
    out_mean = avg_expr(adata, genes, ~cl_mask)
    diff = in_mean - out_mean
    scores.append((s, in_mean, out_mean, diff))

scoredf = pd.DataFrame(scores, columns=["cms","mean_in","mean_out","diff"]).sort_values("diff", ascending=False)
best = scoredf.iloc[0]["cms"]

# coerce label to the expected exact string (e.g., CMS4_like)
label = best if best.endswith("_like") else (best + "_like" if best.upper().startswith("CMS") else f"{best}")

outd = Path(outd)
outd.mkdir(parents=True, exist_ok=True)
(outd/"epcam_cluster.txt").write_text(str(top_cluster) + "\n")
(outd/"epcam_scores.tsv").write_text(scoredf.to_csv(sep="\t", index=False))
(outd/"epcam_cluster_cms.txt").write_text(label + "\n")

print(f"[ok] EPCAM-high cluster={top_cluster}; best CMS={label}")
PY
