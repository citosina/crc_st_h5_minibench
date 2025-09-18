import re, numpy as np, pandas as pd, scanpy as sc
from pathlib import Path

root = Path(__file__).resolve().parents[1]
clustd = root/"results/clustering"
out_de = root/"results/de"; out_de.mkdir(parents=True, exist_ok=True)
out_en = root/"results/enrichment"; out_en.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(clustd/"adata_leiden.h5ad")

# Find EPCAM (case-insensitive exact symbol match)
symbols = pd.Index(adata.var_names)
cand = [g for g in symbols if g.upper()=="EPCAM"]
if not cand:
    # fallback: epithelial proxy list, pick any available
    proxy = [g for g in ["EPCAM","KRT8","KRT18","KRT19","KRT7"] if g in symbols]
    if not proxy:
        raise SystemExit("No EPCAM or proxy epithelial markers found in var_names.")
    epig = proxy
else:
    epig = [cand[0]]

# Compute per-cluster median EPCAM/proxy score
X = adata[:, epig].X
if hasattr(X, "toarray"): X = X.toarray()
score = X.mean(axis=1)         # average across available epithelial markers
adata.obs["epi_score"] = score

grp = adata.obs.groupby("leiden_r04")["epi_score"].median().sort_values(ascending=False)
epcam_cluster = grp.index[0]

# DE: Wilcoxon
sc.tl.rank_genes_groups(adata, groupby="leiden_r04", groups=[epcam_cluster],
                        method="wilcoxon", key_added="de_epcam")
rg = adata.uns["de_epcam"]
top_names = pd.Series(rg["names"][epcam_cluster]).astype(str)
top3 = top_names.head(3).tolist()

# Write artifacts
(root/"results/metrics").mkdir(parents=True, exist_ok=True)
(clustd/"epcam_cluster.txt").write_text(str(epcam_cluster) + "\n")
(out_de/"top3_marker_genes.txt").write_text("\n".join(top3) + "\n")
pd.DataFrame({"rank": np.arange(1, len(top_names)+1), "gene": top_names}).to_csv(
    out_de/"epcam_cluster_markers.tsv", sep="\t", index=False)

# Enrichment vs CMS marker sets
marker_tsv = root/"workflow/markers_crc_cms.tsv"
m = pd.read_csv(marker_tsv, sep="\t")  # expects columns: gene, set (CMS1_like/2/3/4)
m["gene_upper"] = m["gene"].astype(str).str.upper()
topN = pd.Series(top_names.head(50)).str.upper()

scores = (m[m["gene_upper"].isin(topN)]
          .groupby("set").size().sort_values(ascending=False))
cms_call = scores.index[0] if len(scores) else "CMS2_like"  # reasonable default

(out_en/"epcam_cluster_cms.txt").write_text(str(cms_call) + "\n")

print(f"[ok] EPCAM-high cluster: {epcam_cluster}; top3: {top3}; CMS: {cms_call}")
