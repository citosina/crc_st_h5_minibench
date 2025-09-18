import numpy as np, pandas as pd, scanpy as sc
from pathlib import Path

root = Path(__file__).resolve().parents[1]
normd = root/"results/norm"
outd  = root/"results/clustering"; outd.mkdir(parents=True, exist_ok=True)

adata = sc.read_h5ad(normd/"adata_pca.h5ad")
pcs   = int((normd/"pcs_for_80pct.txt").read_text().strip() or "20")

sc.pp.neighbors(adata, n_pcs=pcs)
sc.tl.leiden(adata, resolution=0.4, key_added="leiden_r04")

# Save labels
lab = pd.DataFrame({"barcode": adata.obs_names, "cluster": adata.obs["leiden_r04"].astype(str)})
lab.to_csv(outd/"leiden_res0.4.tsv", sep="\t", index=False)

adata.write_h5ad(outd/"adata_leiden.h5ad")
print(f"[ok] clustering complete; clusters: {lab['cluster'].nunique()}")
