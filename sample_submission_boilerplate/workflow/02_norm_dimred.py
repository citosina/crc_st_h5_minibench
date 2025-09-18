






import os, json, numpy as np, pandas as pd, scanpy as sc
from pathlib import Path
export CRCST_HVG_FLAVOR=seurat
root = Path(__file__).resolve().parents[1]
h5   = root/"results/filtered_feature_bc_matrix.h5"  # created by 00_verify_h5.sh
outd = root/"results/norm"; outd.mkdir(parents=True, exist_ok=True)

adata = sc.read_10x_h5(str(h5))
adata.var_names_make_unique()

# QC: write total feature count for Q1 (genes)
(Path(root/"results/qc").mkdir(parents=True, exist_ok=True))
(Path(root/"results/qc/features_total.txt")
 ).write_text(str(adata.n_vars) + "\n")

# Normalize + HVGs + PCA (Seurat v3 flavor)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat_v3")
adata = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)

# PCs to reach >=80% variance
vr = np.array(adata.uns["pca"]["variance_ratio"])
cum = vr.cumsum()
pcs_for_80 = int(np.searchsorted(cum, 0.80) + 1)

pd.DataFrame({"pc": np.arange(1, len(vr)+1), "variance_ratio": vr,
              "cumulative": cum}).to_csv(outd/"pca_variance.tsv", sep="\t", index=False)
(outd/"pcs_for_80pct.txt").write_text(str(pcs_for_80) + "\n")

adata.write_h5ad(outd/"adata_pca.h5ad")
print(f"[ok] norm/dimred complete; PCs≥80%: {pcs_for_80}")
