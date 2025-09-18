import json, pandas as pd
from pathlib import Path

root = Path(__file__).resolve().parents[1]

# Q1: total genes (features)
q1 = int((root/"results/qc/features_total.txt").read_text().strip())

# Q2: PCs to reach >=80%
q2 = int((root/"results/norm/pcs_for_80pct.txt").read_text().strip())

# Q3: number of Leiden clusters at res=0.4
lab = pd.read_csv(root/"results/clustering/leiden_res0.4.tsv", sep="\t")
q3 = int(lab["cluster"].nunique())

# Q4: top 3 DE genes for EPCAM-high cluster
top3 = (root/"results/de/top3_marker_genes.txt").read_text().strip().splitlines()

# Q5: CMS call for EPCAM-high cluster
q5 = (root/"results/enrichment/epcam_cluster_cms.txt").read_text().strip()

out = {"q1": q1, "q2": q2, "q3": q3, "q4": top3, "q5": q5}
(root/"results").mkdir(exist_ok=True, parents=True)
(root/"results/results.json").write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
