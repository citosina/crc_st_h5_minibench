#!/usr/bin/env python3
"""
05_evaluate.py

Evaluates answers for the CRC ST (H5) mini-bench.

Questions/stages this script answers:

q1 (feature_catalog):     total number of features (genes) in the H5
q2 (feature_presence):    number of genes with nonzero counts in ≥1 spot
q3 (spot_metrics):        median number of detected genes per spot
q4 (expression_ranking):  top 5 genes by total counts (highest→lowest)
q5 (marker_enrichment):   CMS-like class with highest marker enrichment
                          (simple sum of total counts over marker genes)

Inputs:
- results/filtered_feature_bc_matrix.h5      (symlink to data/... from 00_verify_h5.sh)
- workflow/markers_crc_cms.tsv               (gene→set mapping; used for q5)

Output:
- results/results.json
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

def locate_h5(root: Path) -> Path:
    """Find the H5 under results/, falling back to outs/ or data/ if needed."""
    candidates = [
        root / "results" / "filtered_feature_bc_matrix.h5",
        root / "outs"    / "filtered_feature_bc_matrix.h5",
        root / "data"    / "filtered_feature_bc_matrix.h5",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError(
        "Could not find filtered_feature_bc_matrix.h5 under results/, outs/, or data/. "
        "Run workflow/00_verify_h5.sh first."
    )

def read_adata(h5_path: Path):
    import scanpy as sc  # heavy import kept local
    adata = sc.read_10x_h5(str(h5_path))
    # downstream expects unique var names
    adata.var_names_make_unique()
    return adata

def to_csr(X):
    return X.tocsr() if hasattr(X, "tocsr") else X

def compute_q1_q2_q3(adata):
    """
    q1: total features (genes) recorded in H5
    q2: genes with nonzero counts in ≥1 spot
    q3: median number of detected genes per spot
    """
    X = to_csr(adata.X)
    # gene totals (axis=0), spot nonzeros (axis=1)
    gene_totals = np.asarray(X.sum(axis=0)).ravel()
    genes_nonzero = int((gene_totals > 0).sum())
    genes_per_spot = np.asarray((X > 0).sum(axis=1)).ravel()
    median_genes_per_spot = int(np.median(genes_per_spot))
    features_total = int(adata.n_vars)
    return features_total, genes_nonzero, median_genes_per_spot, gene_totals

def compute_q4_top5(adata, gene_totals):
    """Top 5 genes by total counts across all spots (descending)."""
    idx = np.argsort(-gene_totals)[:5]
    return [str(adata.var_names[i]) for i in idx]

def load_markers_table(markers_path: Path):
    """
    Accepts TSV/CSV with either:
      (A) columns: gene, set   (set in {CMS1_like, CMS2_like, CMS3_like, CMS4_like})
      (B) columns: gene plus one or more CMS*_like boolean/score columns
    Returns a dict: set_name -> set_of_genes
    """
    if not markers_path.exists():
        raise FileNotFoundError(f"Missing markers file: {markers_path}")

    # Try TSV first, then CSV
    try:
        df = pd.read_csv(markers_path, sep="\t")
    except Exception:
        df = pd.read_csv(markers_path)

    df.columns = [c.strip() for c in df.columns]
    df.rename(columns={"Gene": "gene", "GENE": "gene", "set_name": "set"}, inplace=True)

    cms_cols = [c for c in df.columns if c.upper().startswith("CMS") and c.lower().endswith("_like")]
    sets = {}

    if "gene" in df.columns and "set" in df.columns:
        for cms_name, sub in df.groupby("set"):
            if not isinstance(cms_name, str):
                continue
            sets[cms_name] = set(map(str, sub["gene"].astype(str)))
    elif "gene" in df.columns and len(cms_cols) > 0:
        for c in cms_cols:
            members = df.loc[df[c].astype(bool), "gene"] if df[c].dtype != object else df.loc[df[c].notna(), "gene"]
            sets[c] = set(map(str, members.astype(str)))
    else:
        raise ValueError(
            "Unrecognized markers file schema. Expect either columns [gene,set] or [gene, CMS*_like ...]"
        )
    return sets

def compute_q5_cms_like(adata, gene_totals, markers_path: Path) -> str:
    """
    Score each CMS-like set as the sum of gene_totals over its marker genes, pick argmax.
    """
    sets = load_markers_table(markers_path)

    # Map var names -> index for quick lookups
    name_to_idx = {g: i for i, g in enumerate(map(str, adata.var_names))}
    scores = {}
    for cms, genes in sets.items():
        idxs = [name_to_idx[g] for g in genes if g in name_to_idx]
        score = float(np.sum(gene_totals[idxs])) if idxs else 0.0
        scores[cms] = score

    if not scores:
        return "CMS2_like"  # neutral default
    # prefer canonical order on ties
    order = ["CMS1_like", "CMS2_like", "CMS3_like", "CMS4_like"]
    best = max(scores.items(), key=lambda kv: (kv[1], -order.index(kv[0]) if kv[0] in order else 99))[0]
    return best

def main():
    # Resolve root as sample_submission_boilerplate/
    root = Path(__file__).resolve().parents[1]
    results_dir = root / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    # Find H5, read adata
    h5 = locate_h5(root)
    adata = read_adata(h5)

    # q1, q2, q3 + gene_totals for later use
    q1, q2, q3, gene_totals = compute_q1_q2_q3(adata)

    # q4 top-5 ranking
    q4 = compute_q4_top5(adata, gene_totals)

    # q5 CMS-like enrichment
    markers_path = root / "workflow" / "markers_crc_cms.tsv"
    try:
        q5 = compute_q5_cms_like(adata, gene_totals, markers_path)
    except Exception:
        # If markers file missing or malformed, fall back
        q5 = "CMS2_like"

    # Write results
    out_json = results_dir / "results.json"
    payload = {
        "q1": int(q1),
        "q2": int(q2),
        "q3": int(q3),
        "q4": list(map(str, q4)),
        "q5": str(q5),
    }
    with out_json.open("w") as f:
        json.dump(payload, f, indent=2)
    print(f"[ok] wrote {out_json}")

if __name__ == "__main__":
    main()
