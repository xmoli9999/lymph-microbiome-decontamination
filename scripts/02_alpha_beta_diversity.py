#!/usr/bin/env python3
"""
Step 2: Alpha and beta diversity analysis (before and after decontamination).

Alpha: Shannon, Simpson, Chao1 (rarefied to 5000 reads)
Beta: Bray-Curtis dissimilarity, PCoA, PERMANOVA (999 permutations)
"""

import pandas as pd
import numpy as np
from scipy.stats import kruskal
from skbio.diversity import alpha_diversity, beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
import os

OUTPUT_DIR = "../results/02_diversity/"
RAREFACTION_DEPTH = 5000

def rarefy(counts, depth):
    """Rarefy a single sample to a given depth."""
    if counts.sum() < depth:
        return None
    probs = counts / counts.sum()
    return np.bincount(
        np.random.choice(len(counts), size=depth, p=probs),
        minlength=len(counts)
    )

def compute_alpha_diversity(otu, meta, depth=RAREFACTION_DEPTH):
    """Compute alpha diversity indices on rarefied data."""
    results = []
    for _, row in meta.iterrows():
        sid = row["sample_id"]
        if sid not in otu.columns:
            continue
        counts = otu[sid].values
        rarefied = rarefy(counts, depth)
        if rarefied is None:
            print(f"  Warning: {sid} has < {depth} reads, excluded from rarefied analysis")
            continue
        results.append({
            "sample_id": sid,
            "group": row["group"],
            "time_point": row.get("time_point", ""),
            "shannon": alpha_diversity("shannon", rarefied),
            "simpson": alpha_diversity("simpson", rarefied),
            "chao1": alpha_diversity("chao1", rarefied),
        })
    return pd.DataFrame(results)

def compute_beta_diversity(otu, meta):
    """Compute Bray-Curtis dissimilarity and PERMANOVA."""
    samples = [s for s in meta["sample_id"] if s in otu.columns]
    data = otu[samples].T.values
    ids = samples

    dm = beta_diversity("braycurtis", data, ids=ids)
    pc = pcoa(dm)

    # PERMANOVA
    groups = meta.set_index("sample_id").loc[samples, "time_point"]
    result = permanova(dm, groups, permutations=999)

    return dm, pc, result

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load decontaminated OTU table
    otu = pd.read_csv("../results/01_decontamination/otu_table_decontaminated.tsv", sep="\t", index_col=0)
    meta = pd.read_csv("../data/sample_metadata.csv")
    bio_meta = meta[meta["sample_type"] != "NegControl"]

    # Alpha diversity
    alpha_df = compute_alpha_diversity(otu, bio_meta)
    alpha_df.to_csv(os.path.join(OUTPUT_DIR, "alpha_diversity.csv"), index=False)

    # Kruskal-Wallis test per metric
    for metric in ["shannon", "simpson", "chao1"]:
        groups = [g[metric].values for _, g in alpha_df.groupby("time_point")]
        stat, p = kruskal(*groups)
        print(f"{metric}: KW H={stat:.3f}, P={p:.4f}")

    # Beta diversity
    dm, pc, permanova_result = compute_beta_diversity(otu, bio_meta)
    pc.samples.to_csv(os.path.join(OUTPUT_DIR, "pcoa_coordinates.csv"))
    print(f"\nPERMANOVA: R²={permanova_result['test statistic']:.3f}, P={permanova_result['p-value']:.4f}")
