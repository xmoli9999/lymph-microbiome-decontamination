#!/usr/bin/env python3
"""
Step 4: CLR-based compositional analysis (approximating ANCOM-BC2).

While ANCOM-BC2 provides a rigorous compositional framework, its performance
is unstable in small-sample settings (n=5 per group). We therefore use a
CLR-based approximation with sample-specific bias correction.

CLR transformation with pseudocount 0.5, bias correction, Welch's t-test.
"""

import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
import os

OUTPUT_DIR = "../results/04_clr/"
PSEUDOCOUNT = 0.5

def clr_transform(counts, pseudocount=PSEUDOCOUNT):
    """Centered log-ratio transformation."""
    log_counts = np.log(counts + pseudocount)
    geometric_mean = log_counts.mean(axis=0)
    clr = log_counts - geometric_mean
    return clr

def estimate_bias(clr_data):
    """Estimate sample-specific bias as mean log-count deviation from grand mean."""
    grand_mean = clr_data.mean().mean()
    sample_bias = clr_data.mean(axis=0) - grand_mean
    return sample_bias

def bias_corrected_clr(counts, pseudocount=PSEUDOCOUNT):
    """Apply CLR transformation with bias correction."""
    clr = clr_transform(counts, pseudocount)
    bias = estimate_bias(clr)
    corrected = clr.subtract(bias, axis=1)
    return corrected

def pairwise_test(clr_data, meta, ref_time="0h", test_time="4h"):
    """Welch's t-test on bias-corrected CLR values for each taxon."""
    ref_samples = meta[meta["time_point"] == ref_time]["sample_id"].tolist()
    test_samples = meta[meta["time_point"] == test_time]["sample_id"].tolist()

    ref_samples = [s for s in ref_samples if s in clr_data.columns]
    test_samples = [s for s in test_samples if s in clr_data.columns]

    results = []
    for taxon in clr_data.index:
        ref_vals = clr_data.loc[taxon, ref_samples].values.astype(float)
        test_vals = clr_data.loc[taxon, test_samples].values.astype(float)

        stat, p = ttest_ind(test_vals, ref_vals, equal_var=False)
        clr_diff = test_vals.mean() - ref_vals.mean()

        results.append({
            "taxon": taxon,
            "comparison": f"{ref_time}_vs_{test_time}",
            "clr_difference": clr_diff,
            "t_statistic": stat,
            "p_value": p,
            "direction": "up" if clr_diff > 0 else "down"
        })

    df = pd.DataFrame(results)

    # FDR correction
    reject, padj, _, _ = multipletests(df["p_value"], method="fdr_bh")
    df["padj"] = padj
    df["fdr_significant"] = reject

    return df

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load decontaminated + aggregated data
    otu = pd.read_csv("../results/01_decontamination/otu_table_decontaminated.tsv", sep="\t", index_col=0)
    meta = pd.read_csv("../data/sample_metadata.csv")
    taxonomy = pd.read_csv("../data/taxonomy.tsv", sep="\t", index_col=0)
    bio_meta = meta[meta["sample_type"] != "NegControl"]

    # Aggregate to family level
    tax_map = taxonomy["family"].fillna("Unclassified")
    otu_fam = otu.copy()
    otu_fam["taxon"] = otu_fam.index.map(lambda x: tax_map.get(x, "Unclassified"))
    family_counts = otu_fam.groupby("taxon").sum()

    # CLR transform with bias correction
    clr_data = bias_corrected_clr(family_counts)

    # Run all pairwise comparisons
    time_points = ["4h", "8h", "12h", "24h"]
    all_results = []

    for tp in time_points:
        res = pairwise_test(clr_data, bio_meta, ref_time="0h", test_time=tp)
        all_results.append(res)
        sig = res[res["fdr_significant"]]
        nom_sig = res[res["p_value"] < 0.05]
        print(f"0h vs {tp}: {len(sig)} FDR-significant, {len(nom_sig)} nominal P<0.05")

    combined = pd.concat(all_results, ignore_index=True)
    combined.to_csv(os.path.join(OUTPUT_DIR, "clr_analysis_results.csv"), index=False)
    print(f"\nResults saved to {OUTPUT_DIR}")
