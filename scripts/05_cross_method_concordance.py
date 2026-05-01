#!/usr/bin/env python3
"""
Step 5: Cross-method concordance assessment.

A taxon-comparison pair is considered concordant if:
  - Nominal P < 0.05 in both DESeq2 and CLR-based methods
  - Same directional change (both up or both down)
"""

import pandas as pd
import os

OUTPUT_DIR = "../results/05_concordance/"

def assess_concordance(deseq2_results, clr_results, p_threshold=0.05):
    """Find concordant taxon-comparison pairs across methods."""
    # Standardize column names
    deseq2 = deseq2_results.rename(columns={"padj": "deseq2_fdr", "pvalue": "deseq2_p"})
    clr = clr_results.rename(columns={"padj": "clr_fdr", "p_value": "clr_p"})

    # Merge on taxon + comparison
    merged = pd.merge(
        deseq2[["taxon", "comparison", "log2FoldChange", "deseq2_p", "deseq2_fdr"]],
        clr[["taxon", "comparison", "clr_difference", "clr_p", "clr_fdr", "direction"]],
        on=["taxon", "comparison"],
        how="inner"
    )

    # Determine DESeq2 direction
    merged["deseq2_direction"] = merged["log2FoldChange"].apply(lambda x: "up" if x > 0 else "down")

    # Concordance criteria
    merged["both_nominal"] = (merged["deseq2_p"] < p_threshold) & (merged["clr_p"] < p_threshold)
    merged["same_direction"] = merged["deseq2_direction"] == merged["direction"]
    merged["concordant"] = merged["both_nominal"] & merged["same_direction"]

    concordant = merged[merged["concordant"]]
    print(f"Total merged pairs: {len(merged)}")
    print(f"Concordant pairs (P<{p_threshold} both + same direction): {len(concordant)}")

    return merged, concordant

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    deseq2 = pd.read_csv("../results/03_deseq2/deseq2_family_results.csv")
    clr = pd.read_csv("../results/04_clr/clr_analysis_results.csv")

    merged, concordant = assess_concordance(deseq2, clr)
    merged.to_csv(os.path.join(OUTPUT_DIR, "cross_method_merged.csv"), index=False)
    concordant.to_csv(os.path.join(OUTPUT_DIR, "concordant_pairs.csv"), index=False)

    # Highlight Enterococcaceae
    entero = concordant[concordant["taxon"].str.contains("Enterococcaceae", case=False, na=False)]
    if not entero.empty:
        print("\n=== Enterococcaceae concordance ===")
        print(entero[["taxon", "comparison", "log2FoldChange", "deseq2_p", "clr_difference", "clr_p"]].to_string())
