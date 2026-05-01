#!/usr/bin/env python3
"""
Step 3: DESeq2 differential abundance analysis using pyDESeq2.

Pairwise comparisons: 0h vs 4h, 8h, 12h, 24h
Taxonomic levels: genus and family
FDR correction: Benjamini-Hochberg
"""

import pandas as pd
import numpy as np
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
import os
import warnings
warnings.filterwarnings("ignore")

OUTPUT_DIR = "../results/03_deseq2/"

def aggregate_taxonomy(otu, taxonomy, level="family"):
    """Aggregate OTU counts to specified taxonomic level."""
    tax_map = taxonomy[level].fillna("Unclassified")
    otu_with_tax = otu.copy()
    otu_with_tax["taxon"] = otu_with_tax.index.map(lambda x: tax_map.get(x, "Unclassified"))
    aggregated = otu_with_tax.groupby("taxon").sum()
    return aggregated

def prevalence_filter_taxa(counts, min_reads=2, min_samples=3):
    """Require >= min_reads in >= min_samples."""
    mask = (counts >= min_reads).sum(axis=1) >= min_samples
    return counts[mask]

def run_deseq2_pairwise(counts, meta, ref_time="0h", test_time="4h"):
    """Run DESeq2 for a single pairwise comparison."""
    # Subset to relevant samples
    relevant = meta[meta["time_point"].isin([ref_time, test_time])]
    samples = [s for s in relevant["sample_id"] if s in counts.columns]
    sub_counts = counts[samples].T
    sub_meta = relevant.set_index("sample_id").loc[samples]

    # Ensure integer counts
    sub_counts = sub_counts.astype(int)

    # Run DESeq2
    dds = DeseqDataSet(
        counts=sub_counts,
        metadata=sub_meta,
        design="~time_point",
        refit_cooks=True
    )
    dds.deseq2()

    # Extract results
    stat = DeseqStats(dds, contrast=["time_point", test_time, ref_time])
    stat.summary()
    results = stat.results_df.copy()
    results["comparison"] = f"{ref_time}_vs_{test_time}"

    return results

def run_all_comparisons(counts, meta, level_name):
    """Run all pairwise comparisons against 0h."""
    time_points = ["4h", "8h", "12h", "24h"]
    all_results = []

    for tp in time_points:
        try:
            res = run_deseq2_pairwise(counts, meta, ref_time="0h", test_time=tp)
            res["level"] = level_name
            all_results.append(res)
            sig = res[res["padj"] < 0.05]
            print(f"  {level_name} 0h vs {tp}: {len(sig)} FDR-significant taxa")
        except Exception as e:
            print(f"  {level_name} 0h vs {tp}: Error - {e}")

    return pd.concat(all_results, ignore_index=True) if all_results else pd.DataFrame()

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load data
    otu = pd.read_csv("../results/01_decontamination/otu_table_decontaminated.tsv", sep="\t", index_col=0)
    meta = pd.read_csv("../data/sample_metadata.csv")
    taxonomy = pd.read_csv("../data/taxonomy.tsv", sep="\t", index_col=0)
    bio_meta = meta[meta["sample_type"] != "NegControl"]

    # Genus level
    print("=== Genus-level analysis ===")
    genus_counts = aggregate_taxonomy(otu, taxonomy, level="genus")
    genus_filtered = prevalence_filter_taxa(genus_counts)
    genus_results = run_all_comparisons(genus_filtered, bio_meta, "genus")
    genus_results.to_csv(os.path.join(OUTPUT_DIR, "deseq2_genus_results.csv"), index=False)

    # Family level
    print("\n=== Family-level analysis ===")
    family_counts = aggregate_taxonomy(otu, taxonomy, level="family")
    family_filtered = prevalence_filter_taxa(family_counts)
    family_results = run_all_comparisons(family_filtered, bio_meta, "family")
    family_results.to_csv(os.path.join(OUTPUT_DIR, "deseq2_family_results.csv"), index=False)

    print(f"\nResults saved to {OUTPUT_DIR}")
