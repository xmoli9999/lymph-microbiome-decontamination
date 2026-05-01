#!/usr/bin/env python3
"""
Step 1: Prevalence-based decontamination + kitome catalog filtering.

Two-pronged decontamination protocol:
  1. Fisher's exact test: prevalence in 6 negative controls vs 60 biological samples
  2. Cross-referencing against published kitome catalogs (78 genera)

Any OTU flagged by either method is removed.
"""

import pandas as pd
import numpy as np
from scipy.stats import fisher_exact
import os

# ============================================================
# Configuration — update these paths for your environment
# ============================================================
OTU_TABLE = "../data/otu_table.tsv"           # rows = OTUs, cols = samples
SAMPLE_META = "../data/sample_metadata.csv"    # must have 'sample_id', 'group', 'sample_type' columns
KITOME_CATALOG = "../data/kitome_catalog.csv"  # must have 'genus' column
TAXONOMY_FILE = "../data/taxonomy.tsv"         # OTU-to-taxonomy mapping
OUTPUT_DIR = "../results/01_decontamination/"

CONTROL_LABEL = "NegControl"    # value in sample_type for negative extraction controls
P_THRESHOLD = 0.05              # Fisher's exact test threshold

# ============================================================
# Load data
# ============================================================
def load_data():
    otu = pd.read_csv(OTU_TABLE, sep="\t", index_col=0)
    meta = pd.read_csv(SAMPLE_META)
    kitome = pd.read_csv(KITOME_CATALOG)
    taxonomy = pd.read_csv(TAXONOMY_FILE, sep="\t", index_col=0)
    return otu, meta, kitome, taxonomy

# ============================================================
# Method 1: Prevalence-based filtering (Fisher's exact test)
# ============================================================
def prevalence_filter(otu, meta, p_threshold=P_THRESHOLD):
    """
    For each OTU, compare prevalence in negative controls vs biological samples
    using Fisher's exact test. Flag OTUs significantly more prevalent in controls.
    """
    control_samples = meta[meta["sample_type"] == CONTROL_LABEL]["sample_id"].tolist()
    bio_samples = meta[meta["sample_type"] != CONTROL_LABEL]["sample_id"].tolist()

    n_ctrl = len(control_samples)
    n_bio = len(bio_samples)

    flagged_otus = []
    results = []

    for otu_id in otu.index:
        ctrl_present = (otu.loc[otu_id, control_samples] > 0).sum()
        bio_present = (otu.loc[otu_id, bio_samples] > 0).sum()

        # 2x2 contingency table
        table = [
            [ctrl_present, n_ctrl - ctrl_present],
            [bio_present, n_bio - bio_present]
        ]
        odds_ratio, p_value = fisher_exact(table, alternative="greater")

        results.append({
            "otu_id": otu_id,
            "ctrl_prevalence": ctrl_present / n_ctrl,
            "bio_prevalence": bio_present / n_bio,
            "odds_ratio": odds_ratio,
            "p_value": p_value,
            "flagged": p_value < p_threshold
        })

        if p_value < p_threshold:
            flagged_otus.append(otu_id)

    return flagged_otus, pd.DataFrame(results)

# ============================================================
# Method 2: Kitome catalog cross-referencing
# ============================================================
def kitome_filter(otu, taxonomy, kitome):
    """
    Flag OTUs whose genus-level classification matches any of the 78 genera
    in published kitome catalogs [Salter 2014, Eisenhofer 2019, Karstens 2019, de Goffau 2019].
    """
    kitome_genera = set(kitome["genus"].str.strip().str.lower())
    flagged_otus = []

    for otu_id in otu.index:
        if otu_id in taxonomy.index:
            genus = str(taxonomy.loc[otu_id, "genus"]).strip().lower()
            if genus in kitome_genera:
                flagged_otus.append(otu_id)

    return flagged_otus

# ============================================================
# Apply decontamination
# ============================================================
def decontaminate(otu, meta, taxonomy, kitome):
    """
    Remove OTUs flagged by either method (union).
    Prioritize specificity over sensitivity.
    """
    prev_flagged, prev_results = prevalence_filter(otu, meta)
    kit_flagged = kitome_filter(otu, taxonomy, kitome)

    all_flagged = list(set(prev_flagged) | set(kit_flagged))

    print(f"Prevalence-flagged OTUs: {len(prev_flagged)}")
    print(f"Kitome-flagged OTUs: {len(kit_flagged)}")
    print(f"Total flagged (union): {len(all_flagged)} / {len(otu)} ({100*len(all_flagged)/len(otu):.1f}%)")

    # Remove flagged OTUs
    bio_samples = meta[meta["sample_type"] != CONTROL_LABEL]["sample_id"].tolist()
    otu_clean = otu.drop(index=all_flagged, errors="ignore")
    otu_clean = otu_clean[bio_samples]

    # Calculate contamination burden per sample type
    for stype in meta["sample_type"].unique():
        if stype == CONTROL_LABEL:
            continue
        samples = meta[meta["sample_type"] == stype]["sample_id"].tolist()
        total_reads = otu[samples].sum().sum()
        contam_reads = otu.loc[otu.index.isin(all_flagged), samples].sum().sum()
        pct = 100 * contam_reads / total_reads if total_reads > 0 else 0
        print(f"  {stype}: {pct:.1f}% reads from contaminants")

    return otu_clean, all_flagged, prev_results

# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    otu, meta, kitome, taxonomy = load_data()
    otu_clean, flagged, prev_results = decontaminate(otu, meta, taxonomy, kitome)

    # Save outputs
    otu_clean.to_csv(os.path.join(OUTPUT_DIR, "otu_table_decontaminated.tsv"), sep="\t")
    prev_results.to_csv(os.path.join(OUTPUT_DIR, "prevalence_test_results.csv"), index=False)
    pd.Series(flagged).to_csv(os.path.join(OUTPUT_DIR, "flagged_otus.csv"), index=False, header=["otu_id"])

    print(f"\nDecontaminated OTU table saved: {len(otu_clean)} OTUs retained")
