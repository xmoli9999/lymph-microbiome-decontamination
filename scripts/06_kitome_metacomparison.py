#!/usr/bin/env python3
"""
Step 6: Kitome meta-comparison.

Compare identified contaminant genera against combined published kitome catalogs:
  - Salter et al. 2014
  - Eisenhofer et al. 2019
  - Karstens et al. 2019
  - de Goffau et al. 2019
"""

import pandas as pd
import os

OUTPUT_DIR = "../results/06_kitome/"

def kitome_metacomparison(our_contaminants, kitome_catalog):
    """Compare our contaminant genera against published kitome references."""
    our_genera = set(our_contaminants["genus"].str.strip().str.lower())
    pub_genera = set(kitome_catalog["genus"].str.strip().str.lower())

    overlap = our_genera & pub_genera
    our_unique = our_genera - pub_genera
    pub_unique = pub_genera - our_genera

    print(f"Our contaminant genera: {len(our_genera)}")
    print(f"Published kitome genera: {len(pub_genera)}")
    print(f"Overlap: {len(overlap)} ({100*len(overlap)/len(our_genera):.1f}% of ours, "
          f"{100*len(overlap)/len(pub_genera):.1f}% of published)")
    print(f"Our unique: {len(our_unique)}")
    print(f"Published unique: {len(pub_unique)}")

    results = {
        "overlap": sorted(overlap),
        "our_unique": sorted(our_unique),
        "published_unique": sorted(pub_unique),
        "summary": {
            "n_our": len(our_genera),
            "n_published": len(pub_genera),
            "n_overlap": len(overlap),
        }
    }
    return results

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load our flagged contaminant genera
    flagged_otus = pd.read_csv("../results/01_decontamination/flagged_otus.csv")
    taxonomy = pd.read_csv("../data/taxonomy.tsv", sep="\t", index_col=0)

    our_contaminants = taxonomy.loc[
        taxonomy.index.isin(flagged_otus["otu_id"]), ["genus"]
    ].dropna().drop_duplicates()

    kitome = pd.read_csv("../data/kitome_catalog.csv")

    results = kitome_metacomparison(our_contaminants, kitome)

    # Save
    pd.DataFrame({"genus": results["overlap"]}).to_csv(
        os.path.join(OUTPUT_DIR, "overlapping_genera.csv"), index=False
    )
    pd.DataFrame({"genus": results["our_unique"]}).to_csv(
        os.path.join(OUTPUT_DIR, "our_unique_genera.csv"), index=False
    )

    # Per-reference breakdown (if catalog has source column)
    if "source" in kitome.columns:
        for src in kitome["source"].unique():
            src_genera = set(kitome[kitome["source"] == src]["genus"].str.strip().str.lower())
            src_overlap = set(our_contaminants["genus"].str.strip().str.lower()) & src_genera
            print(f"  {src}: {len(src_overlap)}/{len(src_genera)} overlap")
