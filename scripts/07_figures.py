#!/usr/bin/env python3
"""
Step 7: Figure generation for manuscript.

Figure 1: Decontamination impact (contamination proportions, Shannon diversity, PCoA)
Figure 2: Volcano plots (genus and family level)
Figure 3: Temporal heatmaps (CLR-transformed abundance)
Figure 5: Kitome Venn diagram
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2
import os

OUTPUT_DIR = "../results/figures/"
plt.rcParams.update({"font.size": 12, "font.family": "Arial"})

def fig1_contamination_impact():
    """Figure 1: Decontamination impact on mesenteric lymph microbiome."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel a: Contamination proportions by sample type
    contam_data = {"Lymph": 49.2, "BALF": 35.8, "Fecal": 12.3}
    ax = axes[0]
    bars = ax.bar(contam_data.keys(), contam_data.values(), color=["#E74C3C", "#F39C12", "#2ECC71"])
    ax.set_ylabel("Contamination (%)")
    ax.set_title("a  Contamination by sample type")
    ax.set_ylim(0, 60)
    for bar, val in zip(bars, contam_data.values()):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1, f"{val}%", ha="center")

    # Panel b: Shannon diversity (placeholder - replace with actual data)
    ax = axes[1]
    ax.set_title("b  Shannon diversity after decontamination")
    ax.set_xlabel("Time point")
    ax.set_ylabel("Shannon index")
    ax.text(0.5, 0.5, "Load from\nalpha_diversity.csv", transform=ax.transAxes, ha="center", va="center")

    # Panel c: PCoA (placeholder)
    ax = axes[2]
    ax.set_title("c  PCoA (Bray-Curtis)")
    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.text(0.5, 0.5, "Load from\npcoa_coordinates.csv", transform=ax.transAxes, ha="center", va="center")

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure1_Decontamination_Impact.png"), dpi=300, bbox_inches="tight")
    plt.close()

def fig2_volcano(deseq2_results, level="family"):
    """Figure 2: Volcano plots of DESeq2 differential abundance."""
    comparisons = deseq2_results["comparison"].unique()
    n_comp = len(comparisons)
    fig, axes = plt.subplots(1, n_comp, figsize=(5*n_comp, 5))
    if n_comp == 1:
        axes = [axes]

    for ax, comp in zip(axes, comparisons):
        data = deseq2_results[deseq2_results["comparison"] == comp].copy()
        data["-log10p"] = -np.log10(data["pvalue"].clip(lower=1e-10))

        # Color by significance
        colors = []
        for _, row in data.iterrows():
            if row["padj"] < 0.05:
                colors.append("#E74C3C")  # FDR significant
            elif row["pvalue"] < 0.05:
                colors.append("#F39C12")  # Nominal
            else:
                colors.append("#BDC3C7")  # NS

        ax.scatter(data["log2FoldChange"], data["-log10p"], c=colors, alpha=0.7, s=30)
        ax.axhline(-np.log10(0.05), ls="--", color="gray", alpha=0.5)
        ax.set_xlabel("log₂ Fold Change")
        ax.set_ylabel("-log₁₀(P)")
        ax.set_title(comp.replace("_", " "))

    plt.suptitle(f"DESeq2 — {level} level", fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"Figure2_Volcano_{level.capitalize()}.png"), dpi=300, bbox_inches="tight")
    plt.close()

def fig3_heatmap(clr_results, meta, level="family"):
    """Figure 3: Temporal heatmap of CLR-transformed abundance."""
    # Pivot to taxon x time_point matrix (mean CLR difference)
    pivot = clr_results.pivot_table(
        index="taxon", columns="comparison", values="clr_difference", aggfunc="mean"
    )

    # Select top 20 most variable
    pivot_var = pivot.var(axis=1).sort_values(ascending=False)
    top_taxa = pivot_var.head(20).index
    plot_data = pivot.loc[top_taxa]

    fig, ax = plt.subplots(figsize=(10, 8))
    sns.heatmap(plot_data, cmap="RdBu_r", center=0, ax=ax,
                xticklabels=True, yticklabels=True, linewidths=0.5)
    ax.set_title(f"CLR-transformed temporal dynamics — {level} level")
    ax.set_xlabel("Comparison")
    ax.set_ylabel("")

    plt.tight_layout()
    plt.savefig(os.path.join(OUTPUT_DIR, f"Figure3_Heatmap_{level.capitalize()}.png"), dpi=300, bbox_inches="tight")
    plt.close()

def fig5_kitome_venn():
    """Figure 5: Kitome meta-comparison Venn diagram."""
    fig, ax = plt.subplots(figsize=(8, 6))
    venn2(
        subsets=(120, 53, 25),  # our_unique, pub_unique, overlap
        set_labels=("This study\n(n=145)", "Published references\n(n=78)"),
        ax=ax
    )
    ax.set_title("Kitome meta-comparison")
    plt.savefig(os.path.join(OUTPUT_DIR, "Figure5_Kitome_Venn.png"), dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Generate figures with available data
    fig1_contamination_impact()

    # Load DESeq2 and CLR results if available
    deseq2_fam = "../results/03_deseq2/deseq2_family_results.csv"
    clr_fam = "../results/04_clr/clr_analysis_results.csv"

    if os.path.exists(deseq2_fam):
        fig2_volcano(pd.read_csv(deseq2_fam), level="family")

    if os.path.exists(clr_fam):
        meta = pd.read_csv("../data/sample_metadata.csv")
        fig3_heatmap(pd.read_csv(clr_fam), meta, level="family")

    fig5_kitome_venn()
    print(f"Figures saved to {OUTPUT_DIR}")
