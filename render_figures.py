#!/usr/bin/env python3
"""Generate Figure 2-5 source CSVs and render publication-ready PNG/PDF files.

This plotting module uses the outputs and analysis definitions from
independent_rewrite_analysis.py. It does not alter any statistical result.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Sequence

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import independent_rewrite_analysis as analysis


def write_csv(path: Path, header, rows) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def read_csv_dict(path: Path):
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def save_figure(fig, output_dir: Path, stem: str) -> None:
    fig.savefig(output_dir / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(output_dir / f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)


def figure2(output_dir: Path) -> None:
    rows = read_csv_dict(output_dir / "Table1_decontamination_rules.csv")
    labels = {
        "prevalence": "Control prevalence",
        "catalogue": "Genus catalog",
        "union": "Union",
        "intersection": "Intersection",
    }
    sample_types = [
        ("lymph_removed_pct", "Lymph 4-24 h"),
        ("fecal_removed_pct", "Feces"),
        ("BALF_removed_pct", "BALF"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.7))
    x = np.arange(len(rows))
    for ax, (column, title) in zip(axes, sample_types):
        values = [float(row[column]) for row in rows]
        bars = ax.bar(x, values)
        ax.set_xticks(x)
        ax.set_xticklabels([labels[row["rule"]] for row in rows], rotation=35, ha="right")
        ax.set_title(title)
        ax.set_ylabel("Reads flagged (%)")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        maximum = max(values)
        margin = max(maximum * 0.12, 0.002)
        ax.set_ylim(0, maximum + margin)
        for bar, value in zip(bars, values):
            label = f"{value:.2f}" if value >= 0.01 else f"{value:.4f}"
            ax.text(
                bar.get_x() + bar.get_width() / 2,
                bar.get_height() + margin * 0.08,
                label,
                ha="center",
                va="bottom",
                fontsize=8,
            )
    fig.suptitle("Reads flagged under four decontamination rules", y=1.02)
    fig.tight_layout()
    save_figure(fig, output_dir, "Figure2_decontamination_rules")


def figure3(output_dir: Path) -> None:
    rows = read_csv_dict(output_dir / "Fig3_baseline_composition.csv")
    group_labels = {
        "baseline_0h": "0 h lymph",
        "negative_controls": "Negative controls",
    }
    groups = [group_labels[row["group"]] for row in rows]
    gut = np.array([100 * float(row["gut_associated_fraction"]) for row in rows])
    reagent = np.array([100 * float(row["reagent_catalogue_fraction"]) for row in rows])
    other = 100 - gut - reagent

    fig, ax = plt.subplots(figsize=(6.2, 4.5))
    x = np.arange(len(groups))
    ax.bar(x, gut, label="Gut-associated genera")
    ax.bar(x, reagent, bottom=gut, label="Reagent/environment catalogue")
    ax.bar(x, other, bottom=gut + reagent, label="Other taxa")
    ax.set_xticks(x)
    ax.set_xticklabels(groups)
    ax.set_ylabel("Pooled read composition (%)")
    ax.set_ylim(0, 100)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
    for index in range(len(groups)):
        ax.text(index, max(gut[index] / 2, 1.2), f"{gut[index]:.2f}%", ha="center", va="center", fontsize=8)
        if reagent[index] > 3:
            ax.text(index, gut[index] + reagent[index] / 2, f"{reagent[index]:.1f}%", ha="center", va="center", fontsize=8)
    fig.tight_layout()
    save_figure(fig, output_dir, "Figure3_baseline_composition")


def create_fig4_source(data, output_dir: Path):
    controls = [item.index for item in data.sample_info if item.is_control]
    baseline = sorted(
        [item.index for item in data.sample_info if item.is_baseline],
        key=lambda index: int(data.sample_info[index].animal),
    )
    indices = baseline + controls
    labels = ["0 h lymph"] * len(baseline) + ["Negative control"] * len(controls)
    no_removal = np.zeros(data.counts.shape[0], dtype=bool)
    relative = analysis.family_relative_abundance(data, no_removal, indices)
    distance = analysis.bray_curtis_matrix(relative)

    n = distance.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centering @ (distance ** 2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gower)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > 1e-12
    positive_values = eigenvalues[positive]
    coordinates = eigenvectors[:, positive] * np.sqrt(positive_values)
    explained = positive_values / positive_values.sum() * 100

    rows = []
    for pos, sample_index in enumerate(indices):
        info = data.sample_info[sample_index]
        rows.append(
            (
                info.sample_id,
                labels[pos],
                info.animal,
                "" if info.timepoint is None else info.timepoint,
                coordinates[pos, 0],
                coordinates[pos, 1],
                explained[0],
                explained[1],
            )
        )
    write_csv(
        output_dir / "Fig4_pcoa_coordinates.csv",
        ["sample_id", "group", "animal", "timepoint_h", "PC1", "PC2", "PC1_variance_pct", "PC2_variance_pct"],
        rows,
    )
    return rows


def figure4(output_dir: Path, rows) -> None:
    groups = sorted(set(row[1] for row in rows))
    markers = {"0 h lymph": "o", "Negative control": "s"}
    fig, ax = plt.subplots(figsize=(6.2, 5.0))
    for group in groups:
        selected = [row for row in rows if row[1] == group]
        ax.scatter(
            [float(row[4]) for row in selected],
            [float(row[5]) for row in selected],
            label=group,
            marker=markers[group],
            s=62,
        )
        for row in selected:
            ax.annotate(row[0], (float(row[4]), float(row[5])), xytext=(4, 4), textcoords="offset points", fontsize=7)
    pc1_var = float(rows[0][6])
    pc2_var = float(rows[0][7])
    ax.set_xlabel(f"PCoA 1 ({pc1_var:.1f}%)")
    ax.set_ylabel(f"PCoA 2 ({pc2_var:.1f}%)")
    ax.axhline(0, linewidth=0.7, alpha=0.35)
    ax.axvline(0, linewidth=0.7, alpha=0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False)
    ax.set_title("Bray-Curtis ordination: baseline lymph and controls")
    fig.tight_layout()
    save_figure(fig, output_dir, "Figure4_pcoa")


def create_fig5_source(data, output_dir: Path):
    controls = [item.index for item in data.sample_info if item.is_control]
    primary = [item.index for item in data.sample_info if item.is_primary_lymph]
    fecal = [item.index for item in data.sample_info if item.sample_type == "FB"]
    balf = [item.index for item in data.sample_info if item.sample_type == "BALF"]
    biological = primary + fecal + balf

    prevalence_mask, _, _ = analysis.prevalence_screen(data, controls, biological)
    catalogue_mask = np.asarray(
        [bool(genus) and genus.lower() in analysis.CONTAMINANT_GENERA_LOWER for genus in data.genus],
        dtype=bool,
    )
    union_mask = prevalence_mask | catalogue_mask
    main_da = analysis.family_da(data, union_mask, primary, include_animal=True)
    ent_row = analysis.find_row(main_da, "Enterococcaceae")
    ent_counts = np.asarray(main_da["counts"])[ent_row]
    ent_clr = np.asarray(main_da["clr"])[ent_row]
    totals = np.asarray(main_da["counts"]).sum(axis=0)

    rows = []
    for position, sample_index in enumerate(primary):
        info = data.sample_info[sample_index]
        relative = float(ent_counts[position] / totals[position]) if totals[position] else 0.0
        rows.append(
            (
                info.sample_id,
                info.animal,
                info.timepoint,
                int(ent_counts[position]),
                int(totals[position]),
                relative,
                float(ent_clr[position]),
            )
        )
    rows.sort(key=lambda row: (int(row[1]), int(row[2])))
    write_csv(
        output_dir / "Fig5_enterococcaceae_trajectory.csv",
        ["sample_id", "animal", "timepoint_h", "raw_family_count", "total_retained_family_reads", "relative_abundance", "clr_value"],
        rows,
    )
    return rows


def figure5(output_dir: Path, rows) -> None:
    animals = sorted(set(row[1] for row in rows), key=int)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for animal in animals:
        selected = sorted([row for row in rows if row[1] == animal], key=lambda row: int(row[2]))
        ax.plot(
            [int(row[2]) for row in selected],
            [float(row[6]) for row in selected],
            marker="o",
            linewidth=1.4,
            label=f"Animal {animal}",
        )
    ax.set_xticks([4, 8, 12, 24])
    ax.set_xlabel("Time after CLP (h)")
    ax.set_ylabel("Enterococcaceae CLR abundance")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, ncol=3)
    ax.set_title("Enterococcaceae trajectories after decontamination")
    fig.tight_layout()
    save_figure(fig, output_dir, "Figure5_enterococcaceae")


def render_all(otu_workbook: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    data = analysis.load_otu_workbook(otu_workbook)
    figure2(output_dir)
    figure3(output_dir)
    fig4_rows = create_fig4_source(data, output_dir)
    figure4(output_dir, fig4_rows)
    fig5_rows = create_fig5_source(data, output_dir)
    figure5(output_dir, fig5_rows)


def parse_arguments(argv: Sequence[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("otu_workbook", type=Path)
    parser.add_argument("--output-dir", type=Path, default=Path("outputs"))
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_arguments(argv)
    render_all(args.otu_workbook.resolve(), args.output_dir.resolve())
    print(f"Figure 2-5 source CSVs and PNG/PDF files written to: {args.output_dir.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
