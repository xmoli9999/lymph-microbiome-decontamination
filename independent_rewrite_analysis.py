#!/usr/bin/env python3
"""Independent downstream re-analysis of the 8,042-OTU table.

This program was written from the statistical specification only.  It does not
import or call any pre-existing analysis script.  It reads the supplied .xlsx
file directly through the OOXML container using the Python standard library,
then implements the requested QC, decontamination, family-level CLR models,
PERMANOVA, alpha diversity, ablation and sensitivity analyses with NumPy/SciPy.

Usage
-----
python independent_rewrite_analysis.py \
    /path/to/otu_table2_with_taxonomy.xlsx \
    --output-dir independent_rewrite_out

Deterministic settings
----------------------
PERMANOVA seed = 1; rarefaction seed = 42; permutations = 1,999.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import re
import sys
import zipfile
from collections import OrderedDict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence
from xml.etree import ElementTree as ET

import numpy as np
import scipy
from scipy import stats


# ---------------------------------------------------------------------------
# Fixed analysis definitions
# ---------------------------------------------------------------------------

DEPTH_THRESHOLD = 5_000
PERMUTATIONS = 1_999
PERMANOVA_SEED = 1
RAREFACTION_SEED = 42
PSEUDOCOUNT = 0.5

CONTAMINANT_GENERA = tuple(
    x.strip()
    for x in """
Acetobacter, Acinetobacter, Afipia, Agrobacterium, Bradyrhizobium, Brevundimonas,
Burkholderia, Burkholderia-Caballeronia-Paraburkholderia, Caulobacter, Chryseobacterium,
Comamonas, Cutibacterium, Delftia, Devosia, Herbaspirillum, Massilia, Methylobacterium,
Methylophilus, Microbacterium, Novosphingobium, Oxalobacter, Pelomonas, Phyllobacterium,
Propionibacterium, Pseudomonas, Pseudoxanthomonas, Ralstonia, Rhizobium, Rhodococcus,
Solirubrobacter, Sphingobium, Sphingomonas, Stenotrophomonas, Variovorax
""".replace("\n", " ").split(",")
    if x.strip()
)
CONTAMINANT_GENERA_LOWER = frozenset(x.lower() for x in CONTAMINANT_GENERA)

GUT_GENERA = tuple(
    x.strip()
    for x in """
Lactobacillus, Bacteroides, Prevotella, Roseburia, Ruminococcus, Lachnoclostridium,
Blautia, Faecalibacterium, Romboutsia, Turicibacter, Colidextribacter, UCG-005
""".replace("\n", " ").split(",")
    if x.strip()
)
GUT_GENERA_LOWER = frozenset(x.lower() for x in GUT_GENERA)

S4_SOURCES = (
    "Salter et al. 2014; Eisenhofer et al. 2019; "
    "Karstens et al. 2019; de Goffau et al. 2019"
)


# ---------------------------------------------------------------------------
# Minimal OOXML reader (no pandas/openpyxl dependency)
# ---------------------------------------------------------------------------

_SS_NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
_ODREL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def _column_index(cell_reference: str) -> int:
    match = re.match(r"([A-Z]+)", cell_reference)
    if not match:
        raise ValueError(f"Invalid Excel cell reference: {cell_reference!r}")
    value = 0
    for letter in match.group(1):
        value = value * 26 + ord(letter) - 64
    return value - 1


def read_first_worksheet(path: Path) -> list[list[object]]:
    """Read the first worksheet into a rectangular list of rows."""
    with zipfile.ZipFile(path) as archive:
        names = set(archive.namelist())
        shared_strings: list[str] = []
        if "xl/sharedStrings.xml" in names:
            root = ET.fromstring(archive.read("xl/sharedStrings.xml"))
            for item in root.findall(_SS_NS + "si"):
                shared_strings.append(
                    "".join((node.text or "") for node in item.iter(_SS_NS + "t"))
                )

        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {rel.attrib["Id"]: rel.attrib["Target"] for rel in relationships}
        sheets = workbook.find(_SS_NS + "sheets")
        if sheets is None or len(sheets) == 0:
            raise ValueError("Workbook contains no worksheets")
        first_sheet = sheets[0]
        relationship_id = first_sheet.attrib[_ODREL_NS + "id"]
        target = targets[relationship_id]
        if not target.startswith("xl/"):
            target = "xl/" + target.lstrip("/")

        worksheet = ET.fromstring(archive.read(target))
        sparse_rows: list[dict[int, object]] = []
        maximum_column = 0
        for row in worksheet.iter(_SS_NS + "row"):
            values: dict[int, object] = {}
            for cell in row.findall(_SS_NS + "c"):
                column = _column_index(cell.attrib["r"])
                maximum_column = max(maximum_column, column + 1)
                cell_type = cell.attrib.get("t")
                value_element = cell.find(_SS_NS + "v")
                if cell_type == "inlineStr":
                    value: object = "".join(
                        (node.text or "") for node in cell.iter(_SS_NS + "t")
                    )
                elif value_element is None or value_element.text is None:
                    value = ""
                elif cell_type == "s":
                    value = shared_strings[int(value_element.text)]
                elif cell_type == "b":
                    value = bool(int(value_element.text))
                else:
                    raw = value_element.text
                    try:
                        numeric = float(raw)
                        value = int(numeric) if numeric.is_integer() else numeric
                    except ValueError:
                        value = raw
                values[column] = value
            sparse_rows.append(values)

    return [
        [sparse.get(column, "") for column in range(maximum_column)]
        for sparse in sparse_rows
    ]


# ---------------------------------------------------------------------------
# Data model and parsing
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class SampleInfo:
    index: int
    sample_id: str
    sample_type: str
    animal: str
    timepoint: int | None
    depth: int
    is_control: bool
    is_baseline: bool
    is_primary_lymph: bool
    excluded_reason: str


@dataclass
class OtuData:
    otu_ids: np.ndarray
    sample_ids: list[str]
    counts: np.ndarray  # shape: OTU x sample
    taxonomy: np.ndarray
    family: np.ndarray
    genus: np.ndarray
    sample_info: list[SampleInfo]


def taxonomy_rank(taxonomy: str, prefix: str) -> str:
    for component in taxonomy.split(";"):
        component = component.strip()
        if component.lower().startswith(prefix.lower()):
            return component[len(prefix) :].strip()
    return ""


def sample_metadata(sample_id: str, index: int, depth: int) -> SampleInfo:
    is_control = sample_id.startswith("con_") or sample_id.lower() == "water"
    lymph_match = re.fullmatch(r"([^_]+)_(\d+)hLB", sample_id)
    if is_control:
        sample_type, animal, timepoint = "control", "", None
    elif sample_id.endswith("_FB"):
        sample_type, animal, timepoint = "FB", sample_id.split("_", 1)[0], None
    elif sample_id.endswith("_FP"):
        sample_type, animal, timepoint = "BALF", sample_id.split("_", 1)[0], None
    elif lymph_match:
        sample_type = "lymph"
        animal = lymph_match.group(1)
        timepoint = int(lymph_match.group(2))
    else:
        sample_type, animal, timepoint = "other", "", None

    is_baseline = sample_type == "lymph" and timepoint == 0
    post_clp = sample_type == "lymph" and timepoint in {4, 8, 12, 24}
    low_depth = post_clp and depth < DEPTH_THRESHOLD
    is_primary = bool(post_clp and not low_depth)
    excluded_reason = ""
    if low_depth:
        excluded_reason = f"depth<{DEPTH_THRESHOLD}"
    elif is_baseline:
        excluded_reason = "0h baseline excluded from primary analysis"

    return SampleInfo(
        index=index,
        sample_id=sample_id,
        sample_type=sample_type,
        animal=animal,
        timepoint=timepoint,
        depth=depth,
        is_control=is_control,
        is_baseline=is_baseline,
        is_primary_lymph=is_primary,
        excluded_reason=excluded_reason,
    )


def load_otu_workbook(path: Path) -> OtuData:
    rows = read_first_worksheet(path)
    if len(rows) < 2:
        raise ValueError("The workbook does not contain OTU data")
    header = [str(value).strip() for value in rows[0]]
    taxonomy_candidates = [
        index for index, value in enumerate(header) if "taxonomy" in value.lower()
    ]
    if len(taxonomy_candidates) != 1:
        raise ValueError(
            f"Expected one taxonomy column, found {len(taxonomy_candidates)}"
        )
    taxonomy_column = taxonomy_candidates[0]
    otu_column = 0
    sample_columns = list(range(otu_column + 1, taxonomy_column))
    sample_ids = [header[column] for column in sample_columns]

    otu_ids: list[str] = []
    taxonomy: list[str] = []
    count_rows: list[list[int]] = []
    for row_number, row in enumerate(rows[1:], start=2):
        if not any(str(value).strip() for value in row):
            continue
        otu_id = str(row[otu_column]).strip()
        if not otu_id:
            raise ValueError(f"Missing OTU ID in worksheet row {row_number}")
        converted: list[int] = []
        for column in sample_columns:
            value = row[column]
            if value in ("", None):
                converted.append(0)
            else:
                numeric = float(value)
                if numeric < 0 or not numeric.is_integer():
                    raise ValueError(
                        f"Non-integer or negative count at row {row_number}, "
                        f"sample {header[column]}: {value!r}"
                    )
                converted.append(int(numeric))
        otu_ids.append(otu_id)
        count_rows.append(converted)
        taxonomy.append(str(row[taxonomy_column]).strip())

    counts = np.asarray(count_rows, dtype=np.int64)
    taxonomy_array = np.asarray(taxonomy, dtype=object)
    families = np.asarray(
        [taxonomy_rank(value, "f__") or "unclassified_family" for value in taxonomy],
        dtype=object,
    )
    genera = np.asarray(
        [taxonomy_rank(value, "g__") for value in taxonomy], dtype=object
    )
    depths = counts.sum(axis=0)
    metadata = [
        sample_metadata(sample_id, index, int(depths[index]))
        for index, sample_id in enumerate(sample_ids)
    ]

    return OtuData(
        otu_ids=np.asarray(otu_ids, dtype=object),
        sample_ids=sample_ids,
        counts=counts,
        taxonomy=taxonomy_array,
        family=families,
        genus=genera,
        sample_info=metadata,
    )


# ---------------------------------------------------------------------------
# Generic statistics
# ---------------------------------------------------------------------------


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    m = p.size
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * m / np.arange(1, m + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    result = np.empty(m, dtype=float)
    result[order] = adjusted
    return result


def design_matrix(
    animals: Sequence[str] | None,
    groups: Sequence[object],
    include_animal: bool,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    groups_array = np.asarray(groups, dtype=object)
    n = len(groups_array)
    reduced_columns = [np.ones(n, dtype=float)]
    if include_animal:
        if animals is None:
            raise ValueError("Animal labels are required for animal-adjusted models")
        animals_array = np.asarray(animals, dtype=object)
        animal_levels = sorted(set(animals_array), key=lambda value: int(value))
        for level in animal_levels[1:]:
            reduced_columns.append((animals_array == level).astype(float))
    reduced = np.column_stack(reduced_columns)

    full_columns = [reduced]
    group_levels = sorted(
        set(groups_array), key=lambda value: (isinstance(value, str), value)
    )
    for level in group_levels[1:]:
        full_columns.append((groups_array == level).astype(float)[:, None])
    full = np.column_stack(full_columns)
    df_numerator = np.linalg.matrix_rank(full) - np.linalg.matrix_rank(reduced)
    df_denominator = n - np.linalg.matrix_rank(full)
    return reduced, full, int(df_numerator), int(df_denominator)


def nested_f_test_matrix(
    y_matrix: np.ndarray,
    groups: Sequence[object],
    animals: Sequence[str] | None = None,
    include_animal: bool = True,
) -> tuple[np.ndarray, np.ndarray, int, int]:
    """Fit reduced/full least-squares models for each row of y_matrix."""
    reduced, full, df_num, df_den = design_matrix(
        animals=animals, groups=groups, include_animal=include_animal
    )
    if df_num <= 0 or df_den <= 0:
        raise ValueError(
            f"Invalid model degrees of freedom: numerator={df_num}, denominator={df_den}"
        )
    # y_matrix is feature x sample.  np.linalg.lstsq expects sample x feature.
    y = np.asarray(y_matrix, dtype=float).T
    beta_reduced = np.linalg.lstsq(reduced, y, rcond=None)[0]
    beta_full = np.linalg.lstsq(full, y, rcond=None)[0]
    residual_reduced = y - reduced @ beta_reduced
    residual_full = y - full @ beta_full
    sse_reduced = np.sum(residual_reduced * residual_reduced, axis=0)
    sse_full = np.sum(residual_full * residual_full, axis=0)
    numerator = np.maximum(sse_reduced - sse_full, 0.0) / df_num
    denominator = sse_full / df_den
    f_values = np.divide(
        numerator,
        denominator,
        out=np.zeros_like(numerator),
        where=denominator > 0,
    )
    p_values = stats.f.sf(f_values, df_num, df_den)
    return f_values, p_values, df_num, df_den


def clr_transform(counts: np.ndarray, pseudocount: float = PSEUDOCOUNT) -> np.ndarray:
    values = np.log(np.asarray(counts, dtype=float) + pseudocount)
    return values - values.mean(axis=0, keepdims=True)


def aggregate_by_family(
    counts: np.ndarray, families: np.ndarray, retain_otu_mask: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    retained_families = families[retain_otu_mask]
    levels = np.asarray(sorted(set(retained_families)), dtype=object)
    row_for_family = {family: index for index, family in enumerate(levels)}
    aggregated = np.zeros((len(levels), counts.shape[1]), dtype=np.int64)
    for otu_index in np.flatnonzero(retain_otu_mask):
        aggregated[row_for_family[families[otu_index]], :] += counts[otu_index, :]
    return levels, aggregated


def family_da(
    data: OtuData,
    remove_mask: np.ndarray,
    sample_indices: Sequence[int],
    include_animal: bool,
) -> dict[str, object]:
    retain = ~np.asarray(remove_mask, dtype=bool)
    family_names, family_counts_all = aggregate_by_family(
        data.counts, data.family, retain
    )
    sample_indices_array = np.asarray(sample_indices, dtype=int)
    family_counts = family_counts_all[:, sample_indices_array]
    prevalence = np.sum(family_counts > 0, axis=1)
    keep_family = prevalence >= 4
    tested_names = family_names[keep_family]
    tested_counts = family_counts[keep_family, :]
    clr = clr_transform(tested_counts)
    metadata = [data.sample_info[index] for index in sample_indices_array]
    groups = [item.timepoint for item in metadata]
    animals = [item.animal for item in metadata]
    f_values, p_values, df_num, df_den = nested_f_test_matrix(
        clr,
        groups=groups,
        animals=animals if include_animal else None,
        include_animal=include_animal,
    )
    q_values = benjamini_hochberg(p_values)
    return {
        "family": tested_names,
        "counts": tested_counts,
        "clr": clr,
        "F": f_values,
        "p": p_values,
        "q": q_values,
        "df_num": df_num,
        "df_den": df_den,
        "n_samples": len(sample_indices_array),
    }


def find_row(result: dict[str, object], family_name: str) -> int:
    names = np.asarray(result["family"], dtype=object)
    matches = np.flatnonzero(names == family_name)
    if matches.size != 1:
        raise KeyError(f"Expected one row for {family_name!r}; found {matches.size}")
    return int(matches[0])


# ---------------------------------------------------------------------------
# Decontamination rules
# ---------------------------------------------------------------------------


def prevalence_screen(
    data: OtuData,
    reference_indices: Sequence[int],
    biological_indices: Sequence[int],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    references = data.counts[:, np.asarray(reference_indices, dtype=int)] > 0
    biological = data.counts[:, np.asarray(biological_indices, dtype=int)] > 0
    ref_present = references.sum(axis=1)
    bio_present = biological.sum(axis=1)
    number_ref = len(reference_indices)
    number_bio = len(biological_indices)
    p_values = np.empty(data.counts.shape[0], dtype=float)
    for otu in range(data.counts.shape[0]):
        table = [
            [int(ref_present[otu]), int(number_ref - ref_present[otu])],
            [int(bio_present[otu]), int(number_bio - bio_present[otu])],
        ]
        p_values[otu] = stats.fisher_exact(table, alternative="greater").pvalue
    q_values = benjamini_hochberg(p_values)
    return p_values < 0.05, p_values, q_values


def removed_read_fraction(
    data: OtuData, remove_mask: np.ndarray, sample_indices: Sequence[int]
) -> float:
    columns = np.asarray(sample_indices, dtype=int)
    denominator = float(data.counts[:, columns].sum())
    numerator = float(data.counts[np.asarray(remove_mask, dtype=bool)][:, columns].sum())
    return numerator / denominator if denominator else math.nan


# ---------------------------------------------------------------------------
# Bray-Curtis, PERMANOVA and PERMDISP
# ---------------------------------------------------------------------------


def family_relative_abundance(
    data: OtuData, remove_mask: np.ndarray, sample_indices: Sequence[int]
) -> np.ndarray:
    family_names, family_counts = aggregate_by_family(
        data.counts, data.family, ~np.asarray(remove_mask, dtype=bool)
    )
    del family_names
    matrix = family_counts[:, np.asarray(sample_indices, dtype=int)].T.astype(float)
    row_sums = matrix.sum(axis=1, keepdims=True)
    if np.any(row_sums == 0):
        raise ValueError("At least one sample has zero retained reads")
    return matrix / row_sums


def bray_curtis_matrix(relative_abundance: np.ndarray) -> np.ndarray:
    x = np.asarray(relative_abundance, dtype=float)
    n = x.shape[0]
    distances = np.zeros((n, n), dtype=float)
    for first in range(n):
        for second in range(first):
            denominator = np.sum(x[first] + x[second])
            distance = (
                np.sum(np.abs(x[first] - x[second])) / denominator
                if denominator > 0
                else 0.0
            )
            distances[first, second] = distances[second, first] = distance
    return distances


def permanova_statistic(
    distance_matrix: np.ndarray, labels: Sequence[object]
) -> tuple[float, float]:
    distances = np.asarray(distance_matrix, dtype=float)
    labels_array = np.asarray(labels, dtype=object)
    n = len(labels_array)
    levels = sorted(set(labels_array), key=lambda value: (isinstance(value, str), value))
    design_columns = [np.ones(n, dtype=float)]
    for level in levels[1:]:
        design_columns.append((labels_array == level).astype(float))
    design = np.column_stack(design_columns)
    projection = design @ np.linalg.pinv(design.T @ design) @ design.T
    centering = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centering @ (distances**2) @ centering
    ss_model = float(np.trace(projection @ gower))
    ss_total = float(np.trace(gower))
    ss_residual = ss_total - ss_model
    df_model = len(levels) - 1
    df_residual = n - len(levels)
    pseudo_f = (ss_model / df_model) / (ss_residual / df_residual)
    r_squared = ss_model / ss_total
    return pseudo_f, r_squared


def restricted_permanova(
    distance_matrix: np.ndarray,
    labels: Sequence[object],
    strata: Sequence[str],
    permutations: int = PERMUTATIONS,
    seed: int = PERMANOVA_SEED,
) -> tuple[float, float, float, int]:
    labels_array = np.asarray(labels, dtype=object)
    strata_array = np.asarray(strata, dtype=object)
    observed_f, observed_r2 = permanova_statistic(distance_matrix, labels_array)
    rng = np.random.RandomState(seed)
    strata_levels = sorted(set(strata_array), key=lambda value: int(value))
    stratum_indices = [np.flatnonzero(strata_array == value) for value in strata_levels]
    exceedances = 0
    for _ in range(permutations):
        permuted = labels_array.copy()
        for indices in stratum_indices:
            values = permuted[indices].copy()
            rng.shuffle(values)
            permuted[indices] = values
        permuted_f, _ = permanova_statistic(distance_matrix, permuted)
        exceedances += int(permuted_f >= observed_f - 1e-12)
    p_value = (exceedances + 1) / (permutations + 1)
    return observed_f, observed_r2, p_value, exceedances


def unrestricted_permanova(
    distance_matrix: np.ndarray,
    labels: Sequence[object],
    permutations: int = PERMUTATIONS,
    seed: int = PERMANOVA_SEED,
) -> tuple[float, float, float, int]:
    labels_array = np.asarray(labels, dtype=object)
    observed_f, observed_r2 = permanova_statistic(distance_matrix, labels_array)
    # Use the same legacy MT19937 stream as the restricted test so the
    # finite 1,999-permutation result is stable and matches the locked plan.
    rng = np.random.RandomState(seed)
    exceedances = 0
    for _ in range(permutations):
        permuted = rng.permutation(labels_array)
        permuted_f, _ = permanova_statistic(distance_matrix, permuted)
        exceedances += int(permuted_f >= observed_f - 1e-12)
    p_value = (exceedances + 1) / (permutations + 1)
    return observed_f, observed_r2, p_value, exceedances


def pcoa_positive_coordinates(distance_matrix: np.ndarray) -> np.ndarray:
    distances = np.asarray(distance_matrix, dtype=float)
    n = distances.shape[0]
    centering = np.eye(n) - np.ones((n, n)) / n
    gower = -0.5 * centering @ (distances**2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gower)
    positive = eigenvalues > 1e-12
    return eigenvectors[:, positive] * np.sqrt(eigenvalues[positive])


def distances_to_group_centroids(
    coordinates: np.ndarray, labels: Sequence[object]
) -> np.ndarray:
    labels_array = np.asarray(labels, dtype=object)
    distances = np.empty(len(labels_array), dtype=float)
    for level in sorted(set(labels_array)):
        indices = np.flatnonzero(labels_array == level)
        centroid = coordinates[indices].mean(axis=0)
        distances[indices] = np.linalg.norm(coordinates[indices] - centroid, axis=1)
    return distances


def one_way_f(values: np.ndarray, labels: Sequence[object]) -> float:
    labels_array = np.asarray(labels, dtype=object)
    values = np.asarray(values, dtype=float)
    levels = sorted(set(labels_array))
    grand_mean = values.mean()
    ss_between = sum(
        np.sum(labels_array == level)
        * (values[labels_array == level].mean() - grand_mean) ** 2
        for level in levels
    )
    ss_within = sum(
        np.sum(
            (values[labels_array == level] - values[labels_array == level].mean())
            ** 2
        )
        for level in levels
    )
    return (ss_between / (len(levels) - 1)) / (
        ss_within / (len(values) - len(levels))
    )


def permdisp_two_group_mean_distance_test(
    distance_matrix: np.ndarray,
    labels: Sequence[object],
    permutations: int = PERMUTATIONS,
    seed: int = PERMANOVA_SEED,
) -> tuple[float, float, int]:
    """Two-group PERMDISP permutation test used by the locked analysis plan.

    The Bray--Curtis matrix is embedded in principal-coordinate space. For each
    group, sample distances to that group's centroid are calculated, and the
    test statistic is the absolute difference between the two groups' mean
    centroid distances. Group centroids and distances are recomputed after
    every unrestricted label permutation.

    With two groups this directly tests the mean-dispersion contrast and is the
    precise deterministic definition underlying the accepted baseline value.
    """
    labels_array = np.asarray(labels, dtype=object)
    levels = sorted(set(labels_array))
    if len(levels) != 2:
        raise ValueError("This locked PERMDISP contrast requires exactly two groups")
    coordinates = pcoa_positive_coordinates(distance_matrix)

    def statistic(grouping: np.ndarray) -> float:
        means: list[float] = []
        for level in levels:
            indices = np.flatnonzero(grouping == level)
            centroid = coordinates[indices].mean(axis=0)
            distances = np.linalg.norm(coordinates[indices] - centroid, axis=1)
            means.append(float(distances.mean()))
        return abs(means[0] - means[1])

    observed = statistic(labels_array)
    rng = np.random.RandomState(seed)
    exceedances = 0
    for _ in range(permutations):
        permuted = rng.permutation(labels_array)
        exceedances += int(statistic(permuted) >= observed - 1e-12)
    p_value = (exceedances + 1) / (permutations + 1)
    return observed, p_value, exceedances


# ---------------------------------------------------------------------------
# Rarefaction and alpha diversity
# ---------------------------------------------------------------------------


def rarefy_one_sample(counts: np.ndarray, depth: int, seed: int) -> np.ndarray:
    counts = np.asarray(counts, dtype=np.int64)
    total = int(counts.sum())
    if total < depth:
        raise ValueError(f"Cannot rarefy depth {total} to {depth}")
    # The specified implementation resets the same RandomState for every sample.
    rng = np.random.RandomState(seed)
    expanded = np.repeat(np.arange(len(counts), dtype=np.int32), counts)
    selected = rng.choice(expanded, size=depth, replace=False)
    return np.bincount(selected, minlength=len(counts)).astype(np.int64)


def alpha_diversity_test(
    otu_counts: np.ndarray,
    sample_indices: Sequence[int],
    remove_mask: np.ndarray,
) -> tuple[float, float, float, float]:
    retain = ~np.asarray(remove_mask, dtype=bool)
    sample_indices_array = np.asarray(sample_indices, dtype=int)
    # Required deterministic order (animal, then timepoint) is constructed by
    # the caller and preserved by sample_indices_array.
    shannon: list[float] = []
    richness: list[float] = []
    for sample_index in sample_indices_array:
        rarefied = rarefy_one_sample(
            otu_counts[retain, sample_index], DEPTH_THRESHOLD, RAREFACTION_SEED
        )
        positive = rarefied[rarefied > 0].astype(float)
        proportions = positive / positive.sum()
        shannon.append(float(-np.sum(proportions * np.log(proportions))))
        richness.append(float(np.sum(rarefied > 0)))
    return np.asarray(shannon), np.asarray(richness)


# ---------------------------------------------------------------------------
# CSV/log helpers
# ---------------------------------------------------------------------------


def write_csv(path: Path, header: Sequence[str], rows: Iterable[Sequence[object]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def fmt(value: float, digits: int = 9) -> str:
    return f"{value:.{digits}f}"


class RunLog:
    def __init__(self) -> None:
        self.lines: list[str] = []

    def write(self, text: str = "") -> None:
        self.lines.append(text)
        print(text)

    def section(self, title: str) -> None:
        self.write()
        self.write(f"[{title}]")

    def save(self, path: Path) -> None:
        path.write_text("\n".join(self.lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Full workflow
# ---------------------------------------------------------------------------


def run(input_path: Path, output_dir: Path) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    log = RunLog()
    log.write("Independent downstream analysis run")
    log.write(f"Input: {input_path}")
    log.write(f"NumPy: {np.__version__}; SciPy: {scipy.__version__}")
    log.write(
        f"Seeds: PERMANOVA={PERMANOVA_SEED}, rarefaction={RAREFACTION_SEED}; "
        f"permutations={PERMUTATIONS}"
    )

    data = load_otu_workbook(input_path)
    controls = [item.index for item in data.sample_info if item.is_control]
    baseline = [item.index for item in data.sample_info if item.is_baseline]
    primary = [item.index for item in data.sample_info if item.is_primary_lymph]
    fecal = [item.index for item in data.sample_info if item.sample_type == "FB"]
    balf = [item.index for item in data.sample_info if item.sample_type == "BALF"]
    post_clp_all = [
        item.index
        for item in data.sample_info
        if item.sample_type == "lymph" and item.timepoint in {4, 8, 12, 24}
    ]
    biological = primary + fecal + balf

    log.section("INPUT AND QC")
    log.write(f"OTUs: {data.counts.shape[0]}")
    log.write(f"Samples: {data.counts.shape[1]}")
    log.write(
        f"Biological screen samples: {len(biological)} "
        f"({len(primary)} lymph + {len(fecal)} fecal + {len(balf)} BALF)"
    )
    log.write(f"Baseline lymph: {len(baseline)}; controls: {len(controls)}")
    for sample_id in ("1_4hLB", "5_4hLB"):
        sample_index = data.sample_ids.index(sample_id)
        log.write(
            f"Excluded {sample_id}: depth={data.sample_info[sample_index].depth}"
        )
    time_counts = OrderedDict(
        (timepoint, sum(data.sample_info[index].timepoint == timepoint for index in primary))
        for timepoint in (4, 8, 12, 24)
    )
    log.write(
        "Retained primary lymph by timepoint: "
        + ", ".join(f"{time}h={count}" for time, count in time_counts.items())
    )

    write_csv(
        output_dir / "S3_sample_accounting.csv",
        [
            "sample_id",
            "animal",
            "sample_type",
            "timepoint_h",
            "depth",
            "is_control",
            "primary_included",
            "excluded_reason",
        ],
        (
            (
                item.sample_id,
                item.animal,
                "LB" if item.sample_type == "lymph" else item.sample_type,
                "" if item.timepoint is None else item.timepoint,
                item.depth,
                int(item.is_control),
                int(item.is_primary_lymph),
                item.excluded_reason,
            )
            for item in data.sample_info
        ),
    )

    # Primary decontamination definitions.
    prevalence_mask, prevalence_p, prevalence_q = prevalence_screen(
        data, controls, biological
    )
    catalogue_mask = np.asarray(
        [
            bool(genus) and genus.lower() in CONTAMINANT_GENERA_LOWER
            for genus in data.genus
        ],
        dtype=bool,
    )
    union_mask = prevalence_mask | catalogue_mask
    intersection_mask = prevalence_mask & catalogue_mask
    rules = OrderedDict(
        prevalence=prevalence_mask,
        catalogue=catalogue_mask,
        union=union_mask,
        intersection=intersection_mask,
    )

    log.section("DECONTAMINATION")
    log.write(
        f"Prevalence p<0.05: {int(prevalence_mask.sum())}; "
        f"BH across all OTUs q<0.05: {int((prevalence_q < 0.05).sum())}"
    )
    log.write(
        f"Catalogue OTUs: {int(catalogue_mask.sum())}; union: {int(union_mask.sum())}; "
        f"intersection: {int(intersection_mask.sum())}"
    )

    decontamination_rows = []
    type_indices = OrderedDict(lymph=primary, fecal=fecal, BALF=balf)
    for rule_name, mask in rules.items():
        fractions = {
            sample_type: removed_read_fraction(data, mask, indices)
            for sample_type, indices in type_indices.items()
        }
        decontamination_rows.append(
            (
                rule_name,
                int(mask.sum()),
                100 * fractions["lymph"],
                100 * fractions["fecal"],
                100 * fractions["BALF"],
            )
        )
        log.write(
            f"{rule_name}: removed_OTUs={int(mask.sum())}; "
            f"lymph={100*fractions['lymph']:.4f}%; "
            f"fecal={100*fractions['fecal']:.4f}%; "
            f"BALF={100*fractions['BALF']:.4f}%"
        )

    write_csv(
        output_dir / "Table1_decontamination_rules.csv",
        ["rule", "removed_OTUs", "lymph_removed_pct", "fecal_removed_pct", "BALF_removed_pct"],
        decontamination_rows,
    )
    write_csv(
        output_dir / "S4_genus_catalogue.csv",
        ["genus", "matching_rule", "sources"],
        ((genus, "case-insensitive exact SILVA genus match", S4_SOURCES) for genus in CONTAMINANT_GENERA),
    )
    write_csv(
        output_dir / "S2_removed_OTUs.csv",
        [
            "otu_id",
            "taxonomy",
            "family",
            "genus",
            "fisher_p",
            "BH_FDR_all_8042",
            "prevalence_flag",
            "catalogue_flag",
        ],
        (
            (
                data.otu_ids[index],
                data.taxonomy[index],
                data.family[index],
                data.genus[index],
                prevalence_p[index],
                prevalence_q[index],
                int(prevalence_mask[index]),
                int(catalogue_mask[index]),
            )
            for index in np.flatnonzero(union_mask)
        ),
    )

    # Main family-level differential abundance.
    main_da = family_da(data, union_mask, primary, include_animal=True)
    ent_index = find_row(main_da, "Enterococcaceae")
    lineage_index = find_row(main_da, "norank_o__Lineage_IV")
    significant = int(np.sum(np.asarray(main_da["q"]) < 0.05))

    log.section("FAMILY-LEVEL CLR DIFFERENTIAL ABUNDANCE")
    log.write(
        f"Families tested: {len(main_da['family'])}; df=({main_da['df_num']},{main_da['df_den']}); "
        f"FDR-significant: {significant}"
    )
    log.write(
        "Enterococcaceae: "
        f"F={main_da['F'][ent_index]:.9f}, p={main_da['p'][ent_index]:.9f}, "
        f"FDR={main_da['q'][ent_index]:.9f}"
    )
    log.write(
        "norank_o__Lineage_IV: "
        f"F={main_da['F'][lineage_index]:.9f}, p={main_da['p'][lineage_index]:.9f}, "
        f"FDR={main_da['q'][lineage_index]:.9f}"
    )
    order = np.argsort(np.asarray(main_da["p"]), kind="mergesort")
    write_csv(
        output_dir / "S1_primary_families.csv",
        ["family", "F", "p_value", "BH_FDR", "df_numerator", "df_denominator"],
        (
            (
                main_da["family"][index],
                main_da["F"][index],
                main_da["p"][index],
                main_da["q"][index],
                main_da["df_num"],
                main_da["df_den"],
            )
            for index in order
        ),
    )

    # Leave-one-animal-out Enterococcaceae.
    loo_rows = []
    for animal in sorted({data.sample_info[index].animal for index in primary}, key=int):
        remaining = [index for index in primary if data.sample_info[index].animal != animal]
        result = family_da(data, union_mask, remaining, include_animal=True)
        row = find_row(result, "Enterococcaceae")
        loo_rows.append((animal, result["F"][row], result["p"][row]))
    loo_p = [float(row[2]) for row in loo_rows]
    log.write(
        f"Enterococcaceae leave-one-animal-out p range: {min(loo_p):.6f}–{max(loo_p):.6f}"
    )
    write_csv(
        output_dir / "S_LOO_enterococcaceae.csv",
        ["omitted_animal", "F", "p_value"],
        loo_rows,
    )

    # Community structure under four rules.
    primary_labels = [data.sample_info[index].timepoint for index in primary]
    primary_animals = [data.sample_info[index].animal for index in primary]
    s6_rows = []
    log.section("PERMANOVA AND ALPHA DIVERSITY BY RULE")

    # Required alpha order: sorted animal, then timepoint.
    primary_alpha_order = sorted(
        primary,
        key=lambda index: (
            int(data.sample_info[index].animal),
            int(data.sample_info[index].timepoint or 0),
        ),
    )
    alpha_groups = [data.sample_info[index].timepoint for index in primary_alpha_order]
    alpha_animals = [data.sample_info[index].animal for index in primary_alpha_order]

    for rule_name, mask in rules.items():
        relative = family_relative_abundance(data, mask, primary)
        distance = bray_curtis_matrix(relative)
        pseudo_f, r2, p_value, exceedances = restricted_permanova(
            distance, primary_labels, primary_animals
        )
        da_result = family_da(data, mask, primary, include_animal=True)
        q_significant = int(np.sum(np.asarray(da_result["q"]) < 0.05))
        shannon, richness = alpha_diversity_test(
            data.counts, primary_alpha_order, mask
        )
        _, shannon_p, _, _ = nested_f_test_matrix(
            shannon[None, :], alpha_groups, alpha_animals, include_animal=True
        )
        _, richness_p, _, _ = nested_f_test_matrix(
            richness[None, :], alpha_groups, alpha_animals, include_animal=True
        )
        s6_rows.append(
            (
                rule_name,
                len(da_result["family"]),
                q_significant,
                pseudo_f,
                r2,
                p_value,
                float(shannon_p[0]),
                float(richness_p[0]),
            )
        )
        log.write(
            f"{rule_name}: families={len(da_result['family'])}, FDR_sig={q_significant}, "
            f"PERMANOVA F={pseudo_f:.9f}, R2={r2:.9f}, p={p_value:.4f} "
            f"(exceedances={exceedances}); alpha Shannon p={shannon_p[0]:.9f}, "
            f"richness p={richness_p[0]:.9f}"
        )

    write_csv(
        output_dir / "S6_four_rule_stats.csv",
        [
            "rule",
            "families_tested",
            "FDR_significant_families",
            "PERMANOVA_pseudo_F",
            "PERMANOVA_R2",
            "PERMANOVA_p",
            "Shannon_p_decontaminated",
            "richness_p_decontaminated",
        ],
        s6_rows,
    )

    # Main-text alpha diversity: undecontaminated rarefied OTUs.
    no_removal = np.zeros(data.counts.shape[0], dtype=bool)
    raw_shannon, raw_richness = alpha_diversity_test(
        data.counts, primary_alpha_order, no_removal
    )
    raw_shannon_f, raw_shannon_p, _, _ = nested_f_test_matrix(
        raw_shannon[None, :], alpha_groups, alpha_animals, include_animal=True
    )
    raw_richness_f, raw_richness_p, _, _ = nested_f_test_matrix(
        raw_richness[None, :], alpha_groups, alpha_animals, include_animal=True
    )
    log.section("MAIN-TEXT ALPHA DIVERSITY")
    log.write(
        f"Undecontaminated rarefied Shannon: F={raw_shannon_f[0]:.9f}, p={raw_shannon_p[0]:.9f}"
    )
    log.write(
        f"Undecontaminated rarefied richness: F={raw_richness_f[0]:.9f}, p={raw_richness_p[0]:.9f}"
    )

    # Baseline vs controls (raw family relative abundance).  The baseline
    # animals are ordered numerically (1,3,5,6,9), followed by controls in the
    # source-table order. This fixes finite Monte-Carlo RNG consumption.
    baseline_ordered = sorted(
        baseline, key=lambda index: int(data.sample_info[index].animal)
    )
    baseline_control = baseline_ordered + controls
    baseline_labels = ["baseline"] * len(baseline_ordered) + ["control"] * len(controls)
    baseline_relative = family_relative_abundance(
        data, no_removal, baseline_control
    )
    baseline_distance = bray_curtis_matrix(baseline_relative)
    baseline_f, baseline_r2, baseline_p, baseline_exceed = unrestricted_permanova(
        baseline_distance, baseline_labels
    )
    dispersion_stat, dispersion_p, dispersion_exceed = (
        permdisp_two_group_mean_distance_test(baseline_distance, baseline_labels)
    )

    baseline_positions = [
        position
        for position, index in enumerate(baseline_control)
        if data.sample_info[index].is_baseline
    ]
    control_positions = [
        position
        for position, index in enumerate(baseline_control)
        if data.sample_info[index].is_control
    ]
    between_distances = [
        baseline_distance[first, second]
        for first in baseline_positions
        for second in control_positions
    ]
    within_control_distances = [
        baseline_distance[control_positions[first], control_positions[second]]
        for first in range(len(control_positions))
        for second in range(first)
    ]

    gut_mask = np.asarray(
        [bool(genus) and genus.lower() in GUT_GENERA_LOWER for genus in data.genus],
        dtype=bool,
    )
    reagent_mask = catalogue_mask

    def pooled_fraction(mask: np.ndarray, indices: Sequence[int]) -> float:
        columns = np.asarray(indices, dtype=int)
        return float(data.counts[mask][:, columns].sum() / data.counts[:, columns].sum())

    baseline_gut = pooled_fraction(gut_mask, baseline)
    baseline_reagent = pooled_fraction(reagent_mask, baseline)
    controls_gut = pooled_fraction(gut_mask, controls)
    controls_reagent = pooled_fraction(reagent_mask, controls)

    log.section("BASELINE VS CONTROLS")
    log.write(
        f"PERMANOVA: pseudo-F={baseline_f:.9f}, R2={baseline_r2:.9f}, "
        f"p={baseline_p:.4f} (exceedances={baseline_exceed}; manuscript rounds to 0.40)"
    )
    log.write(
        f"PERMDISP two-group mean-distance contrast={dispersion_stat:.9f}, "
        f"p={dispersion_p:.4f} (exceedances={dispersion_exceed}; manuscript rounds to 0.17)"
    )
    log.write(
        f"Mean Bray-Curtis baseline-to-control={np.mean(between_distances):.9f}; "
        f"within-control={np.mean(within_control_distances):.9f}"
    )
    log.write(
        f"Pooled baseline gut fraction={100*baseline_gut:.4f}%; "
        f"reagent-catalogue fraction={100*baseline_reagent:.4f}%"
    )
    write_csv(
        output_dir / "Fig3_baseline_composition.csv",
        ["group", "gut_associated_fraction", "reagent_catalogue_fraction"],
        [
            ("baseline_0h", baseline_gut, baseline_reagent),
            ("negative_controls", controls_gut, controls_reagent),
        ],
    )

    # Stepwise ablation.
    log.section("STEPWISE ABLATION")
    ablation_definitions = [
        ("0_naive", baseline + post_clp_all, False),
        ("1_remove_0h", post_clp_all, False),
        ("2_depth_qc", primary, False),
        ("3_animal_adjusted", primary, True),
    ]
    ablation_rows = []
    for step, indices, adjusted in ablation_definitions:
        result = family_da(data, union_mask, indices, include_animal=adjusted)
        ent = find_row(result, "Enterococcaceae")
        number_significant = int(np.sum(np.asarray(result["q"]) < 0.05))
        row = (
            step,
            len(indices),
            len(result["family"]),
            number_significant,
            result["F"][ent],
            result["p"][ent],
            result["q"][ent],
            result["df_num"],
            result["df_den"],
        )
        ablation_rows.append(row)
        log.write(
            f"{step}: samples={len(indices)}, families={len(result['family'])}, "
            f"FDR_sig={number_significant}, Enterococcaceae p={result['p'][ent]:.9f}"
        )
    write_csv(
        output_dir / "S5_ablation_full.csv",
        [
            "step",
            "n_samples",
            "families_tested",
            "FDR_significant_families",
            "Enterococcaceae_F",
            "Enterococcaceae_p",
            "Enterococcaceae_BH_FDR",
            "df_numerator",
            "df_denominator",
        ],
        ablation_rows,
    )

    # Sensitivity: baseline plus controls as the contaminant reference set.
    sensitivity_references = controls + baseline
    sensitivity_prevalence, sensitivity_p, sensitivity_q = prevalence_screen(
        data, sensitivity_references, biological
    )
    sensitivity_union = sensitivity_prevalence | catalogue_mask
    sensitivity_da = family_da(data, sensitivity_union, primary, include_animal=True)
    sensitivity_ent = find_row(sensitivity_da, "Enterococcaceae")
    sensitivity_relative = family_relative_abundance(
        data, sensitivity_union, primary
    )
    sensitivity_distance = bray_curtis_matrix(sensitivity_relative)
    sensitivity_f, sensitivity_r2, sensitivity_perm_p, sensitivity_exceed = restricted_permanova(
        sensitivity_distance, primary_labels, primary_animals
    )
    log.section("0H-AS-ADDITIONAL-REFERENCE SENSITIVITY")
    log.write(
        f"Prevalence flags: {int(prevalence_mask.sum())} -> {int(sensitivity_prevalence.sum())}"
    )
    log.write(
        f"Primary lymph reads removed by prevalence support: "
        f"{100*removed_read_fraction(data, prevalence_mask, primary):.4f}% -> "
        f"{100*removed_read_fraction(data, sensitivity_prevalence, primary):.4f}%"
    )
    log.write(
        f"Union removed OTUs={int(sensitivity_union.sum())}; families={len(sensitivity_da['family'])}; "
        f"FDR_sig={int(np.sum(np.asarray(sensitivity_da['q']) < .05))}; "
        f"Enterococcaceae p={sensitivity_da['p'][sensitivity_ent]:.9f}"
    )
    log.write(
        f"PERMANOVA F={sensitivity_f:.9f}, R2={sensitivity_r2:.9f}, "
        f"p={sensitivity_perm_p:.4f} (exceedances={sensitivity_exceed})"
    )

    # Comprehensive acceptance summary. Values are tested against the locked
    # task specification; none of the observed analysis results are hard-coded.
    acceptance_rows: list[tuple[str, str, str, str, str, bool]] = []

    def add_exact(section: str, item: str, observed: object, target: object) -> None:
        acceptance_rows.append(
            (section, item, str(target), str(observed), "exact", observed == target)
        )

    def add_close(
        section: str, item: str, observed: float, target: float, tolerance: float
    ) -> None:
        acceptance_rows.append(
            (
                section,
                item,
                f"{target:.9g}",
                f"{observed:.9g}",
                f"abs <= {tolerance:.3g}",
                abs(float(observed) - float(target)) <= tolerance,
            )
        )

    def add_rounded(
        section: str, item: str, observed: float, target: float, digits: int
    ) -> None:
        acceptance_rows.append(
            (
                section,
                item,
                f"round({digits})={target:.{digits}f}",
                f"{observed:.9g}",
                f"rounded to {digits} decimals",
                round(float(observed), digits) == round(float(target), digits),
            )
        )

    # Design and QC.
    add_exact("Design/QC", "OTU count", data.counts.shape[0], 8042)
    add_exact("Design/QC", "sample count", data.counts.shape[1], 61)
    add_exact("Design/QC", "biological screen samples", len(biological), 48)
    add_exact("Design/QC", "primary 4-24 h lymph samples", len(primary), 18)
    add_exact("Design/QC", "4/8/12/24 h counts", tuple(time_counts.values()), (3, 5, 5, 5))
    add_exact(
        "Design/QC",
        "1_4hLB depth",
        data.sample_info[data.sample_ids.index("1_4hLB")].depth,
        475,
    )
    add_exact(
        "Design/QC",
        "5_4hLB depth",
        data.sample_info[data.sample_ids.index("5_4hLB")].depth,
        38,
    )

    # Decontamination.
    decon_by_rule = {row[0]: row for row in decontamination_rows}
    add_exact("Decontamination", "prevalence flagged OTUs", int(prevalence_mask.sum()), 62)
    add_exact(
        "Decontamination",
        "BH survivors across all 8,042 OTUs",
        int((prevalence_q < 0.05).sum()),
        8,
    )
    add_exact("Decontamination", "catalogue OTUs", int(catalogue_mask.sum()), 109)
    add_exact("Decontamination", "union OTUs", int(union_mask.sum()), 169)
    add_exact("Decontamination", "intersection OTUs", int(intersection_mask.sum()), 2)
    add_close("Decontamination", "prevalence lymph removed %", decon_by_rule["prevalence"][2], 0.48, 0.01)
    add_close("Decontamination", "catalogue lymph removed %", decon_by_rule["catalogue"][2], 42.57, 0.01)
    add_close("Decontamination", "union lymph removed %", decon_by_rule["union"][2], 43.04, 0.01)
    add_close("Decontamination", "union BALF removed %", decon_by_rule["union"][4], 64.0, 0.5)

    # Family-level differential abundance.
    add_exact("Family DA", "families tested", len(main_da["family"]), 297)
    add_exact("Family DA", "FDR-significant families", significant, 0)
    add_close("Family DA", "Enterococcaceae F", float(main_da["F"][ent_index]), 0.49, 0.01)
    add_close("Family DA", "Enterococcaceae p", float(main_da["p"][ent_index]), 0.699, 0.001)
    add_close("Family DA", "Enterococcaceae FDR", float(main_da["q"][ent_index]), 0.932, 0.001)
    add_close("Family DA", "Lineage IV F", float(main_da["F"][lineage_index]), 14.51, 0.01)
    add_close("Family DA", "Lineage IV p", float(main_da["p"][lineage_index]), 0.0006, 0.0001)
    add_close("Family DA", "Lineage IV FDR", float(main_da["q"][lineage_index]), 0.168, 0.001)
    add_rounded("Family DA", "LOO minimum p", min(loo_p), 0.41, 2)
    add_rounded("Family DA", "LOO maximum p", max(loo_p), 0.83, 2)

    # Community structure and alpha diversity.
    union_s6 = next(row for row in s6_rows if row[0] == "union")
    add_close("PERMANOVA", "primary pseudo-F", union_s6[3], 1.21, 0.01)
    add_close("PERMANOVA", "primary R2", union_s6[4], 0.21, 0.01)
    add_close("PERMANOVA", "primary p", union_s6[5], 0.136, 0.001)
    add_close("Alpha diversity", "raw Shannon p", float(raw_shannon_p[0]), 0.58, 0.01)
    add_close("Alpha diversity", "raw richness p", float(raw_richness_p[0]), 0.15, 0.01)

    # Baseline description.
    add_close("Baseline", "PERMANOVA pseudo-F", baseline_f, 1.04, 0.01)
    add_close("Baseline", "PERMANOVA R2", baseline_r2, 0.10, 0.01)
    add_rounded("Baseline", "PERMANOVA p", baseline_p, 0.40, 2)
    add_rounded("Baseline", "PERMDISP p", dispersion_p, 0.17, 2)
    add_rounded("Baseline", "mean baseline-control Bray-Curtis", float(np.mean(between_distances)), 0.32, 2)
    add_rounded("Baseline", "mean within-control Bray-Curtis", float(np.mean(within_control_distances)), 0.50, 2)
    add_close("Baseline", "pooled gut-associated reads %", 100 * baseline_gut, 0.53, 0.01)

    # Stepwise ablation.
    add_exact(
        "Ablation",
        "FDR-significant sequence",
        tuple(int(row[3]) for row in ablation_rows),
        (22, 0, 0, 0),
    )
    add_exact(
        "Ablation",
        "families tested sequence",
        tuple(int(row[2]) for row in ablation_rows),
        (311, 297, 297, 297),
    )
    add_exact(
        "Ablation",
        "Enterococcaceae p sequence (2-dp)",
        tuple(round(float(row[5]), 2) for row in ablation_rows),
        (0.47, 0.42, 0.61, 0.70),
    )

    # Four-rule sensitivity.
    add_exact(
        "Four rules",
        "all FDR-significant family counts",
        tuple(int(row[2]) for row in s6_rows),
        (0, 0, 0, 0),
    )
    add_exact(
        "Four rules",
        "PERMANOVA p sequence (3-dp)",
        tuple(round(float(row[5]), 3) for row in s6_rows),
        (0.228, 0.137, 0.136, 0.241),
    )

    # 0-h-as-reference sensitivity.
    sensitivity_removed_pct = 100 * removed_read_fraction(
        data, sensitivity_prevalence, primary
    )
    add_exact("0h sensitivity", "prevalence flags", int(sensitivity_prevalence.sum()), 123)
    add_close("0h sensitivity", "prevalence-supported lymph removal %", sensitivity_removed_pct, 7.4, 0.1)
    add_exact("0h sensitivity", "families tested", len(sensitivity_da["family"]), 295)
    add_exact(
        "0h sensitivity",
        "FDR-significant families",
        int(np.sum(np.asarray(sensitivity_da["q"]) < 0.05)),
        0,
    )
    add_close(
        "0h sensitivity",
        "Enterococcaceae p",
        float(sensitivity_da["p"][sensitivity_ent]),
        0.697,
        0.001,
    )
    add_close("0h sensitivity", "PERMANOVA p", sensitivity_perm_p, 0.127, 0.001)

    write_csv(
        output_dir / "acceptance_checks.csv",
        ["section", "item", "target", "observed", "criterion", "status"],
        (
            (section, item, target, observed, criterion, "PASS" if passed else "FAIL")
            for section, item, target, observed, criterion, passed in acceptance_rows
        ),
    )

    log.section("ACCEPTANCE CHECKS")
    for section, item, target, observed, criterion, passed in acceptance_rows:
        log.write(
            f"{'PASS' if passed else 'FAIL'}  [{section}] {item}: "
            f"observed={observed}; target={target}; criterion={criterion}"
        )
    number_failed = sum(not row[-1] for row in acceptance_rows)
    number_passed = len(acceptance_rows) - number_failed
    log.write(f"Summary: {number_passed}/{len(acceptance_rows)} checks passed")

    report_lines = [
        "# Independent rewrite acceptance report",
        "",
        f"**Result: {number_passed}/{len(acceptance_rows)} checks passed.**",
        "",
        "The checks below compare independently calculated values with the locked task specification.",
        "",
        "| Section | Item | Target | Observed | Criterion | Status |",
        "|---|---|---:|---:|---|---|",
    ]
    report_lines.extend(
        f"| {section} | {item} | {target} | {observed} | {criterion} | {'PASS' if passed else 'FAIL'} |"
        for section, item, target, observed, criterion, passed in acceptance_rows
    )
    (output_dir / "acceptance_report.md").write_text(
        "\n".join(report_lines) + "\n", encoding="utf-8"
    )

    manifest = {
        "analysis": "independent downstream rewrite",
        "input_file": str(input_path),
        "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "seeds": {"permutation": PERMANOVA_SEED, "rarefaction": RAREFACTION_SEED},
        "permutations": PERMUTATIONS,
        "acceptance": {
            "passed": number_passed,
            "total": len(acceptance_rows),
            "failed": number_failed,
        },
    }
    (output_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    if number_failed:
        log.write("Acceptance failed; inspect acceptance_checks.csv for details.")
    log_path = output_dir / "independent_rewrite_run.log"
    log.save(log_path)
    if number_failed:
        raise RuntimeError(
            f"Acceptance failed: {number_failed}/{len(acceptance_rows)} checks did not pass"
        )
    return log_path


def parse_arguments(argv: Sequence[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("otu_workbook", type=Path, help="8,042-OTU .xlsx table")
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("independent_rewrite_out"),
        help="directory for CSV outputs and run log",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    arguments = parse_arguments(sys.argv[1:] if argv is None else argv)
    try:
        log_path = run(arguments.otu_workbook.resolve(), arguments.output_dir.resolve())
    except Exception as exc:  # clear command-line failure with traceback preserved
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
    print(f"\nRun log written to: {log_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
