# Design-aware reanalysis of mesenteric-lymph 16S data reveals instability of apparent temporal associations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21756011.svg)](https://doi.org/10.5281/zenodo.21756011)

Code and data for:

> **Design-aware reanalysis of mesenteric-lymph 16S data reveals instability of apparent temporal associations.**

This repository contains the **author-written, independent implementation** of the downstream
analysis. `independent_rewrite_analysis.py` was written from the statistical specification only
(it reads the OTU-with-taxonomy table with the Python standard library and implements the QC,
decontamination, family-level CLR models, PERMANOVA, alpha diversity, ablation and sensitivity
analyses with NumPy/SciPy); it does not import or reuse any prior script. `render_figures.py`
draws Figures 2–5 from the analysis outputs. The upstream steps
(Trimmomatic → FLASH → UCHIME → UPARSE → SILVA 138) were performed by the sequencing provider
and are not re-run here; the analysis begins from the 8,042-OTU table.

## Reproduce everything

```bash
pip install -r requirements.txt
python run_all.py data/otu_table2_with_taxonomy.xlsx --output-dir outputs
```

This regenerates every reported number, Supplementary Tables S1–S6 and Figures 2–5 (PNG + PDF).
The run is deterministic (PERMANOVA permutation seed = 1, rarefaction seed = 42, 1,999
permutations); **49/49 acceptance checks pass** against the pre-specified targets and two reruns
are byte-identical (`outputs/acceptance_report.md`, `outputs/determinism_check.txt`).

## Files

- `independent_rewrite_analysis.py` — author-written analysis (numbers, tables S1–S6, figure-source CSVs).
- `render_figures.py` — renders Figures 2–5 from the analysis outputs.
- `run_all.py` — runs the analysis and the figures end to end.
- `run_analysis.sh` — convenience wrapper.
- `requirements.txt` — pinned Python dependencies.
- `data/otu_table2_with_taxonomy.xlsx` — OTU-with-taxonomy table (8,042 OTUs, 61 samples).
- `SHA256SUMS.txt` — checksums of the delivered files.

## Data availability

Raw 16S rRNA reads: NCBI SRA, BioProject **PRJNA1464423** (public).
The complete, permanently archived bundle (this code, the input table, all outputs, Figures 2–5,
the run log, acceptance and determinism checks) is on Zenodo, **v2.0.0**, DOI
**10.5281/zenodo.21756011** (concept DOI; always resolves to the latest version).

## Citation / version

Please cite the paper together with the Zenodo DOI (concept DOI 10.5281/zenodo.21756011, or the
v2.0.0 version DOI 10.5281/zenodo.21802368) and, if needed, the exact release tag / commit hash.
