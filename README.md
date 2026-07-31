# Mesenteric lymph 16S — decontamination reanalysis (rat CLP sepsis)

Reproducible analysis for the manuscript:
**"Baseline mesenteric lymph is indistinguishable from reagent blanks and no robust
temporal 16S signal survives design-appropriate reanalysis in experimental sepsis."**

This repository reanalyses 16S rRNA (V4) profiles from a rat cecal ligation and
puncture (CLP) sepsis model, with correct handling of contamination, sequencing
depth, batch structure and the repeated-measures design. Running one script
regenerates every number, table and figure in the paper.

## What the analysis shows

- Read-level "contamination fraction" is **rule-dependent**: the six negative
  controls support removing only ~0.5% of lymph reads, whereas a genus blocklist
  removes ~40%.
- **Baseline (0 h) mesenteric lymph is near-sterile and compositionally
  indistinguishable from the negative controls** (gut genera ~0.5%; Bray–Curtis
  distance to controls ~0.32) — the clearest biological observation in the data.
- In the design-appropriate 4–24 h window (5 animals, repeated measures),
  **no bacterial family is FDR-significant** (0/297) and there is no alpha- or
  beta-diversity time effect. Enterococcaceae, emphasised in earlier work, is not
  significant (F = 0.49, p = 0.70); its apparent 4 h peak was driven by a single
  38-read sample removed at quality control.
- The null result is **robust to decontamination strategy** (sensitivity analysis
  using the near-blank 0 h samples as extra controls).

## Data

| File | Description |
|------|-------------|
| `data/otu_taxon.xls` | Majorbio OTU count table (8,042 OTUs × 66 samples) with SILVA-138 taxonomy in the final `taxonomy` column. Tab-separated despite the `.xls` extension. |
| `data/sample_metadata.csv` | Per-sample annotation: `sample_id, sample_type, group, animal_id, timepoint_h, batch, notes`. |

Raw sequencing reads (FASTQ) are in the NCBI Sequence Read Archive under
**BioProject PRJNA1464423**.

### Design notes (read before reusing the data)

- **Mesenteric lymph** was collected longitudinally from **5 animals** (IDs
  1, 3, 5, 6, 9) at 0, 4, 8, 12, 24 h — i.e. a repeated-measures design, **not**
  25 independent animals. Analyses must treat `animal_id` as a block / random
  effect.
- The **first-batch 0 h lymph** (`*_0hLB`) was contaminated (unflushed cannula)
  and re-collected as **`LB_1`–`LB_5`** (LB_1→animal 1, LB_2→3, LB_3→5, LB_4→6,
  LB_5→9), each with a paired extraction blank (`con_1`–`con_5`). Because `LB_1–5`
  are a **separate extraction batch**, 0 h is confounded with reagent batch and is
  used only descriptively; the primary analysis is the within-batch 4→24 h window.
  The original `*_0hLB` columns are retained in the table but excluded from the
  primary analysis (see `notes` in the metadata).
- **Faecal (`*_FB`)** and **BALF (`*_FP`)** samples span three condition groups
  (normal control 21/23/24/25/26; sepsis-no-drainage 11/13/14/15/16;
  sepsis-drainage 1/3/5/6/9).

## Reproduce

```bash
pip install -r requirements.txt
python reanalysis_pipeline.py
```

Outputs are written to `reanalysis_out/`:

- `decontamination_rules.csv` — OTUs removed and read% by specimen type, for four rules
- `differential_4_24h.csv` — family-level differential abundance (CLR + animal-adjusted F-test, BH-FDR)
- `fig1_decontam_rules.png` — decontamination is rule-dependent
- `fig2_enterococcaceae.png` — no Enterococcaceae time trend (4–24 h)
- `fig3_pcoa_specimen.png` — 0 h baseline clusters with the negative controls
- `fig4_baseline_blank.png` — baseline lymph is blank-like (almost no gut taxa)

The console also prints the beta-diversity PERMANOVA, alpha-diversity tests and the
decontamination-strategy sensitivity analysis.

## Methods (as implemented in `reanalysis_pipeline.py`)

1. **QC** — prespecified read-depth floor of 5,000 reads (drops `5_4hLB`=38 and
   `1_4hLB`=475); rarefaction to 5,000 for diversity.
2. **Decontamination** — control-prevalence (one-sided Fisher exact, 6 controls vs
   biological samples) ∪ genus-catalogue (reagent contaminant genera); reported
   alongside prevalence-only, catalogue-only and intersection.
3. **Differential abundance** — family-level CLR, animal-adjusted nested F-test for
   the timepoint effect, Benjamini–Hochberg FDR.
4. **Beta diversity** — Bray–Curtis PERMANOVA, permutations restricted within animal.
5. **Sensitivity** — repeat using the near-blank 0 h samples as additional controls.

## Dependencies

Python 3.10+, and only `pandas`, `numpy`, `scipy` (≥1.11) and `matplotlib`
(see `requirements.txt`). No `statsmodels`/`scikit-bio` required.

## Ethics

All animal procedures were approved by the Animal Ethics Committee of Shanghai Jiao
Tong University (approval no. **A2016091**).

## License

Code: MIT (see `LICENSE`). Data files: CC BY 4.0.
