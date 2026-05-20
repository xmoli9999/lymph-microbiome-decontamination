# Rigorous decontamination reveals a reproducible Enterococcaceae signal in mesenteric lymph during experimental sepsis

## Overview

This repository contains the analysis code and processed metadata for the manuscript:

> **Rigorous decontamination reveals a reproducible Enterococcaceae signal in mesenteric lymph during experimental sepsis**
>
> Wang ZX, Jin LZ, Gao ZG, Dong HS, Yuan CL, Li Y
>
> *Scientific Reports* (2026)

## Data Availability

- **Raw sequencing data**: NCBI Sequence Read Archive (SRA) under BioProject [PRJNA1464423](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1464423)
- **Processed OTU tables**: Zenodo (DOI: [10.5281/zenodo.19942976](https://doi.org/10.5281/zenodo.19942976))

## Repository Structure

```
├── README.md
├── LICENSE
├── requirements.txt
├── scripts/
│   ├── 01_decontamination.py           # Prevalence-based decontamination + kitome filtering
│   ├── 02_alpha_beta_diversity.py      # Alpha/beta diversity (pre/post decontamination)
│   ├── 03_deseq2_analysis.py           # DESeq2 differential abundance (pyDESeq2)
│   ├── 04_clr_compositional.py         # CLR-based compositional analysis
│   ├── 05_cross_method_concordance.py  # Dual-method concordance assessment
│   ├── 06_kitome_metacomparison.py     # Kitome catalog cross-referencing
│   └── 07_figures.py                   # Figure generation
└── data/
    └── sample_metadata.csv             # Sample metadata and group assignments
```

## Analysis Pipeline

1. **Decontamination** (`01_decontamination.py`): Two-pronged protocol combining Fisher's exact test (prevalence in 6 negative controls vs. 60 biological samples) with cross-referencing against published kitome catalogs (Salter et al. 2014; Eisenhofer et al. 2019; Karstens et al. 2019; de Goffau et al. 2019).

2. **Diversity Analysis** (`02_alpha_beta_diversity.py`): Shannon, Simpson, and Chao1 indices (rarefied to 5,000 reads); Bray–Curtis dissimilarity with PERMANOVA (999 permutations).

3. **DESeq2 Analysis** (`03_deseq2_analysis.py`): Negative binomial GLM via pyDESeq2 v0.4. Pairwise comparisons (0 h vs. 4, 8, 12, 24 h) at genus and family levels with Benjamini–Hochberg FDR correction.

4. **CLR Compositional Analysis** (`04_clr_compositional.py`): CLR transformation (pseudocount 0.5), sample-specific bias correction, Welch's t-test with FDR correction.

5. **Cross-Method Concordance** (`05_cross_method_concordance.py`): A taxon–comparison pair is concordant if nominal P < 0.05 in both DESeq2 and CLR with directional agreement.

6. **Kitome Meta-comparison** (`06_kitome_metacomparison.py`): Comparison of identified contaminant genera against combined published kitome catalogs.

7. **Figure Generation** (`07_figures.py`): All manuscript figures (Figures 1–5).

## Dependencies

```
Python 3.10+
pandas >= 1.5
numpy >= 1.23
scipy >= 1.9
scikit-bio >= 0.5.8
pydeseq2 >= 0.4
matplotlib >= 3.6
seaborn >= 0.12
```

## Usage

```bash
pip install -r requirements.txt

python scripts/01_decontamination.py
python scripts/02_alpha_beta_diversity.py
python scripts/03_deseq2_analysis.py
python scripts/04_clr_compositional.py
python scripts/05_cross_method_concordance.py
python scripts/06_kitome_metacomparison.py
python scripts/07_figures.py
```

## Ethics

All animal procedures were approved by the Institutional Animal Care and Use Committee of Shanghai Fourth People's Hospital (Approval No. A2015035).

## License

MIT License

## Contact

Yan Li (criticalcare@163.com) — Shanghai Fourth People's Hospital, Tongji University
