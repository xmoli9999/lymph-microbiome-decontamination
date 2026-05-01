# Decontamination Reveals a Convergent Microbiome–Functional Signature of Gut–Lymph–Lung Translocation in Sepsis

## Overview

This repository contains the analysis code for the manuscript:

> **Decontamination Reveals a Convergent Microbiome–Functional Signature of Gut–Lymph–Lung Translocation in Sepsis**
>
> Submitted to *Microbiome*

## Data Availability

Raw 16S rRNA gene sequencing data have been deposited in the NCBI Sequence Read Archive (SRA) under BioProject accession **[PRJNA XXXXX]**.

## Repository Structure

```
├── README.md
├── scripts/
│   ├── 01_decontamination.py        # Prevalence-based decontamination + kitome filtering
│   ├── 02_alpha_beta_diversity.py   # Alpha/beta diversity analysis (pre/post decontamination)
│   ├── 03_deseq2_analysis.py        # DESeq2 differential abundance (pyDESeq2)
│   ├── 04_clr_compositional.py      # CLR-based compositional analysis
│   ├── 05_cross_method_concordance.py  # Dual-method concordance assessment
│   ├── 06_kitome_metacomparison.py  # Kitome catalog cross-referencing
│   └── 07_figures.py                # Figure generation
├── data/
│   ├── kitome_catalog.csv           # Combined 78-genus kitome reference catalog
│   └── sample_metadata.csv          # Sample metadata and group assignments
└── requirements.txt
```

## Dependencies

```
Python 3.10+
pandas >= 1.5
numpy >= 1.23
scipy >= 1.9
scikit-bio >= 0.5.8
pyDESeq2 >= 0.4
matplotlib >= 3.6
seaborn >= 0.12
```

## Analysis Pipeline

1. **Decontamination** (`01_decontamination.py`): Two-pronged protocol combining Fisher's exact test (prevalence in 6 negative controls vs. 60 biological samples) with cross-referencing against published kitome catalogs [Salter et al. 2014; Eisenhofer et al. 2019; Karstens et al. 2019; de Goffau et al. 2019].

2. **Diversity Analysis** (`02_alpha_beta_diversity.py`): Shannon, Simpson, and Chao1 indices (rarefied to 5,000 reads); Bray–Curtis dissimilarity with PERMANOVA (999 permutations).

3. **DESeq2 Analysis** (`03_deseq2_analysis.py`): Negative binomial GLM via pyDESeq2 v0.4. Pairwise comparisons (0 h vs. 4, 8, 12, 24 h) at genus and family levels. Benjamini–Hochberg FDR correction.

4. **CLR Compositional Analysis** (`04_clr_compositional.py`): CLR transformation (pseudocount 0.5), sample-specific bias correction, Welch's t-test with FDR correction.

5. **Cross-Method Concordance** (`05_cross_method_concordance.py`): Taxon–comparison pairs concordant if nominal P < 0.05 in both DESeq2 and CLR with directional agreement.

6. **Kitome Meta-comparison** (`06_kitome_metacomparison.py`): Comparison of identified contaminant genera against combined published kitome catalogs.

## Usage

```bash
# Install dependencies
pip install -r requirements.txt

# Run the full pipeline
python scripts/01_decontamination.py
python scripts/02_alpha_beta_diversity.py
python scripts/03_deseq2_analysis.py
python scripts/04_clr_compositional.py
python scripts/05_cross_method_concordance.py
python scripts/06_kitome_metacomparison.py
python scripts/07_figures.py
```

## Ethics

All animal procedures were approved by the Institutional Animal Care and Use Committee of Shanghai (Approval No. A2015035).

## License

This project is licensed under the MIT License.

## Contact

[Corresponding Author] - [corresponding.author@institution.edu]
