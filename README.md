# Dexamethasone-Induced Differential Gene Expression Analysis

## Overview

This project investigates transcriptional changes in human airway smooth muscle cells following treatment with dexamethasone, a synthetic glucocorticoid commonly used in the treatment of inflammatory conditions such as asthma.

The objective is to identify differentially expressed genes between untreated and dexamethasone-treated cells using an RNA-sequencing (RNA-seq) analysis pipeline.

## Biological Question

Which genes respond to dexamethasone treatment in human airway smooth muscle cells?

Understanding these transcriptional responses can provide insights into glucocorticoid-mediated regulation of gene expression and help identify potential therapeutic targets.

## Dataset

- **Study:** Himes et al. (2014)
- **GEO Accession:** GSE52778
- **Publication:** RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Mesenchymal Cell Differentiation
- **Journal:** PLOS ONE
- **DOI:** 10.1371/journal.pone.0099625
- **Experimental groups:** 4 untreated and 4 dexamethasone-treated samples
- **Cell type:** Human airway smooth muscle cells
- **Sequencing platform:** Illumina HiSeq 2000

Dataset: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52778

## Analysis Workflow

The project follows a reproducible bulk RNA-seq workflow:

1. **Quality Control:** Assess raw sequencing reads using FastQC and MultiQC.
2. **Read Trimming:** Remove adapters and low-quality bases using fastp.
3. **Read Alignment:** Align reads to the human reference genome using STAR.
4. **Gene Quantification:** Generate gene-level read counts using featureCounts.
5. **Differential Expression Analysis:** Identify differentially expressed genes using DESeq2.
6. **Visualization:** Generate volcano plots, heatmaps and expression plots using ggplot2.

## Tools and Technologies

- Linux/Unix command-line tools
- Conda
- SRA Toolkit
- FastQC
- fastp
- STAR
- SAMtools
- Subread / featureCounts
- R
- DESeq2
- ggplot2
- MultiQC

## Expected Biological Findings

The original publication reported glucocorticoid-responsive genes including:

- CRISPLD2
- DUSP1
- KLF15

These genes will be examined as part of the differential expression analysis.

## Repository Structure

```text
dexamethasone-rnaseq-analysis/
├── data/
│   ├── raw/
│   └── reference/
├── results/
│   ├── qc/
│   ├── trimmed/
│   ├── aligned/
│   ├── counts/
│   ├── de/
│   └── figures/
├── scripts/
├── envs/
├── README.md
├── .gitignore
└── environment.yml
```

## Reproducibility

The analysis scripts, software environment specifications and workflow documentation will be maintained in this repository.

Raw sequencing data and large intermediate files will not be committed to GitHub. They can be retrieved from the public GEO/SRA archives.

## References

Himes, B. E., et al. (2014). RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Mesenchymal Cell Differentiation. PLOS ONE.

https://doi.org/10.1371/journal.pone.0099625
