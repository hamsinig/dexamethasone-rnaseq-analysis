# Transcriptomic Analysis of Dexamethasone Response in Human Airway Smooth Muscle Cells

## Project Progress and Methods Report

**Dataset:** GSE52778  
**Study:** Himes et al. (2014)  
**Analysis stage:** Data acquisition, reference preparation and read quality control  
**Status:** Pre-alignment

---

# 1. Introduction

Glucocorticoids are widely used anti-inflammatory agents in the treatment of airway diseases such as asthma. Understanding how these compounds influence gene expression in airway smooth muscle cells can provide insights into their molecular mechanisms of action.

This project investigates transcriptional changes associated with dexamethasone treatment in human airway smooth muscle cells using publicly available RNA sequencing data.

The analysis is based on the study by Himes et al. (2014), which investigated glucocorticoid-responsive genes in human airway smooth muscle cells and identified genes including CRISPLD2, DUSP1 and KLF15 as responsive to dexamethasone.

The primary objective is to develop a reproducible bulk RNA-seq analysis workflow to identify differentially expressed genes between untreated and dexamethasone-treated samples.

## 1.1 Research question

Which genes exhibit altered expression following dexamethasone treatment in human airway smooth muscle cells?

## 1.2 Objectives

1. Obtain and characterize publicly available RNA-seq data.
2. Establish a reproducible computational environment and project structure.
3. Prepare a compatible reference genome and gene annotation.
4. Assess raw sequencing quality.
5. Perform adapter trimming and quality filtering.
6. Evaluate the effectiveness of preprocessing.
7. Align the processed reads to the reference genome.
8. Generate gene-level expression counts.
9. Identify differentially expressed genes associated with dexamethasone treatment.

---

# 2. Dataset Selection and Experimental Design

## 2.1 Dataset description

The dataset selected for this analysis is GSE52778, available through the NCBI Gene Expression Omnibus (GEO).

The dataset originates from:

Himes BE et al. (2014). *RNA-Seq Transcriptome Profiling Identifies CRISPLD2 as a Glucocorticoid Responsive Gene that Modulates Cytokine Function in Airway Smooth Muscle Cells.*

The study used RNA sequencing to investigate transcriptional responses to different treatments in primary human airway smooth muscle cell lines.

The original experiment included four treatment conditions:

- Untreated control
- Albuterol
- Dexamethasone
- Albuterol combined with dexamethasone

For this project, the analysis focuses on the untreated and dexamethasone-treated conditions.

## 2.2 Sample selection

Eight sequencing runs were selected, representing four paired airway smooth muscle cell lines.

| Cell line | Untreated sample | Dexamethasone sample |
|---|---|---|
| N61311 | SRR1039508 | SRR1039509 |
| N052611 | SRR1039512 | SRR1039513 |
| N080611 | SRR1039516 | SRR1039517 |
| N061011 | SRR1039520 | SRR1039521 |

This paired experimental design allows cell-line-specific variation to be considered during downstream differential expression analysis.

The planned DESeq2 design is:

`~ cellLine + treatment`

This accounts for differences between cell lines while estimating the treatment-associated expression changes.

**Current progress:** The eight sample accessions and their treatment assignments have been documented in `samples.csv`. Raw sequencing data acquisition and preprocessing have been completed for SRR1039508 only.

---

# 3. Computational Environment and Project Organization

## 3.1 Computational environment

The analysis was configured on an Apple MacBook Air with an Apple M4 processor and 16 GB RAM.

A dedicated Conda environment named `rnaseq` was created to isolate the bioinformatics dependencies.

The environment includes:

| Tool | Purpose |
|---|---|
| SRA Toolkit | Download and convert sequencing data |
| FastQC | Raw and processed read quality assessment |
| fastp | Adapter trimming and quality filtering |
| HISAT2 | Reference-guided read alignment |
| samtools | Alignment file processing |
| featureCounts | Gene-level read quantification |
| MultiQC | Aggregation of quality-control reports |
| Python | Workflow scripting and analysis |
| pandas | Numerical data processing |
| Matplotlib | Data visualization |

The project uses Git and GitHub for version control.

GitHub repository:

https://github.com/hamsinig/dexamethasone-rnaseq-analysis

## 3.2 Project directory structure

The project was organized into separate directories for raw data, reference files, scripts and analysis outputs.

```text
dexamethasone-rnaseq-analysis/
│
├── data/
│   ├── raw/
│   └── reference/
│
├── results/
│   ├── qc/
│   ├── trimmed/
│   ├── aligned/
│   ├── counts/
│   ├── de/
│   └── figures/
│
├── scripts/
│   ├── compare_fastqc.py
│   └── summarize_fastqc.py
│
├── envs/
│
├── samples.csv
├── environment.yml
├── README.md
└── .gitignore
```

Large sequencing and reference files are stored on an external SSD to avoid exhausting internal storage. Executable scripts and project configuration files remain in the GitHub repository.

---

# 4. Raw Sequencing Data Acquisition

## 4.1 Downloading sequencing data

The raw sequencing data was obtained from the NCBI Sequence Read Archive (SRA).

The SRA Toolkit was used to download the first selected sample, SRR1039508, using `prefetch`.

The downloaded SRA archive was converted into paired-end FASTQ files using `fasterq-dump` with the `--split-files` option.

The resulting files were:

- `SRR1039508_1.fastq`
- `SRR1039508_2.fastq`

The FASTQ files were subsequently compressed using gzip.

## 4.2 Raw data characteristics

The first sample contained:

| Parameter | Value |
|---|---:|
| SRA accession | SRR1039508 |
| Sequencing layout | Paired-end |
| Number of read pairs | 22,935,521 |
| Total individual reads | 45,871,042 |
| Read length | 63 bp |
| Total bases | 2,889,875,646 |

The resulting compressed FASTQ files were retained as the original raw sequencing data and were not modified during subsequent preprocessing.

---

# 5. Reference Genome and Annotation Preparation

## 5.1 Reference genome selection

The original study used the hg19 human reference genome for read alignment.

To maintain compatibility with the original study, the Ensembl GRCh37.75 reference assembly was selected for this reanalysis.

The following reference files were downloaded:

- `Homo_sapiens.GRCh37.75.dna.primary_assembly.fa`
- `Homo_sapiens.GRCh37.75.gtf`

The FASTA file contains the reference nucleotide sequences, while the GTF file provides gene, transcript and exon annotation information.

## 5.2 Reference annotation processing

The GTF annotation was examined to understand its structure and verify chromosome naming compatibility with the reference FASTA.

The annotation includes information such as:

- Gene identifiers
- Transcript identifiers
- Gene names
- Genomic coordinates
- Strand orientation
- Exon boundaries

HISAT2 helper files were generated from the GTF annotation:

- `splice_sites.txt`
- `exons.txt`

These files contain splice-junction and exon information that can assist splice-aware alignment.

## 5.3 HISAT2 index construction

A HISAT2 index was constructed using the GRCh37 reference FASTA.

The resulting index was tested using a small artificial FASTQ file to verify that HISAT2 could successfully load and process reads against the reference.

The standard genome index was retained, with splice-site information intended to be supplied separately during alignment.

**Current status:** Reference genome, annotation, helper files and HISAT2 index have been prepared and tested.

---

# 6. Raw Read Quality Assessment

## 6.1 FastQC analysis

FastQC was used to assess the quality of the raw paired-end reads before preprocessing.

Both R1 and R2 were examined independently.

The following FastQC modules were evaluated:

1. Per base sequence quality
2. Per tile sequence quality
3. Per sequence quality scores
4. Per base sequence content
5. Per sequence GC content
6. Per base N content
7. Sequence length distribution
8. Sequence duplication levels
9. Overrepresented sequences
10. Adapter content
11. Basic statistics

## 6.2 Initial quality observations

The raw sequencing data showed generally high per-base quality.

Mean quality scores began around Q39–40 and declined towards approximately Q34 near the ends of the reads.

Other observations included:

- Localized per-tile quality variation.
- Initial nucleotide composition fluctuations.
- Low adapter content according to the FastQC Adapter Content module.
- Small amounts of overrepresented adapter or primer-associated sequences.
- Elevated sequence duplication estimates.
- A localized increase in ambiguous nucleotide content in R2 around position 9.

The initial data was considered suitable for conservative preprocessing rather than aggressive trimming.

---

# 7. Read Preprocessing Using fastp

## 7.1 Preprocessing strategy

The raw paired-end reads were processed using fastp v1.4.0.

The following filtering parameters were applied:

| Parameter | Setting |
|---|---|
| Adapter detection | Automatic paired-end detection |
| Qualified base threshold | Q20 |
| Maximum unqualified-base percentage | 20% |
| Maximum ambiguous bases per read | 5 |
| Minimum retained read length | 30 bp |
| Threads | 4 |

No sliding-window trimming, arbitrary removal of the first ten bases, or automatic duplicate removal was performed.

The conservative approach was selected because the raw reads already exhibited relatively high quality.

## 7.2 Preprocessing results

The initial dataset contained 22,935,521 paired-end fragments.

After preprocessing, 21,977,189 read pairs were retained.

This corresponds to a read retention rate of 95.82%.

A total of 1,916,664 individual reads failed filtering.

### Table 1. Sequencing statistics before and after preprocessing

| Parameter | Before | After |
|---|---:|---:|
| Total individual reads | 45,871,042 | 43,954,378 |
| Read pairs | 22,935,521 | 21,977,189 |
| Total bases | 2,889,875,646 | 2,766,453,310 |
| Mean R1 read length | 63 bp | 62 bp |
| Mean R2 read length | 63 bp | 62 bp |
| Q20 bases | 98.08% | 99.23% |
| Q30 bases | 96.05% | 97.59% |
| GC content | 49.98% | 49.86% |
| Read retention | — | 95.82% |

### Table 2. Preprocessing summary

| Parameter | Result |
|---|---:|
| Individual reads passed filtering | 43,954,378 |
| Individual reads failed filtering | 1,916,664 |
| Reads with adapter sequences trimmed | 260,067 |
| Adapter-associated bases trimmed | 4,488,299 |
| fastp duplication rate | 7.07% |
| Estimated insert-size peak | 95 bp |

The increase in Q20 and Q30 proportions indicates an improvement in the quality distribution of the retained reads.

GC content remained relatively stable, suggesting that preprocessing did not substantially alter the overall nucleotide composition.

---

# 8. Post-Preprocessing Quality Assessment

## 8.1 FastQC comparison

FastQC was run again on the trimmed paired-end reads.

To evaluate the effect of preprocessing, the numerical data from the raw and trimmed FastQC ZIP reports were extracted using Python.

A custom script, `compare_fastqc.py`, was developed to generate side-by-side visualizations of the FastQC modules.

The comparison was generated separately for R1 and R2.

A multipage PDF containing the visual comparisons was saved as:

`SRR1039508_full_fastqc_comparison.pdf`

## 8.2 Numerical comparison

### Table 3. Per-base sequencing quality

| Metric | R1 Raw | R1 Trimmed | R2 Raw | R2 Trimmed |
|---|---:|---:|---:|---:|
| Mean per-base quality (Phred) | 38.24 | 38.68 | 38.04 | 38.61 |
| Mean per-base median quality (Phred) | 39.65 | 39.68 | 39.73 | 39.76 |
| Final-base mean quality (Phred) | 33.87 | 34.63 | 33.87 | 34.77 |
| Minimum per-base mean quality (Phred) | 33.87 | 34.63 | 33.87 | 34.42 |

Both reads showed modest improvements in mean per-base quality following preprocessing.

### Table 4. GC content and ambiguous bases

| Metric | R1 Raw | R1 Trimmed | R2 Raw | R2 Trimmed |
|---|---:|---:|---:|---:|
| Weighted mean GC (%) | 50.02 | 49.91 | 49.92 | 49.82 |
| Maximum per-base N content (%) | 0.327 | 0.025 | 10.684 | 10.611 |
| Position of maximum N content | 22 | 20 | 9 | 9 |

R1 showed a substantial reduction in maximum ambiguous-base content.

However, R2 continued to exhibit elevated N content around position 9, indicating that the issue was not resolved by the selected preprocessing parameters.

### Table 5. Maximum adapter content

| Adapter type | R1 Raw | R1 Trimmed | R2 Raw | R2 Trimmed |
|---|---:|---:|---:|---:|
| Illumina Universal Adapter | 0.2046 | 0.0018 | 0.2082 | 0.0422 |
| Illumina Small RNA 3′ Adapter | 0.0003 | 0.0002 | 0.0003 | 0.0002 |
| Illumina Small RNA 5′ Adapter | 0.0000 | 0.0000 | 0.0000 | 0.0000 |
| Nextera Transposase Sequence | 0.0002 | 0.0002 | 0.0002 | 0.0002 |
| PolyA | 0.4792 | 0.2410 | 0.6066 | 0.3324 |
| PolyG | 0.0208 | 0.0018 | 0.0374 | 0.0017 |

Adapter content decreased for several adapter types, particularly the Illumina Universal Adapter, PolyA and PolyG sequences.

### Table 6. Sequence duplication estimates

| Metric | R1 Raw | R1 Trimmed | R2 Raw | R2 Trimmed |
|---|---:|---:|---:|---:|
| FastQC deduplicated percentage | 49.70% | 49.39% | 52.06% | 51.51% |

FastQC continued to flag sequence duplication. These estimates should not automatically be interpreted as PCR duplication because highly expressed RNA transcripts naturally generate repeated sequences.

No automatic deduplication was performed.

## 8.3 FastQC module status comparison

| Module | R1 Raw | R1 Trimmed | R2 Raw | R2 Trimmed |
|---|---|---|---|---|
| Adapter Content | PASS | PASS | PASS | PASS |
| Basic Statistics | PASS | PASS | PASS | PASS |
| Overrepresented Sequences | WARN | WARN | WARN | WARN |
| Per Base N Content | PASS | PASS | WARN | WARN |
| Per Base Sequence Content | PASS | PASS | PASS | PASS |
| Per Base Sequence Quality | PASS | PASS | PASS | PASS |
| Per Sequence GC Content | PASS | PASS | PASS | PASS |
| Per Sequence Quality Scores | PASS | PASS | PASS | PASS |
| Per Tile Sequence Quality | WARN | WARN | FAIL | FAIL |
| Sequence Duplication Levels | FAIL | FAIL | WARN | WARN |
| Sequence Length Distribution | PASS | WARN | PASS | WARN |

The overall FastQC status profile remained broadly similar following preprocessing.

The sequence length distribution changed from PASS to WARN because filtering introduced variable read lengths. This does not necessarily indicate poor-quality data.

Per-tile sequence quality issues persisted, particularly in R2.

## 8.4 Interpretation

Overall, preprocessing improved the quality distribution of the retained reads while preserving approximately 95.82% of the original sequencing data.

The proportion of Q20 and Q30 bases increased, adapter-associated signals decreased, and GC composition remained stable.

However, elevated N content and per-tile quality abnormalities persisted in R2. These issues will be documented and considered during downstream alignment assessment.

Sequence duplication estimates remained elevated, but no duplicate removal was performed because duplication in bulk RNA-seq can reflect biological transcript abundance.

---

# 9. Current Project Status

At this stage, the following tasks have been completed:

| Task | Status |
|---|---|
| GitHub repository setup | Completed |
| Conda environment setup | Completed |
| Dataset selection | Completed |
| Sample metadata preparation | Completed |
| Download of SRR1039508 | Completed |
| FASTQ conversion and compression | Completed |
| GRCh37 reference genome download | Completed |
| GTF annotation preparation | Completed |
| HISAT2 index construction and test | Completed |
| Raw FastQC assessment for SRR1039508 | Completed |
| fastp preprocessing for SRR1039508 | Completed |
| Post-trimming FastQC assessment | Completed |
| Numerical QC comparison | Completed |
| Sample-specific QC report | Completed |

The remaining seven samples have not yet been processed through this workflow.

## 9.1 Next steps

The planned continuation of the project is:

1. Acquire and preprocess the remaining seven sequencing samples.
2. Generate quality-control summaries for all eight samples.
3. Perform splice-aware alignment using HISAT2.
4. Sort and index alignment files using samtools.
5. Quantify gene-level read counts using featureCounts.
6. Perform differential expression analysis using DESeq2.
7. Investigate significant genes and biological pathways.

No alignment or differential expression results are reported at this stage.

---

# 10. Reproducibility and Project Files

The project is maintained in a GitHub repository:

https://github.com/hamsinig/dexamethasone-rnaseq-analysis

Important scripts developed so far:

- `scripts/compare_fastqc.py` — generates side-by-side raw and trimmed FastQC visualizations.
- `scripts/summarize_fastqc.py` — extracts numerical QC measurements and generates summary reports.

Important sample-specific output files:

- `SRR1039508_full_fastqc_comparison.pdf`
- `SRR1039508_qc_comparison.csv`
- `SRR1039508_fastp_summary.csv`
- `SRR1039508_qc_report.md`

Raw sequencing data, processed FASTQ files and reference files are stored separately on the external SSD.

The analysis environment is documented in `environment.yml` to support reproducibility.

---

# 11. Conclusion

A reproducible computational workflow has been established for investigating dexamethasone-associated transcriptional changes in human airway smooth muscle cells.

The dataset and experimental design were characterized, the GRCh37 reference genome and annotation were prepared, and an initial paired-end sequencing sample was successfully processed through quality assessment and preprocessing.

The preprocessing strategy retained most of the sequencing reads while improving overall quality metrics and reducing several adapter-associated signals. Some localized quality issues persisted, particularly in the second read, and will be considered during alignment assessment.

The project is currently at the pre-alignment stage. Processing the remaining samples, generating gene-level counts and conducting differential expression analysis will provide the basis for evaluating the transcriptional response to dexamethasone.
