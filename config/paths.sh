
#!/bin/bash

# ==========================================
# PROJECT DIRECTORY — MAC
# ==========================================

export PROJECT_DIR="$HOME/Desktop/projects/compbio/bioskillslab/dexamethasone-rnaseq-analysis"

# ==========================================
# LARGE DATA DIRECTORY — EXTERNAL SSD
# ==========================================

export DATA_DIR="/Volumes/ghamsini/Hamsini_Documents/data/dexamethasone-rnaseq"

# Input sequencing data
export RAW_DIR="$DATA_DIR/raw"
export FASTQ_DIR="$RAW_DIR/fastq"
export TRIMMED_DIR="$DATA_DIR/trimmed"

# Reference genome files
export REFERENCE_DIR="$DATA_DIR/reference"
export HISAT2_INDEX="$REFERENCE_DIR/hisat2_index/grch37"
export GTF="$REFERENCE_DIR/Homo_sapiens.GRCh37.75.gtf"

# ==========================================
# SMALL REFERENCE FILES — MAC
# ==========================================

export SPLICE_SITES="$PROJECT_DIR/data/reference/splice_sites.txt"
export EXONS="$PROJECT_DIR/data/reference/exons.txt"

# ==========================================
# LARGE COMPUTATIONAL FILES — EXTERNAL SSD
# ==========================================

export ALIGNED_DIR="$DATA_DIR/aligned"
export TMP_DIR="$DATA_DIR/tmp"

# ==========================================
# PROJECT OUTPUTS — MAC
# ==========================================

export REPORT_DIR="$PROJECT_DIR/reports"
export RESULTS_DIR="$PROJECT_DIR/results"

export FIGURE_DIR="$RESULTS_DIR/figures"
export TABLE_DIR="$RESULTS_DIR/tables"
export COUNTS_DIR="$RESULTS_DIR/counts"
export LOG_DIR="$PROJECT_DIR/logs"

