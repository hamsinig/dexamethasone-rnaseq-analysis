
#!/bin/bash

# Project directory on Mac
export PROJECT_DIR="$HOME/Desktop/projects/compbio/bioskillslab/dexamethasone-rnaseq-analysis"

# Large data directory on external SSD
export DATA_DIR="/Volumes/ghamsini/Hamsini_Documents/data/dexamethasone-rnaseq"

# Input data
export RAW_DIR="$DATA_DIR/raw"
export FASTQ_DIR="$DATA_DIR/raw/fastq"
export TRIMMED_DIR="$DATA_DIR/trimmed"

# Reference files
export REFERENCE_DIR="$DATA_DIR/reference"
export HISAT2_INDEX="$REFERENCE_DIR/hisat2_index/grch37"
export GTF="$REFERENCE_DIR/Homo_sapiens.GRCh37.75.gtf"
export SPLICE_SITES="$REFERENCE_DIR/splice_sites.txt"

# Large computational outputs on SSD
export ALIGNED_DIR="$DATA_DIR/aligned"
export COUNTS_DIR="$DATA_DIR/counts"
export TMP_DIR="$DATA_DIR/tmp"
export LOG_DIR="$DATA_DIR/logs"

# Small outputs and documentation on Mac
export REPORT_DIR="$PROJECT_DIR/reports"
export RESULTS_DIR="$PROJECT_DIR/results"
export FIGURE_DIR="$RESULTS_DIR/figures"
export TABLE_DIR="$RESULTS_DIR/tables"

