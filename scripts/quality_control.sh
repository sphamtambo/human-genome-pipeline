#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/sample.sh

echo "Running FastQC..."
log_message "Running FastQC..."

fastqc -t $THREADS -o $QC_DIR $INPUT_DIR/$FORWARD_READS $INPUT_DIR/$REVERSE_READS

log_message "Completed FastQC analysis."
