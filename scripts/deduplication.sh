#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/sample.sh
source ../configs/reference.sh

echo "Running Deduplication..."
log_message "Running Deduplication..."

# Mark duplicates
gatk MarkDuplicatesSpark \
	-I $ALIGN_DIR/${SAMPLE_NAME}.sam \
	-O $ALIGN_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam

log_message "Marked duplicates for sample: $SAMPLE_NAME"
