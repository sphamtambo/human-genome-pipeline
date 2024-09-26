#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/reference.sh
source ../configs/sample.sh

echo "Collecting Metrics..."
log_message "Collecting Metrics..."

# Collect Alignment and Insert Size Metrics
gatk CollectAlignmentSummaryMetrics \
	R=$REFERENCE_GENOME \
	I=$ALIGN_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam \
	O=$METRICS_DIR/${SAMPLE_NAME}_alignment_metrics.txt

log_message "Collected Alignment Metrics for sample: $SAMPLE_NAME"

gatk CollectInsertSizeMetrics \
	R=$REFERENCE_GENOME \
	I=$ALIGN_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam \
	O=$METRICS_DIR/${SAMPLE_NAME}_insert_size_metrics.txt \
	H=$METRICS_DIR/${SAMPLE_NAME}_insert_size_histogram.pdf

log_message "Collected Insert Size Metrics for sample: $SAMPLE_NAME"
