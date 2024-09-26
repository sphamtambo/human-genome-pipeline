#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/sample.sh
source ../configs/reference.sh

echo "Running BQSR..."
log_message "Running BQSR..."

# Base Quality Score Recalibration (BQSR)
gatk BaseRecalibrator \
	-I $ALIGN_DIR/${SAMPLE_NAME}_dedup.bam \
	-R $REFERENCE_GENOME \
	--known-sites $KNOWN_SITES_VCF \
	-O $ALIGN_DIR/${SAMPLE_NAME}_recal_data.table

log_message "Completed BQSR BaseRecalibrator for sample: $SAMPLE_NAME"

gatk ApplyBQSR \
	-R $REFERENCE_GENOME \
	-I $ALIGN_DIR/${SAMPLE_NAME}_dedup.bam \
	--bqsr-recal-file $ALIGN_DIR/${SAMPLE_NAME}_recal_data.table \
	-O $ALIGN_DIR/${SAMPLE_NAME}_recal.bam

log_message "Completed BQSR ApplyBQSR for sample: $SAMPLE_NAME"
