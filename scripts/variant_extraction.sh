#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/reference.sh
source ../configs/filters.sh

echo "Running Variant Annotation..."
log_message "Running Variant Annotation..."

# Extract Fields from VCF to a Table
gatk VariantsToTable \
	-V $VCF_DIR/${SAMPLE_NAME}_annotated_snps.vcf \
	-F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O $VCF_DIR/${SAMPLE_NAME}_snps_table.txt
log_message "Extract Fields from VCF to a Table (SNP) for sample: $SAMPLE_NAME"

gatk VariantsToTable \
	-V $VCF_DIR/${SAMPLE_NAME}_annotated_indels.vcf \
	-F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O $VCF_DIR/${SAMPLE_NAME}_indels_table.txt
log_message "Extract Fields from VCF to a Table (INDELs) for sample: $SAMPLE_NAME"
