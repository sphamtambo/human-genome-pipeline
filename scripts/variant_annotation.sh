#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/reference.sh
source ../configs/filters.sh

echo "Running Variant Annotation..."
log_message "Running Variant Annotation..."

# Step 12: Annotate Variants using Funcotator
gatk Funcotator \
	--variant $VCF_DIR/${SAMPLE_NAME}_genotype_filtered_snps.vcf \
	--reference $REFERENCE_GENOME \
	--ref-version $FUNCOTATOR_REF_VERSION \
	--data-sources-path $FUNCOTATOR_DATA_SOURCES \
	--output $VCF_DIR/${SAMPLE_NAME}_annotated_snps.vcf \
	--output-file-format $FUNCOTATOR_OUTPUT_FORMAT
log_message "Annotated varaint (SNP) for sample: $SAMPLE_NAME"

gatk Funcotator \
	--variant $VCF_DIR/${SAMPLE_NAME}_genotype_filtered_indels.vcf \
	--reference $REFERENCE_GENOME \
	--ref-version $FUNCOTATOR_REF_VERSION \
	--data-sources-path $FUNCOTATOR_DATA_SOURCES \
	--output $VCF_DIR/${SAMPLE_NAME}_annotated_indels.vcf \
	--output-file-format $FUNCOTATOR_OUTPUT_FORMAT
log_message "Annotated varaint (INDELS) for sample: $SAMPLE_NAME"
