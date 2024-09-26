#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/reference.sh

echo "Running Variant Calling..."
log_message "Running Variant Calling..."

# Variant Calling with HaplotypeCaller
gatk HaplotypeCaller \
	-R $REFERENCE_GENOME \
	-I $ALIGN_DIR/${SAMPLE_NAME}_sorted_dedup_reads.bam \
	-O $VCF_DIR/${SAMPLE_NAME}_raw_variants.vcf
log_message "Completed Variant Calling for sample: $SAMPLE_NAME"

# Select Variants (SNP)
gatk SelectVariants \
	-R $REFERENCE_GENOME \
	-V $VCF_DIR/${SAMPLE_NAME}_raw_variants.vcf \
	--select-type SNP \
	-O $VCF_DIR/${SAMPLE_NAME}_raw_snps.vcf
log_message "Selected variant (SNP) for sample: $SAMPLE_NAME"

# Select Variants (INDELS)
gatk SelectVariants \
	-R $REFERENCE_GENOME \
	-V $VCF_DIR/${SAMPLE_NAME}_raw_variants.vcf \
	--select-type INDEL \
	-O $VCF_DIR/${SAMPLE_NAME}_raw_indels.vcf
log_message "Selected variant (SNP) for sample: $SAMPLE_NAME"
