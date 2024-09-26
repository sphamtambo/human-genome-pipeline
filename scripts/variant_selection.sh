#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/reference.sh
source ../configs/filters.sh

echo "Running Variant Selection..."
log_message "Running Variant Selection..."

# Select Variants that Passed Filters (SNP)
gatk SelectVariants \
	-V $VCF_DIR/${SAMPLE_NAME}_filtered_snps.vcf \
	--exclude-filtered \
	-O $VCF_DIR/${SAMPLE_NAME}_pass_snps.vcf
log_message "Selected variant that passed filters (SNP) for sample: $SAMPLE_NAME"

# Select Variants that Passed Filters (INDELs)
gatk SelectVariants \
	-V $VCF_DIR/${SAMPLE_NAME}_filtered_indels.vcf \
	--exclude-filtered \
	-O $VCF_DIR/${SAMPLE_NAME}_pass_indels.vcf
log_message "Selected variant that passed filters (INDELS) for sample: $SAMPLE_NAME"

# Exclude Variants that Failed Genotype Filters (SNP)
gatk VariantFiltration \
	-R $REFERENCE_GENOME \
	-V $VCF_DIR/${SAMPLE_NAME}_pass_snps.vcf \
	--genotype-filter-expression "GQ < 20" \
	--genotype-filter-name "lowGQ" \
	-O $VCF_DIR/${SAMPLE_NAME}_genotype_filtered_snps.vcf
log_message "Excluded variants tha failed genotype filters (SNPs) for sample: $SAMPLE_NAME"

# Exclude Variants that Failed Genotype Filters (INDELS)
gatk VariantFiltration \
	-R $REFERENCE_GENOME \
	-V $VCF_DIR/${SAMPLE_NAME}_pass_indels.vcf \
	--genotype-filter-expression "GQ < 20" \
	--genotype-filter-name "lowGQ" \
	-O $VCF_DIR/${SAMPLE_NAME}_genotype_filtered_indels.vcf
log_message "Excluded variants tha failed genotype filters (INDELS) for sample: $SAMPLE_NAME"
