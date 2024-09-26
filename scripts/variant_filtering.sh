#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/reference.sh
source ../configs/filters.sh

echo "Running Variant Filtering..."
log_message "Running Variant Filtering..."

# Filter SNPS
gatk VariantFiltration \
	-R $REFERENCE_GENOME \
	-V $VCF_DIR/${SAMPLE_NAME}_raw_snps.vcf \
	-O $VCF_DIR/${SAMPLE_NAME}_filtered_snps.vcf \
	-filter-name "QD_filter" --filter-expression "${SNP_QD_FILTER_EXPR}" \
	-filter-name "FS_filter" --filter-expression "${SNP_FS_FILTER_EXPR}" \
	-filter-name "MQ_filter" --filter-expression "${SNP_MQ_FILTER_EXPR}" \
	-filter-name "SOR_filter" --filter-expression "${SNP_SOR_FILTER_EXPR}" \
	-filter-name "MQRankSum_filter" --filter-expression "${SNP_MQRankSum_FILTER_EXPR}" \
	-filter-name "ReadPosRankSum_filter" --filter-expression "${SNP_ReadPosRankSum_FILTER_EXPR}" \
	-genotype-filter-expression "${DP_GENOTYPE_FILTER_EXPR}" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "${GQ_GENOTYPE_FILTER_EXPR}" \
	-genotype-filter-name "GQ_filter"
log_message "Applied SNP filter for sample: $SAMPLE_NAME"

# Filter INDELs
gatk VariantFiltration \
	-R $REFERENCE_GENOME \
	-V $VCF_DIR/${SAMPLE_NAME}_raw_indels.vcf \
	-O $VCF_DIR/${SAMPLE_NAME}_filtered_indels.vcf \
	-filter-name "QD_filter" --filter-expression "${INDEL_QD_FILTER_EXPR}" \
	-filter-name "FS_filter" --filter-expression "${INDEL_FS_FILTER_EXPR}" \
	-filter-name "SOR_filter" --filter-expression "${INDEL_SOR_FILTER_EXPR}" \
	-genotype-filter-expression "${DP_GENOTYPE_FILTER_EXPR}" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "${GQ_GENOTYPE_FILTER_EXPR}" \
	-genotype-filter-name "GQ_filter"
log_message "Applied INDEL filter for sample: $SAMPLE_NAME"
