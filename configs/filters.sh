#!/usr/bin/env bash

# SNP Filter Expressions
SNP_QD_FILTER_EXPR="QD < 2.0"
SNP_FS_FILTER_EXPR="FS > 60.0"
SNP_MQ_FILTER_EXPR="MQ < 40.0"
SNP_SOR_FILTER_EXPR="SOR > 4.0"
SNP_MQRankSum_FILTER_EXPR="MQRankSum < -12.5"
SNP_ReadPosRankSum_FILTER_EXPR="ReadPosRankSum < -8.0"

# INDEL Filter Expressions
INDEL_QD_FILTER_EXPR="QD < 2.0"
INDEL_FS_FILTER_EXPR="FS > 200.0"
INDEL_SOR_FILTER_EXPR="SOR > 10.0"

# Genotype Filter Expressions (common for both SNPs and INDELs)
DP_GENOTYPE_FILTER_EXPR="DP < 10"
GQ_GENOTYPE_FILTER_EXPR="GQ < 10"
