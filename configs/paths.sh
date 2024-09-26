#!/usr/bin/env bash

# Directories
INPUT_DIR="../data"
QC_DIR="../qc"
ALIGN_DIR="../alignment"
METRICS_DIR="../metrics"
VCF_DIR="../vcf"
LOG_DIR="../logs"
REFERENCE_DIR="../reference"
FUNCOTATOR_DATA_SOURCES_DIR="../funcotator"

# Ensure directories exist
mkdir -p $INPUT_DIR $QC_DIR $REFERENCE_DIR $ALIGN_DIR $METRICS_DIR $VCF_DIR $LOG_DIR $FUNCOTATOR_DATA_SOURCES_DIR

# Logging settings
LOG_LEVEL="INFO"
LOG_FILE="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"

# Logging function
log_message() {
	local message="$1"
	echo "[$(date +'%Y-%m-%d %H:%M:%S')] [$LOG_LEVEL] $message" | tee -a $LOG_FILE
}
