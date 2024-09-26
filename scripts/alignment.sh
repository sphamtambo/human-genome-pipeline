#!/usr/bin/env bash

# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/sample.sh
source ../configs/reference.sh

echo "Running Alignment..."
log_message "Running Alignment..."

# Index the reference genome with BWA
bwa index $REFERENCE_GENOME
log_message "Indexed reference genome with BWA: $REFERENCE_GENOME"

# Alignment with BWA
bwa mem -M -t $THREADS -R "$RG_STRING" $REFERENCE_GENOME $INPUT_DIR/$FORWARD_READS $INPUT_DIR/$REVERSE_READS >$ALIGN_DIR/${SAMPLE_NAME}.sam

log_message "Completed alignment of sample: $SAMPLE_NAME"
