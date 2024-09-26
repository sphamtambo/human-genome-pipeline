#!/usr/bin/env bash
# Load configurations
source ../configs/config.sh
source ../configs/paths.sh
source ../configs/sample.sh
source ../configs/reference.sh
source ../configs/download.sh

echo "Checking for downloaded files..."
log_message "Checking for downloaded files..."

# Check if forward reads are already downloaded
if [ -f "$INPUT_DIR/$FORWARD_READS" ]; then
	log_message "Forward reads already downloaded: $FORWARD_READS"
else
	echo "Downloading forward reads..."
	log_message "Downloading forward reads: $FORWARD_READS"
	wget -O $INPUT_DIR/$FORWARD_READS "$FORWARD_READ_URL"
	log_message "Downloaded forward reads: $FORWARD_READS"
fi

# Check if reverse reads are already downloaded
if [ -f "$INPUT_DIR/$REVERSE_READS" ]; then
	log_message "Reverse reads already downloaded: $REVERSE_READS"
else
	echo "Downloading reverse reads..."
	log_message "Downloading reverse reads: $REVERSE_READS"
	wget -O $INPUT_DIR/$REVERSE_READS "$REVERSE_READ_URL"
	log_message "Downloaded reverse reads: $REVERSE_READS"
fi

# Check if reference genome is already downloaded
# if [ -f "$REFERENCE_GENOME" ]; then
# 	log_message "Reference genome already downloaded: $REFERENCE_GENOME"
# else
# 	echo "Downloading reference genome..."
# 	log_message "Downloading reference genome: $REFERENCE_GENOME.gz"
# 	wget -O ${REFERENCE_GENOME}.gz "$REFERENCE_GENOME_URL"
# 	log_message "Downloaded reference genome: ${REFERENCE_GENOME}.gz"
# 	echo "Unzipping reference genome..."
# 	log_message "Unzipping reference genome: ${REFERENCE_GENOME}"
# 	gunzip ${REFERENCE_GENOME}.gz
# 	log_message "Unzipped reference genome: ${REFERENCE_GENOME}"
# fi
#
# # Check if reference genome index is already created
# if [ -f "$REFERENCE_INDEX" ]; then
# 	log_message "Reference genome index already created: $REFERENCE_INDEX"
# else
# 	echo "Creating reference genome index..."
# 	log_message "Creating reference genome index: $REFERENCE_INDEX"
# 	samtools faidx $REFERENCE_GENOME
# 	log_message "Created reference genome index: $REFERENCE_INDEX"
# fi
#
# Check if reference genome dictionary is already created
# if [ -f "$REFERENCE_DICT" ]; then
# 	log_message "Reference genome dictionary already created: $REFERENCE_DICT"
# else
# 	echo "Creating reference genome dictionary..."
# 	log_message "Creating reference genome dictionary: $REFERENCE_DICT"
# 	gatk CreateSequenceDictionary -R $REFERENCE_GENOME -O $REFERENCE_DICT
# 	log_message "Created reference genome dictionary: $REFERENCE_DICT"
# fi
#
# if [ -f "$KNOWN_SITES_VCF" ]; then
# 	log_message "Known sites VCF already downloaded: $KNOWN_SITES_VCF"
# else
# 	echo "Downloading known sites VCF..."
# 	log_message "Downloading known sites VCF: $KNOWN_SITES_VCF"
# 	wget -O $KNOWN_SITES_VCF "$DBSNP_VCF_URL"
# 	log_message "Downloaded known sites VCF: $KNOWN_SITES_VCF"
# fi
#
# if [ -f "$KNOWN_SITES_VCF_IDX" ]; then
# 	log_message "Known sites VCF index already downloaded: $KNOWN_SITES_VCF_IDX"
# else
# 	echo "Downloading known sites VCF index..."
# 	log_message "Downloading known sites VCF index: $KNOWN_SITES_VCF_IDX"
# 	wget -O $KNOWN_SITES_VCF_IDX "$DBSNP_VCF_IDX_URL"
# 	log_message "Downloaded known sites VCF index: $KNOWN_SITES_VCF_IDX"
# fi

# Check if Funcotator data sources are already downloaded
if [ -d "$FUNCOTATOR_DATA_SOURCES_DIR" ] && [ "$(ls -A $FUNCOTATOR_DATA_SOURCES_DIR)" ]; then
	log_message "Funcotator data sources already downloaded: $FUNCOTATOR_DATA_SOURCES_DIR"
else
	echo "Downloading funcotator..."
	log_message "Downloading Funcotator data sources for germline annotation..."
	gatk FuncotatorDataSourceDownloader \
		--germline \
		--validate-integrity \
		--extract-after-download \
		--hg38 \
		--output $FUNCOTATOR_DATA_SOURCES_DIR
	log_message "Funcotator germline data sources downloaded to $FUNCOTATOR_DATA_SOURCES_DIR"
fi

echo "All files ready to use."
log_message "All files ready to use."
