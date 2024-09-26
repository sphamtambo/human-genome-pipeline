# GATK Pipeline for Whole Genome Sequencing

This pipeline implements the GATK (Genome Analysis Toolkit) Best Practices workflow for germline SNP and Indel discovery in whole genome sequencing (WGS) data. It provides a comprehensive solution for variant calling and filtration, designed to ensure high-quality results from raw sequencing data to filtered variant calls.

The pipeline follows these key steps, adhering to GATK best practices:

1. Quality Control of Raw Sequencing Data (FastQC)
2. Mapping to Reference Genome (BWA-MEM)
3. Marking Duplicates (MarkDuplicates)
4. Base Quality Score Recalibration (BQSR)
5. Calling Variants (HaplotypeCaller)
6. Variant Filtration (Hard Filters)
7. Variant Annotation (Funcotator)

This workflow is optimized for high-throughput whole genome sequencing data and incorporates the latest recommendations from the GATK team to maximize the accuracy and reliability of variant calls.

## Prerequisites

- Bash shell
- Required bioinformatics tools:
  - FastQC
  - BWA
  - samtools
  - GATK

## Usage

Follow these steps to set up and run the GATK pipeline:

1. Clone the repository:

   ```bash
   git clone https://github.com/sphamtambo/human-genome-pipeline.git
   cd human-genome-pipeline
   ```

2. Set up the configuration files:

   - Navigate to the `config` directory:
     ```bash
     cd config
     ```
   - Edit the following configuration files according to your needs:
     - `config.sh`: Set general pipeline parameters (e.g., number of threads)
     - `paths.sh`: Define paths for input data, output directories, and reference files
     - `sample.sh`: Specify sample-specific details (e.g., sample name, FASTQ file names)
     - `reference.sh`: Provide information about the reference genome
     - `filters.sh`: Define SNP and INDEL filter expressions
     - `download.sh`: Specify URLs or paths for downloading necessary data (if not providing your own)
     - `funcotator.sh`: Configure settings for functional annotation, including the data source

3. Prepare your input data:
   You have two options for providing input data:

   a. Use your own data:

   - Place your FASTQ files in the input directory specified in `paths.sh`
   - Ensure your reference genome and other required files are in the locations specified in `reference.sh`
   - If using Funcotator, place your Funcotator data source in the directory specified in `funcotator.sh`

   b. Use downloaded data:

   - In `download.sh`, provide URLs for the input FASTQ files, reference genome, and any other required files
   - If using Funcotator, provide the URL for the Funcotator data source in `funcotator.sh`
   - Run the download script to fetch the necessary data:
     ```bash
     bash scripts/download_data.sh
     ```

4. Run the pipeline steps:

   - Navigate to the `scripts` directory:
     ```bash
     cd ../scripts
     ```
   - Execute each script in the following order, adjusting as necessary for your specific analysis needs:

     a. Quality Control:

     ```bash
     ./quality_control.sh
     ```

     b. Alignment:

     ```bash
     ./alignment.sh
     ```

     c. Deduplication:

     ```bash
     ./deduplication.sh
     ```

     d. Base Quality Score Recalibration:

     ```bash
     ./base_quality_score_recalibration.sh
     ```

     e. Alignment Metrics:

     ```bash
     ./alignment_metrics.sh
     ```

     f. Variant Calling:

     ```bash
     ./variant_calling.sh
     ```

     g. Variant Filtering:

     ```bash
     ./variant_filtering.sh
     ```

     h. Variant Selection:

     ```bash
     ./variant_selection.sh
     ```

     i. Variant Annotation:

     ```bash
     ./variant_annotation.sh
     ```

     j. Variant Extraction:

     ```bash
     ./variant_extraction.sh
     ```

   - Monitor the output of each script for progress updates and any error messages.
   - Ensure each step completes successfully before moving to the next one.

5. Retrieve results:
   - Once all steps are complete, find your results in the output directory specified in `paths.sh`
   - Key output files typically include:
     - Aligned BAM files
     - VCF files containing called variants
     - Various quality control and metrics reports
     - Annotated variant files

## Logging

Each script logs its activities to a separate log file in the directory specified in `paths.sh`. Check these logs for detailed information about each step's execution.
