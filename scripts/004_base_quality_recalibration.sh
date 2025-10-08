#!/usr/bin/bash -l
#SBATCH -p batch                                   # Partition to submit to
#SBATCH -J variant_cal                             # Job name
#SBATCH -n 16                                      # Number of CPU cores requested
#SBATCH -o ../error_reports/base_cal_build_model.out   # Standard output log
#SBATCH -e ../error_reports/base_cal_build_model.err   # Standard error log

# Load required GATK module
module load gatk/4.4.0.0

# Define directories to Hg38 annotation files ## First analyses was based on these files
#INPUT_DIR="../../variant_calling/entire_genome/alignment"
#OUTPUT_DIR="../../variant_calling/entire_genome/base_quality_recalibration_with_dups"
#REFERENCE="../RB1_variant_calling/entire_genome/reference/hg38/Homo_sapiens_assembly38.fasta"
# KNOWN_SITES="../RB1_variant_calling/entire_genome/reference/Homo_sapiens_assembly38.dbsnp138.vcf"

# Define directories   ## Second batch was based on these after learning primers were designed based on hg19 reference genome
INPUT_DIR="../../variant_calling/entire_genome/hg19/aligned_bam_files"
OUTPUT_DIR="../../variant_calling/entire_genome/hg19/base_quality_recalibration"
REFERENCE="../../RB1_variant_calling/entire_genome/hg19_analysis/reference_genome/GCF_000001405.13_GRCh37_genomic.fna"
KNOWN_SITES="../../RB1_variant_calling/entire_genome/hg19_analysis/annotaion_files_hg19/Homo_sapiens_assembly19.dbsnp138.vcf"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over all BAM files in the input directory
for BAM in "$INPUT_DIR"/*.bam; do
    BASENAME=$(basename "$BAM" .bam)
    
    gatk BaseRecalibrator \
        -I "$BAM" \
        -R "$REFERENCE" \
        --known-sites "$KNOWN_SITES" \
        -O "$OUTPUT_DIR/${BASENAME}_recal_data.table"
done
