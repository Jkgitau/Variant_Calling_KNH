#!/usr/bin/bash -l
#SBATCH -p batch                                   # Partition to submit to
#SBATCH -J variant_cal                             # Job name
#SBATCH -n 16                                      # Number of CPU cores requested
#SBATCH -o ../error_reports/Apply_bqsr.out         # Standard output log
#SBATCH -e ../error_reports/Apply_bqsr.err         # Standard error log


# Load required GATK module
module load gatk/4.4.0.0


# Define directories
INPUT_DIR="../../variant_calling/entire_genome/marked_duplicates"
OUTPUT_DIR="../../variant_calling/entire_genome/base_quality_recalibration_step_2"
REFERENCE="../RB1_variant_calling/entire_genome/reference/hg38/Homo_sapiens_assembly38.fasta"
BQSR_DIR="../variant_calling/entire_genome/base_quality_recalibration"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through all BAM files
for bam_file in "$INPUT_DIR"/*.sorted_dedup.bam; do
    # Extract the sample base (remove the .sorted_dedup.bam part)
    sample_base=$(basename "$bam_file" .sorted_dedup.bam)

    echo "Processing sample: $sample_base"

    gatk ApplyBQSR \
        -I "$bam_file" \
        -R "$REFERENCE" \
        --bqsr-recal-file "$BQSR_DIR/${sample_base}.sorted_recal_data.table" \
        -O "$OUTPUT_DIR/${sample_base}.sorted_dedup.bqsr.bam"
done

echo "All samples processed!"
