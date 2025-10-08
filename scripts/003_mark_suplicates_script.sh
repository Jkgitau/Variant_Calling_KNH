#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J mark-dup
#SBATCH -n 6
#SBATCH -o ../error_reports/mark_dup.out    # Standard output log
#SBATCH -e ../error_reports/mark_dups.err    # Standard error log

# Define Picard JAR path
PICARD_JAR="/export/apps/picard/2.8.2/picard.jar"

# Define directories
INPUT_DIR="/var/scratch/global/emurungi/variant_calling/entire_genome/hg19/aligned_bam_files"
OUTPUT_DIR="/var/scratch/global/emurungi/variant_calling/entire_genome/hg19/duplicates_marked"

# Create output directory if it does not exist
mkdir -p "$OUTPUT_DIR"


# Run MarkDuplicates on each BAM file
for BAM in "$INPUT_DIR"/*.bam; do
    BASENAME=$(basename "$BAM" .bam)
    java -jar "$PICARD_JAR" MarkDuplicates \
        INPUT="$BAM" \
        OUTPUT="$OUTPUT_DIR/${BASENAME}_dedup.bam" \
        METRICS_FILE="$OUTPUT_DIR/${BASENAME}_metrics.txt" \
        CREATE_INDEX=true \
        VALIDATION_STRINGENCY=SILENT
done
