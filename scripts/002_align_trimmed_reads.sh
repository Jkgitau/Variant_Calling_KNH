#!/usr/bin/bash -l
#SBATCH -p batch                                  # Partition to submit to
#SBATCH -J align_and_sort                         # Job name
#SBATCH -n 12                                     # Number of CPU cores requested
#SBATCH -o align_and_sort.out      # Standard output log
#SBATCH -e align_and_sort.err      # Standard error log

# Load modules
#module load bwa/0.7.18
module load bwa/0.7.19
module load samtools/1.17

# Define paths
#REF="../../RB1_variant_calling/entire_genome/reference/hg38/Homo_sapiens_assembly38.fasta"
#OUTPUT_DIR="../variant_calling/entire_genome/alignment"
#INPUT_DIR="../RB1_variant_calling/RB_data/trimmed_fastq_data"

REF="/home/emurungi/gitau/Doctorat/RB1_variant_calling/entire_genome/hg19_analysis/reference_genome/GCF_000001405.13_GRCh37_genomic.fna"
OUTPUT_DIR="/var/scratch/global/emurungi/variant_calling/entire_genome/hg19/aligned_bam_files"
INPUT_DIR="/home/emurungi/gitau/Doctorat/RB1_variant_calling/RB_data/trimmed_fastq_data"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop over all R1 FASTQ files
for R1 in "$INPUT_DIR"/*_trimmed_R1.fastq.gz; do
    # Get base sample name (e.g., C001_S1)
    SAMPLE=$(basename "$R1" _L001_trimmed_R1.fastq.gz)

    # Construct R2 file path
    R2="$INPUT_DIR/${SAMPLE}_L001_trimmed_R2.fastq.gz"

    # Check if R2 file exists
    if [[ ! -f "$R2" ]]; then
        echo "Missing R2 for $SAMPLE, skipping..."
        continue
    fi

    # Output sorted BAM file
    OUT_BAM="$OUTPUT_DIR/${SAMPLE}.sorted.bam"

    # Define read group
    RG="@RG\tID:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA\tPM:MiSeq\tSM:${SAMPLE}"

    echo "Aligning and sorting $SAMPLE..."

    # Align, convert to BAM, and sort
    bwa mem -t 12 -M -R "$RG" "$REF" "$R1" "$R2" | \
    samtools view -b -h | \
    samtools sort -@ 12 -o "$OUT_BAM"

    # Optionally index the BAM file
    samtools index "$OUT_BAM"

    echo "Finished processing $SAMPLE."
done



