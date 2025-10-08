#!/usr/bin/bash -l
#SBATCH -p batch                                  # Partition to submit to
#SBATCH -J fastqc                                 # Job name
#SBATCH -n 8                                      # Number of CPU cores requested
#SBATCH -o ../error_reports/fastqc.out   # standard output log
#SBATCH -e ../error_reports/fastqc.err   # standard error log

# Load modules
fastqc/0.11.9

# Define file directories
INPUT_DIR="../../../RB_data/trimmed_fastq_data/copy_of_trimmed_reads/"
OUTPUT_DIR="/var/scratch/global/emurungi/variant_calling/entire_genome/hg19/fastq_files"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through FASTQ files
for file in "$INPUT_DIR"/*.fastq.gz
do
    fastqc "$file" -o "$OUTPUT_DIR"
done

