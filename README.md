###  Variant Calling Analysis of 70 Kenyan Retinoblastoma Kids
A repository with Code, and results from the Kenya RB1 Study


### Note:
* To re-use these Analysis scripts, you have to update your file paths accordingly.
* ake a copy of raw Fastq files before starting manipulating the files (Back-up).


------------------------------------------------------------------------------------------------
### Some useful Links 

High confidence SNPS and annotation files: https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false

Reference Analysis steps (Melbrone Univesity): https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false


------------------------------------------------------------------------------------------------

#### Step 1: Quality Assessment of Fastq Files
Read More here: https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/QC_raw_data.html

We used FastQC to assess sequencing quality, including potential errors, duplicate reads, GC content, and overrepresented sequences. Reads with a Phred score below 20 were removed before further processing. After this cleaning step, the data were re-evaluated with FastQC to confirm that the retained reads were of high quality.

```
#!/usr/bin/bash -l                                # Shebang
#SBATCH -p batch                                  # Partition to submit to
#SBATCH -J fastqc                                 # Job name
#SBATCH -n 8                                      # Number of CPU cores requested

# Load modules
fastqc/0.11.9

# Define file directories
INPUT_DIR="../../../RB_data/trimmed_fastq_data/copy_of_trimmed_reads/"
OUTPUT_DIR="../../variant_calling/entire_genome/hg19/fastq_files"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Loop through FASTQ files
for file in "$INPUT_DIR"/*.fastq.gz
do
    fastqc "$file" -o "$OUTPUT_DIR"
done
```

#### Step 2: Aligning reads to the reference genome
Clean reads were aligned to the hg19 reference genome to determine their genomic origin. Although this was a targeted study, reads were first mapped to the whole genome, and variants were later restricted to chromosome 13 exons using a BED file.

Importance: ensures accurate alignment and variant calling in the regions of interest.

```
#!/usr/bin/bash -l                 #Shebang
#SBATCH -p batch                   # Partition to submit to
#SBATCH -J align_and_sort          # Job name
#SBATCH -n 12                      # Number of CPU cores requested
#SBATCH -o align_and_sort.out      # Standard output log
#SBATCH -e align_and_sort.err      # Standard error log

# Load required modules
#module load bwa/0.7.18
module load bwa/0.7.19
module load samtools/1.17

# Define paths to hg38 reference genome
#REF="../..reference/RB1_variant_calling/entire_genome/reference/hg38/Homo_sapiens_assembly38.fasta"
#OUTPUT_DIR="..results/variant_calling/entire_genome/alignment"
#INPUT_DIR="..input/RB1_variant_calling/RB_data/trimmed_fastq_data"

REF="../..reference/reference/RB1_variant_calling/entire_genome/hg19_analysis/reference_genome/GCF_000001405.13_GRCh37_genomic.fna"
OUTPUT_DIR="..results/variant_calling/entire_genome/hg19/aligned_bam_files"
INPUT_DIR="../..input/RB1_variant_calling/RB_data/trimmed_fastq_data"

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
```

#### Step 3: Marking duplicates
Duplicate reads, often introduced during PCR amplification in targeted sequencing, can artificially inflate coverage and bias variant calling. Marking duplicates ensures that only unique molecules contribute to variant detection, improving the accuracy and reliability of downstream analyses.

Get GATK standard high confidence SNPS and annotation files (here)[ https://console.cloud.google.com/storage/browser/gcp-public-data--broad-references/hg19/v0;tab=objects?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&prefix=&forceOnObjectsSortingFiltering=false]

```
#!/usr/bin/bash -l
#SBATCH -p batch
#SBATCH -J mark-dup
#SBATCH -n 6
#SBATCH -o ../error_reports/mark_dup.out    # Standard output log
#SBATCH -e ../error_reports/mark_dups.err    # Standard error log

# Define Picard JAR path
PICARD_JAR="/export/apps/picard/2.8.2/picard.jar"

# Define directories
INPUT_DIR="/input_files/variant_calling/entire_genome/hg19/aligned_bam_files"
OUTPUT_DIR="/results/variant_calling/entire_genome/hg19/duplicates_marked"

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
```

#### Step 4a: Base Quality Recalibration
Base Quality Score Recalibration (BQSR) corrects systematic errors made by the sequencer when estimating base call accuracy. In targeted variant calling, this step is important because even small biases in quality scores can lead to false positives or missed variants in the focused regions of interest, thereby improving the precision of variant detection.

```
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
```

