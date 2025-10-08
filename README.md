###  Variant Calling Analysis of 70 Kenyan Retinoblastoma Kids
A repository with Code, and results from the Kenya RB1 Study


### Note:
* To re-use these Analysis scripts, you have to update your file paths accordingly.
* ake a copy of raw Fastq files before starting manipulating the files (Back-up).


------------------------------------------------------------------------------------------------
### These are important links used in this Analysis

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

