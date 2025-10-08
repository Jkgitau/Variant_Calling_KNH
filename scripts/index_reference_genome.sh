#!/usr/bin/bash -l
#SBATCH -p batch                                  # Partition to submit to
#SBATCH -J indexing                         # Job name
#SBATCH -n 8                                      # Number of CPU cores requested


# Load modules
module load bwa/0.7.19

# Specify the indexing command
bwa index ../reference_genome/GCF_000001405.13_GRCh37_genomic.fna
