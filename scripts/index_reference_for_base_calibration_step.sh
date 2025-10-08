#!/bin/bash


# Load required modules
module load samtools/1.9
module load gatk/4.4.0.0

# Use samtools to index the reference genome, generating neeed .fai file
samtools faidx samtools faidx Homo_sapiens_assembly19.fasta

# Use GATK to create sequence dictionary
gatk CreateSequenceDictionary \
     -R Homo_sapiens_assembly19.fasta \
     -O Homo_sapiens_assembly19.dict



