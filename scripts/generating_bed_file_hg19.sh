#!/bin/bash

# Where did I get hg19 bed file?
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.13_GRCh37/

# from the above link download a gtf file (GCF_000001405.13_GRCh37_genomic.gtf)

# Extract all entries with 'exon'
grep -w "RB1" GCF_000001405.13_GRCh37_genomic.gtf | grep "exon" > RB1_exons.gtf

# Use awk to grab the 27 exons that code for mRNA
awk -F'\t' '$3=="exon" && $9 ~ /gene_id "RB1"/ {
    match($9, /exon_number "([0-9]+)"/, a);
    print $1"\t"$4"\t"$5"\t"a[1]
}' RB1_exons.gtf > RB1_exons_coordinates.bed  ## Confirm these coordinates by oppening the gft file and checking the coordinates of the exons
                                            ## This file will be passed to sarek to restrict variant calling to these regions


