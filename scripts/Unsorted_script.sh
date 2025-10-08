
# Path to fastq raw Fastq files
/home/emurungi/gitau/Doctorat/RB1_variant_calling/RB_data/raw_fastq_data

Path to BWA indexed reference genome
/home/emurungi/gitau/Doctorat/RB1_variant_calling/entire_genome/hg19_analysis/reference_genome
pattern: GCF_000001405.13_GRCh37_genomic.fna.*

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



---------------------------------------------------------------------------------------------------------------------------------
# Command to calculate sequencing coverage per sample using bedtools - to get the mean reads per exon 
bedtools coverage -a RB1_exons_coordinates_hg19.bed -b ../aligned_bam_files/C001_S1.sorted.bam -mean > coverage_C001.mean.txt

# Once we have the coverage, we use  to plot and save the results  ## This can be automated using loops to process one file at a time
library(ggplot2)
library(viridis)

# read mean coverage
mean_cov <- read.table("coverage_C001.mean.txt", header=FALSE, fill=TRUE)
colnames(mean_cov) <- c("Chromosome","Start coordinate","End coordinate","Exon number","Mean depth")

# order exons by genomic coordinate
mean_cov$Exon_Number <- factor(mean_cov$Exon_Number, levels=unique(mean_cov$Exon_Number[order(mean_cov$Start_coordinate)]))

# Line Plot
ggplot(mean_cov, aes(x=Exon_Number, y=Mean_Depth, group=1)) +
  geom_line(color="black", size=1) +
  geom_point(color="darkred", size=2) +
  theme_minimal() +
  labs(title="Mean sequencing depth across RB1 exons in sample C001",
       x="RBI Exons", y="Mean Sequencing Coverage (Sample C001)") +
  theme(axis.text.x = element_text(angle=45, hjust=1))

# BarPlots

ggplot(mean_cov, aes(x=factor(Exon_Number), y=Mean_Depth)) +
  geom_col(fill=viridis(1, option="D"), color="red") +  # color-blind friendly green-blue
  labs(title="Mean sequencing depth across RB1 exons in sample C001",
       x="RB1 Exons", y="Mean Sequencing Coverage (Sample C001)") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=14),
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    plot.title = element_text(hjust = 0.5, size=16, face="bold")
  )


# Save the coverage dataframe for C001 sample
write_tsv(mean_cov, file = 'C001_Mean_coverage.tsv')

