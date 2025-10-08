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
