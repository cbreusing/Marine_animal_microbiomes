library(ggplot2)

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/microbiome")
data <- read.table("differential_taxa.txt", header=T)

col <- c("Migrator" = "paleturquoise4", "Non-Migrator" = "palevioletred4")

pdf("Differential_abundance_summary.pdf", width=8, height=10)
ggplot(subset(data, Count > 0), aes(x = Group, y = Family, size = Count)) + geom_point(aes(fill=Group), shape=21, color="gray32") + scale_size_continuous(range = c(2,10)) + theme_classic() + scale_fill_manual(values = col, name="Group") + theme(text = element_text(size = 15)) + labs(x = "", y = "Family") + scale_y_discrete(limits = rev) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()