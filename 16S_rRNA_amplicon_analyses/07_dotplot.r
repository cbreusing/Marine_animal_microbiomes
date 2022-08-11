library(ggplot2)
library(gridExtra)
library(cowplot)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/microbiome")
data <- read.table("differential_taxa.txt", header=T)

col <- c("Migrator" = "paleturquoise4", "Non-Migrator" = "palevioletred4")

pdf("Differential_abundance_summary.pdf", width=8, height=10)
ggplot(subset(data, Count > 0), aes(x = Group, y = Family, size = Count)) + geom_point(aes(fill=Group), shape=21, color="gray32") + scale_size_continuous(range = c(2,10)) + theme_classic() + scale_fill_manual(values = col, name="Group") + theme(text = element_text(size = 15)) + labs(x = "", y = "Family") + scale_y_discrete(limits = rev) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
dev.off()

data2 <- read.table("differential_taxa_function.txt", header=T)

data2$Function <- factor(data2$Function, levels = c("Phototrophy", "Photoheterotrophy", "Chemoheterotrophy", "Aerobic Chemoheterotrophy", "Fermentation", "Chitinolysis", "Hydrocarbon Degradation", "Aliphatic Non-Methane Hydrocarbon Degradation", "Methylotrophy", "Methanol Oxidation", "Nitrogen Respiration", "Nitrate Respiration", "Nitrate Reduction", "Animal Parasites or Symbionts", "Intracellular Parasites", "Miscellaneous"))
col2 <- c("Phototrophy" = "midnightblue", "Photoheterotrophy" = "midnightblue", "Chemoheterotrophy" = "lightskyblue", "Aerobic Chemoheterotrophy" = "lightskyblue", "Fermentation" = "lightskyblue", "Chitinolysis" = "lightskyblue", "Hydrocarbon Degradation" = "lightskyblue", "Aliphatic Non-Methane Hydrocarbon Degradation" = "lightskyblue", "Methylotrophy" = "lightskyblue", "Methanol Oxidation" = "lightskyblue", "Nitrogen Respiration" = "darkseagreen", "Nitrate Respiration" = "darkseagreen", "Nitrate Reduction" = "darkseagreen", "Animal Parasites or Symbionts" = "darkorange", "Intracellular Parasites" = "darkorange", "Miscellaneous" = "papayawhip")

pdf("Differential_functions_summary.pdf", width=20, height=12)
plot1 <- ggplot(subset(data, Count > 0), aes(x = Group, y = Family, size = Count)) + geom_point(aes(fill=Group), shape=21, color="gray32") + scale_size_continuous(range = c(2,10)) + theme_classic() + scale_fill_manual(values = col, name="Group") + theme(text = element_text(size = 15)) + labs(x = "", y = "Family") + scale_y_discrete(limits = rev) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
plot2 <- ggplot(subset(data2, Count > 0), aes(x = Function, y = Family, size = Count)) + geom_point(aes(fill=Function), shape=21, color="gray32") + scale_size_continuous(range = c(2,10)) + theme_classic() + scale_fill_manual(values = col2, name="Function") + theme(text = element_text(size = 15)) + labs(x = "", y = "") + scale_y_discrete(limits = rev) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), axis.text.y = element_blank())
plot_grid(plot1, plot2, align = "h")
dev.off()