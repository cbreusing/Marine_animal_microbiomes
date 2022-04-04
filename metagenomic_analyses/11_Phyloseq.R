library(ggplot2)
library(plyr)
library(dplyr)
library(devtools)
library(phyloseq)
library(ape)
library(vegan)
library(RColorBrewer)
library(circlize)
library(rbiom)
library(viridis)
library(psych)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Marine_animal_microbiomes_UCSD/metagenomics/SHOGUN/utree")

# Import final ASV biom and mapping files
biomFile <- import_biom("taxatable.utree.species.biom", parseFunction = parse_taxonomy_default)
mapFile <- import_qiime_sample_data("Smithsonian_metadata.txt")
marmic <- merge_phyloseq(biomFile, mapFile)

colnames(tax_table(marmic)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove low abundance ASVs and all samples with less than 1000 reads
marmic = prune_samples(sample_sums(marmic) > 1000, marmic)
marmic = filter_taxa(marmic, function(x) sum(x) > 10, TRUE)
marmic = filter_taxa(marmic, function (x) {sum(x > 0) > 1}, TRUE)

sample_data(marmic)$Depth <- factor(sample_data(marmic)$Depth, levels = c("242", "308", "817-862", "1120-1150"))
sample_data(marmic)$Migration <- factor(sample_data(marmic)$Migration, levels = c("Yes", "No", "Unknown"))

# Transform counts to proportional abundances, this has recently been shown to be the most appropriate normalization technique for ecological microbiome analyses such as alpha and beta diversity (when using abundance based distance methods such as Bray-Curtis or weighted UniFrac)
marmictrans <- transform_sample_counts(marmic, function(x) x/sum(x))

otu <- t(otu_table(marmic))
spec <- specaccum(otu)

pdf("Specaccum.filtered.pdf")
plot(spec, ci.type="poly", col = "black", lwd=2, ci.lty=0, ci.col="grey", ci = 1.96, xlab = "# individuals", ylab = "# OTUs")
dev.off()

pdf("Fractional_abundance_a.filtered.pdf", width=24, height=14)
plot_bar(marmictrans, x= "Sample", fill = "Phylum") + facet_grid(. ~ Group, scales = "free", space = "free") + geom_bar(aes(color=Phylum, fill=Phylum), stat='identity', position='stack') + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_text(size = 15, face="bold")) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom") + scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"))
dev.off()

pdf("Fractional_abundance_b.filtered.pdf", width=24, height=14)
plot_bar(marmictrans, x= "Sample", fill = "Class") + facet_grid(. ~ Group, scales = "free", space = "free") + geom_bar(aes(color=Class, fill=Class), stat='identity', position='stack') + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_text(size = 15, face="bold")) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom") + scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"))
dev.off()

# Group
colors <- c("Fish" = "cornflowerblue", "Isopod" = "pink1", "Krill" = "red1", "Polychaete" = "goldenrod1", "Shrimp" = "maroon4", "Cephalopod" = "midnightblue")
# Species
colors2 <- c("Acanthamunnopsis milleri" = "mistyrose", "Cyclothone atraria" = "lightblue1", "Cyclothone signata" = "lightsteelblue2", "Munneurycope murrayi" = "pink2", "Poeobius meseres" = "lightgoldenrod1", "Eusergestes similis" = "maroon4", "Stenobrachius leucopsarus" = "royalblue1", "Euphausia pacifica" = "red1", "Tomopteris sp." = "darkgoldenrod1", "Vampyroteuthis infernalis" = "midnightblue")
# Migration
colors3 <- c("Yes" = "paleturquoise4", "No" = "palevioletred4", "Unknown" = "darkgrey")
# Diet
colors4 <- c("Detritus" = "sandybrown", "Phytoplankton" = "darkgreen", "Zooplankton" = "saddlebrown")
# Depth
colors5 <- c("242" = "#eff3ff", "308" = "#bdd7e7", "452" = "#6baed6", "817-862" = "#2171b5", "1120-1150" = "#013270")
# Reads
colors6 <- c("<1000" = "#eff3ff", ">1000" = "#013270")


# PCoA ordination plots
PCoA2 <- ordinate(marmictrans, "PCoA", "bray")
p2 <- plot_ordination(marmictrans, PCoA2)

pdf("PCoA_Bray_a.filtered.pdf")
p2 + scale_fill_manual(values = colors, name="Group") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Group)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_b.filtered.pdf")
p2 + scale_fill_manual(values = colors2, name="Species") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + theme(legend.text = element_text(face = "italic"))
dev.off()

pdf("PCoA_Bray_c.filtered.pdf")
p2 + scale_fill_manual(values = colors3, name="Diel Vertical Migration") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Migration)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_d.filtered.pdf")
p2 + scale_fill_manual(values = colors4, name="Diet") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Diet)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_e.filtered.pdf")
p2 + scale_fill_manual(values = colors5, name="Median Depth (m)") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_f.filtered.pdf")
p2 + scale_fill_manual(values = colors6, name="Reads") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Reads)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_g.filtered.pdf")
p2 + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species), shape = Migration, color = factor(Species)), size=3) + scale_shape_manual(values = c(19, 17, 15), name="Diel Vertical Migration") + scale_fill_manual(values = colors2, name="Species", aesthetics = c("colour", "fill")) + theme(text = element_text(size = 15)) + theme(legend.text = element_text(face = "italic"))
dev.off()

pdf("PCoA_Bray_h.filtered.pdf")
p2 + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth), shape = Diet, color = factor(Depth)), size=3) + scale_shape_manual(values = c(19, 17, 15), name="Diet") + scale_fill_manual(values = colors5, name="Median Depth (m)", aesthetics = c("colour", "fill")) + theme(text = element_text(size = 15))
dev.off()


# Traditional PERMANOVA
marmicbray <- distance(marmictrans, method = "bray")
sampledf <- data.frame(sample_data(marmic))

sink("PERMANOVA.txt")
# Adonis test
adonis(marmicbray ~ Species, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Diet, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Migration, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Depth, data = sampledf, add = "cailliez", permutations = 999)

beta1 <- betadisper(marmicbray, sampledf$Species)
permutest(beta1)

beta2 <- betadisper(marmicbray, sampledf$Diet)
permutest(beta2)

beta3 <- betadisper(marmicbray, sampledf$Migration)
permutest(beta3)

beta4 <- betadisper(marmicbray, sampledf$Depth)
permutest(beta4)
sink()
