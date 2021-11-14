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

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/diet")

# Import final ASV biom and mapping files
biomFile <- import_biom("zotu-table-diet_tax.biom", parseFunction = parse_taxonomy_default)
mapFile <- import_qiime_sample_data("Smithsonian_metadata.txt")
biomMapFile <- merge_phyloseq(biomFile, mapFile)

# Import representative sequences and remove non-ASV information
repsetFile <- Biostrings::readDNAStringSet("dna-sequences.fasta")
names(repsetFile) <- gsub("\\s.+$", "", names(repsetFile))

# Import tree file
treefile <- read_tree("tree.nwk")
new_tre <- ape::multi2di(treefile)

biomMapTree <- merge_phyloseq(biomMapFile, new_tre)

# Create full phyloseq object
phyloseq <- merge_phyloseq(biomMapTree, repsetFile)
colnames(tax_table(phyloseq)) = c("Domain", "Category", "Phylum", "Class", "Order", "Family", "Genus", "Species")
#colnames(tax_table(phyloseq)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove low abundance ASVs and all samples with less than 1000 reads (recommended); for diet I included all samples with >500 reads
marmic = prune_samples(sample_sums(phyloseq) > 500, phyloseq)
marmic = filter_taxa(marmic, function(x) sum(x) > 10, TRUE)

#sample_data(marmic)$Depth <- factor(sample_data(marmic)$Depth, levels = c("239-244", "278-337", "448-455", "500-900", "814-909", "540-1700", "500-1800"))
sample_data(marmic)$Depth <- factor(sample_data(marmic)$Depth, levels = c("239-244", "278-337", "448-455", "500-900", "814-909", "500-1800"))
sample_data(marmic)$Migration <- factor(sample_data(marmic)$Migration, levels = c("Yes", "No", "Unknown"))

# Transform counts to proportional abundances, this has recently been shown to be the most appropriate normalization technique for ecological microbiome analyses such as alpha and beta diversity (when using abundance based distance methods such as Bray-Curtis or weighted UniFrac)
marmictrans <- transform_sample_counts(marmic, function(x) x/sum(x))

otu <- t(otu_table(marmic))
spec <- specaccum(otu)

pdf("Specaccum.pdf")
plot(spec, ci.type="poly", col = "black", lwd=2, ci.lty=0, ci.col="grey", ci = 1.96, xlab = "# individuals", ylab = "# ASVs")
dev.off()

pdf("Fractional_abundance_a.pdf", width=20, height=14)
plot_bar(marmictrans, x= "Sample", fill = "Phylum") + facet_grid(. ~ Group, scales = "free", space = "free") + geom_bar(aes(color=Phylum, fill=Phylum), stat='identity', position='stack') + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_text(size = 15, face="bold")) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom") + scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("Fractional_abundance_b.pdf", width=20, height=14)
plot_bar(marmictrans, x= "Sample", fill = "Class") + facet_grid(. ~ Group, scales = "free", space = "free") + geom_bar(aes(color=Class, fill=Class), stat='identity', position='stack') + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_text(size = 15, face="bold")) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom") + scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

pdf("Fractional_abundance_c.pdf", width=20, height=14)
plot_bar(marmictrans, x= "Sample", fill = "Category") + facet_grid(. ~ Group, scales = "free", space = "free") + geom_bar(aes(color=Category, fill=Category), stat='identity', position='stack') + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_text(size = 15, face="bold")) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom") + scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()


# Group
#colors <- c("Fish" = "cornflowerblue", "Isopod" = "pink1", "Krill" = "red1", "Polychaete" = "goldenrod1", "Shrimp" = "maroon4", "Squid" = "midnightblue", "Tunicate" = "cyan4")
colors <- c("Fish" = "cornflowerblue", "Isopod" = "pink1", "Krill" = "red1", "Polychaete" = "goldenrod1", "Shrimp" = "maroon4", "Tunicate" = "cyan4")
# Species
#colors2 <- c("Acanthamunnopsis" = "mistyrose", "Cyclothone" = "lightblue1", "Krill" = "red1", "Munneurycope" = "pink2", "Myctophid" = "royalblue1", "Poeobius" = "lightgoldenrod1", "Sergestes" = "maroon4", "Tomopteris" = "darkgoldenrod1", "Vampyroteuthis" = "midnightblue", "Vitreosalpa" = "cyan4")
colors2 <- c("Acanthamunnopsis" = "mistyrose", "Cyclothone" = "lightblue1", "Krill" = "red1", "Munneurycope" = "pink2", "Myctophid" = "royalblue1", "Poeobius" = "lightgoldenrod1", "Sergestes" = "maroon4", "Vitreosalpa" = "cyan4")
# Migration
colors3 <- c("Yes" = "paleturquoise4", "No" = "palevioletred4", "Unknown" = "darkgrey")
# Diet
colors4 <- c("Detritus" = "sandybrown", "Mixed" = "darkkhaki", "Phytoplankton" = "darkgreen", "Zooplankton" = "saddlebrown")
# Depth
#colors5 <- c("239-244" = "#EFF3FF", "278-337" = "#C6DBEF", "448-455" = "#9ECAE1", "500-900" = "#6BAED6", "814-909" = "#4292C6", "540-1700" = "#2171B5", "500-1800" = "#084594")
colors5 <- c("239-244" = "#EFF3FF", "278-337" = "#C6DBEF", "448-455" = "#9ECAE1", "500-900" = "#6BAED6", "814-909" = "#4292C6", "500-1800" = "#084594")


# PCoA ordination plots
PCoA1 <- ordinate(marmictrans, "PCoA", "unifrac", weighted = TRUE)
p1 <- plot_ordination(marmictrans, PCoA1)

pdf("PCoA_Unifrac_a.pdf")
p1 + scale_fill_manual(values = colors, name="Group") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Group)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Unifrac_b.pdf")
p1 + scale_fill_manual(values = colors2, name="Species") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Unifrac_c.pdf")
p1 + scale_fill_manual(values = colors3, name="Diel Migration") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Migration)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Unifrac_d.pdf")
p1 + scale_fill_manual(values = colors4, name="Diet") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Diet)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Unifrac_e.pdf")
p1 + scale_fill_manual(values = colors5, name="Depth (m)") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()


PCoA2 <- ordinate(marmictrans, "PCoA", "bray")
p2 <- plot_ordination(marmictrans, PCoA2)

pdf("PCoA_Bray_a.pdf")
p2 + scale_fill_manual(values = colors, name="Group") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Group)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_b.pdf")
p2 + scale_fill_manual(values = colors2, name="Species") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_c.pdf")
p2 + scale_fill_manual(values = colors3, name="Diel Migration") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Migration)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_d.pdf")
p2 + scale_fill_manual(values = colors4, name="Diet") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Diet)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_e.pdf")
p2 + scale_fill_manual(values = colors5, name="Depth (m)") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()


# Rarefy taxa for unweighted UniFrac
marmicrf <- rarefy_even_depth(marmic, sample.size = min(sample_sums(marmic)), rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

PCoA3 <- ordinate(marmicrf, "PCoA", "unifrac", weighted = FALSE)
p3 <- plot_ordination(marmicrf, PCoA3)

pdf("PCoA_uwUnifrac_a.pdf")
p3 + scale_fill_manual(values = colors, name="Group") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Group)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_uwUnifrac_b.pdf")
p3 + scale_fill_manual(values = colors2, name="Species") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_uwUnifrac_c.pdf")
p3 + scale_fill_manual(values = colors3, name="Diel Migration") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Migration)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_uwUnifrac_d.pdf")
p3 + scale_fill_manual(values = colors4, name="Diet") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Diet)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_uwUnifrac_e.pdf")
p3 + scale_fill_manual(values = colors5, name="Depth (m)") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()


# Alpha diversity
pdf("Alpha_diversity.pdf")
plot_richness(marmictrans, x="Species", measures=c("Shannon", "Simpson")) + scale_fill_manual(values = colors, name="Group") + theme_bw() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Group)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()

# PERMANOVA
marmicuf <- UniFrac(marmictrans, weighted=TRUE)
marmicbray <- distance(marmictrans, method = "bray")
sampledf <- data.frame(sample_data(marmic))

sink("PERMANOVA.txt")
# Adonis test
adonis(marmicuf ~ Group, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicuf ~ Species, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicuf ~ Diet, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicuf ~ Migration, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicuf ~ Depth, data = sampledf, add = "cailliez", permutations = 999)

adonis(marmicbray ~ Group, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Species, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Diet, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Migration, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Depth, data = sampledf, add = "cailliez", permutations = 999)

# Dispersion test
beta1 <- betadisper(marmicuf, sampledf$Group)
permutest(beta1)

beta2 <- betadisper(marmicuf, sampledf$Species)
permutest(beta2)

beta3 <- betadisper(marmicuf, sampledf$Diet)
permutest(beta3)

beta4 <- betadisper(marmicuf, sampledf$Migration)
permutest(beta4)

beta5 <- betadisper(marmicuf, sampledf$Depth)
permutest(beta5)

beta1b <- betadisper(marmicbray, sampledf$Group)
permutest(beta1b)

beta2b <- betadisper(marmicbray, sampledf$Species)
permutest(beta2b)

beta3b <- betadisper(marmicbray, sampledf$Diet)
permutest(beta3b)

beta4b <- betadisper(marmicbray, sampledf$Migration)
permutest(beta4b)

beta5b <- betadisper(marmicbray, sampledf$Depth)
permutest(beta5b)
sink()

# RDA
pred <- subset(sampledf, select=-c(X.SampleID, BarcodeSequence, LinkerPrimerSequence, Description, Group))

rda1 <- capscale(marmicbray ~ Species + Diet + Migration + Depth, data = pred, na.action = na.exclude, add = "cailliez")
rda2 <- capscale(marmicuf ~ Species + Diet + Migration + Depth, data = pred, na.action = na.exclude, add = "cailliez")

RsquareAdj(rda1)
vif.cca(rda1) #Should be below 10

signif.full1 <- anova.cca(rda1, parallel=getOption("mc.cores"))
signif.full1
signif.axis1 <- anova.cca(rda1, by="axis", parallel=getOption("mc.cores"))
signif.axis1
signif.term1 <- anova.cca(rda1, by="term", parallel=getOption("mc.cores"))
signif.term1
signif.margin1 <- anova.cca(rda1, by="margin", parallel=getOption("mc.cores"))
signif.margin1

RsquareAdj(rda2)
vif.cca(rda2) #Should be below 10

signif.full2 <- anova.cca(rda2, parallel=getOption("mc.cores"))
signif.full2
signif.axis2 <- anova.cca(rda2, by="axis", parallel=getOption("mc.cores"))
signif.axis2
signif.term2 <- anova.cca(rda2, by="term", parallel=getOption("mc.cores"))
signif.term2
signif.margin2 <- anova.cca(rda2, by="margin", parallel=getOption("mc.cores"))
signif.margin2

levels(pred$Species) <- c("Acanthamunnopsis", "Cyclothone", "Krill", "Munneurycope", "Myctophid", "Poeobius", "Sergestes", "Tomopteris", "Vampyroteuthis", "Vitreosalpa")

pdf("Microbiome_RDA_Bray.pdf")
plot(rda1, type="n", scaling=3)
points(rda1, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=colors2[pred$Species])
text(rda1, scaling=3, display="bp", labels = NULL, col="black", cex=0.5)
legend("bottomleft", legend=levels(pred$Species), bty="n", col="gray32", pch=21, cex=1, pt.bg=colors2)
dev.off()

pdf("Microbiome_RDA_Unifrac.pdf")
plot(rda2, type="n", scaling=3)
points(rda2, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=colors2[pred$Species])
text(rda2, scaling=3, display="bp", labels = NULL, col="black", cex=0.5)
legend("topleft", legend=levels(pred$Species), bty="n", col="gray32", pch=21, cex=1, pt.bg=colors2)
dev.off()
