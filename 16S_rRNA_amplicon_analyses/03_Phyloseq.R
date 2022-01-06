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

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/microbiome")

# Import final ASV biom and mapping files
biomFile <- import_biom("zotu-table-microbiome_tax.biom", parseFunction = parse_taxonomy_default)
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
colnames(tax_table(phyloseq)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

# Remove low abundance ASVs and all samples with less than 1000 reads (recommended)
marmic = prune_samples(sample_sums(phyloseq) > 1000, phyloseq)
marmic = filter_taxa(marmic, function(x) sum(x) > 10, TRUE)
marmic = filter_taxa(marmic, function (x) {sum(x > 0) > 1}, TRUE)

sample_data(marmic)$Depth <- factor(sample_data(marmic)$Depth, levels = c("239-244", "278-337", "448-455", "500-900", "814-909", "540-1700", "500-1800"))
sample_data(marmic)$Migration <- factor(sample_data(marmic)$Migration, levels = c("Yes", "No", "Unknown"))

# Transform counts to proportional abundances, this has recently been shown to be the most appropriate normalization technique for ecological microbiome analyses such as alpha and beta diversity (when using abundance based distance methods such as Bray-Curtis or weighted UniFrac)
marmictrans <- transform_sample_counts(marmic, function(x) x/sum(x))

otu <- t(otu_table(marmic))
spec <- specaccum(otu)

pdf("Specaccum.pdf")
plot(spec, ci.type="poly", col = "black", lwd=2, ci.lty=0, ci.col="grey", ci = 1.96, xlab = "# individuals", ylab = "# ASVs")
dev.off()

pdf("Fractional_abundance_a.pdf", width=20, height=14)
plot_bar(marmictrans, x= "Sample", fill = "Phylum") + facet_grid(. ~ Group, scales = "free", space = "free") + geom_bar(aes(color=Phylum, fill=Phylum), stat='identity', position='stack') + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_text(size = 15, face="bold")) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom") + scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"))
dev.off()

pdf("Fractional_abundance_b.pdf", width=20, height=14)
plot_bar(marmictrans, x= "Sample", fill = "Class") + facet_grid(. ~ Group, scales = "free", space = "free") + geom_bar(aes(color=Class, fill=Class), stat='identity', position='stack') + ylab("Fractional abundance") + theme_bw() + theme(axis.title = element_text(size=15, face="bold")) + theme(axis.text = element_text(size=13)) + theme(legend.text = element_text(size = 13)) + theme(legend.title = element_text(size = 15, face="bold")) + theme(strip.text.x = element_text(size = 15, face="bold")) + theme(legend.key = element_rect(size = 0.4)) + theme(legend.position="bottom") + scale_color_viridis(discrete = TRUE) + scale_fill_viridis(discrete = TRUE) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, face = "italic"))
dev.off()

# Group
colors <- c("Fish" = "cornflowerblue", "Isopod" = "pink1", "Krill" = "red1", "Polychaete" = "goldenrod1", "Shrimp" = "maroon4", "Squid" = "midnightblue", "Tunicate" = "cyan4")
# Species
colors2 <- c("Acanthamunnopsis milleri" = "mistyrose", "Cyclothone atraria" = "lightblue1", "Cyclothone signata" = "lightsteelblue2", "Munneurycope murrayi" = "pink2", "Poeobius meseres" = "lightgoldenrod1", "Eusergestes similis" = "maroon4", "Stenobrachius leucopsarus" = "royalblue1", "Euphausia pacifica" = "red1", "Tomopteris sp." = "darkgoldenrod1", "Vampyroteuthis infernalis" = "midnightblue", "Vitreosalpa sp." = "cyan4")
# Migration
colors3 <- c("Yes" = "paleturquoise4", "No" = "palevioletred4", "Unknown" = "darkgrey")
# Diet
colors4 <- c("Detritus" = "sandybrown", "Phytoplankton" = "darkgreen", "Zooplankton" = "saddlebrown")
# Depth
colors5 <- c("239-244" = "#EFF3FF", "278-337" = "#C6DBEF", "448-455" = "#9ECAE1", "500-900" = "#6BAED6", "814-909" = "#4292C6", "540-1700" = "#2171B5", "500-1800" = "#084594")


# PCoA ordination plots
PCoA1 <- ordinate(marmictrans, "PCoA", "unifrac", weighted = TRUE)
p1 <- plot_ordination(marmictrans, PCoA1)

pdf("PCoA_Unifrac_a.pdf")
p1 + scale_fill_manual(values = colors, name="Group") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Group)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Unifrac_b.pdf")
p1 + scale_fill_manual(values = colors2, name="Species") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + theme(legend.text = element_text(face = "italic"))
dev.off()

pdf("PCoA_Unifrac_c.pdf")
p1 + scale_fill_manual(values = colors3, name="Diel Vertical Migration") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Migration)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Unifrac_d.pdf")
p1 + scale_fill_manual(values = colors4, name="Diet") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Diet)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Unifrac_e.pdf")
p1 + scale_fill_manual(values = colors5, name="Depth (m)") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Unifrac_f.pdf")
p1 + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species), shape = Migration, color = factor(Species)), size=3) + scale_shape_manual(values = c(19, 17, 15), name="Diel Vertical Migration") + scale_fill_manual(values = colors2, name="Species", aesthetics = c("colour", "fill")) + theme(text = element_text(size = 15)) + theme(legend.text = element_text(face = "italic"))
dev.off()

pdf("PCoA_Unifrac_g.pdf")
p1 + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth), shape = Diet, color = factor(Depth)), size=3) + scale_shape_manual(values = c(19, 17, 15), name="Diet") + scale_fill_manual(values = colors5, name="Depth (m)", aesthetics = c("colour", "fill")) + theme(text = element_text(size = 15))
dev.off()


PCoA2 <- ordinate(marmictrans, "PCoA", "bray")
p2 <- plot_ordination(marmictrans, PCoA2)

pdf("PCoA_Bray_a.pdf")
p2 + scale_fill_manual(values = colors, name="Group") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Group)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_b.pdf")
p2 + scale_fill_manual(values = colors2, name="Species") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + theme(legend.text = element_text(face = "italic"))
dev.off()

pdf("PCoA_Bray_c.pdf")
p2 + scale_fill_manual(values = colors3, name="Diel Vertical Migration") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Migration)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_d.pdf")
p2 + scale_fill_manual(values = colors4, name="Diet") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Diet)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_e.pdf")
p2 + scale_fill_manual(values = colors5, name="Depth (m)") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_Bray_f.pdf")
p2 + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species), shape = Migration, color = factor(Species)), size=3) + scale_shape_manual(values = c(19, 17, 15), name="Diel Vertical Migration") + scale_fill_manual(values = colors2, name="Species", aesthetics = c("colour", "fill")) + theme(text = element_text(size = 15)) + theme(legend.text = element_text(face = "italic"))
dev.off()

pdf("PCoA_Bray_g.pdf")
p2 + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth), shape = Diet, color = factor(Depth)), size=3) + scale_shape_manual(values = c(19, 17, 15), name="Diet") + scale_fill_manual(values = colors5, name="Depth (m)", aesthetics = c("colour", "fill")) + theme(text = element_text(size = 15))
dev.off()

# Rarefy taxa for unweighted UniFrac
marmicrf <- rarefy_even_depth(marmic, sample.size = min(sample_sums(marmic)), rngseed = TRUE, replace = TRUE, trimOTUs = TRUE, verbose = TRUE)

PCoA3 <- ordinate(marmicrf, "PCoA", "unifrac", weighted = FALSE)
p3 <- plot_ordination(marmicrf, PCoA3)

pdf("PCoA_uwUnifrac_a.pdf")
p3 + scale_fill_manual(values = colors, name="Group") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Group)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_uwUnifrac_b.pdf")
p3 + scale_fill_manual(values = colors2, name="Species") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15)) + theme(legend.text = element_text(face = "italic"))
dev.off()

pdf("PCoA_uwUnifrac_c.pdf")
p3 + scale_fill_manual(values = colors3, name="Diel Vertical Migration") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Migration)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_uwUnifrac_d.pdf")
p3 + scale_fill_manual(values = colors4, name="Diet") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Diet)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_uwUnifrac_e.pdf")
p3 + scale_fill_manual(values = colors5, name="Depth (m)") + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth)), shape=21, color="gray32", size=3) + theme(text = element_text(size = 15))
dev.off()

pdf("PCoA_uwUnifrac_f.pdf")
p3 + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Species), shape = Migration, color = factor(Species)), size=3) + scale_shape_manual(values = c(19, 17, 15), name="Diel Vertical Migration") + scale_fill_manual(values = colors2, name="Species", aesthetics = c("colour", "fill")) + theme(text = element_text(size = 15)) + theme(legend.text = element_text(face = "italic"))
dev.off()

pdf("PCoA_uwUnifrac_g.pdf")
p3 + theme_classic() + geom_point(position=position_jitter(width=0, height=0), aes(fill = factor(Depth), shape = Diet, color = factor(Depth)), size=3) + scale_shape_manual(values = c(19, 17, 15), name="Diet") + scale_fill_manual(values = colors5, name="Depth (m)", aesthetics = c("colour", "fill")) + theme(text = element_text(size = 15))
dev.off()


# Traditional PERMANOVA
marmicuf <- UniFrac(marmictrans, weighted=TRUE)
marmicbray <- distance(marmictrans, method = "bray")
sampledf <- data.frame(sample_data(marmic))

sink("PERMANOVA.txt")
# Adonis test
adonis(marmicuf ~ Species, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicuf ~ Diet, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicuf ~ Migration, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicuf ~ Depth, data = sampledf, add = "cailliez", permutations = 999)

adonis(marmicbray ~ Species, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Diet, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Migration, data = sampledf, add = "cailliez", permutations = 999)
adonis(marmicbray ~ Depth, data = sampledf, add = "cailliez", permutations = 999)

# Dispersion test
beta1 <- betadisper(marmicuf, sampledf$Species)
permutest(beta1)

beta2 <- betadisper(marmicuf, sampledf$Diet)
permutest(beta2)

beta3 <- betadisper(marmicuf, sampledf$Migration)
permutest(beta3)

beta4 <- betadisper(marmicuf, sampledf$Depth)
permutest(beta4)

beta1b <- betadisper(marmicbray, sampledf$Species)
permutest(beta1b)

beta2b <- betadisper(marmicbray, sampledf$Diet)
permutest(beta2b)

beta3b <- betadisper(marmicbray, sampledf$Migration)
permutest(beta3b)

beta4b <- betadisper(marmicbray, sampledf$Depth)
permutest(beta4b)
sink()

