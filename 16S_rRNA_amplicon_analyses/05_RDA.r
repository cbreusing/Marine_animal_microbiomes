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
library(tidyverse)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/microbiome")

biomFile <- import_biom("zotu-table-microbiome_tax.biom", parseFunction = parse_taxonomy_default)
mapFile <- import_qiime_sample_data("Smithsonian_RDA.txt")
biomMapFile <- merge_phyloseq(biomFile, mapFile)

repsetFile <- Biostrings::readDNAStringSet("dna-sequences.fasta")
names(repsetFile) <- gsub("\\s.+$", "", names(repsetFile))

treefile <- read_tree("tree.nwk")
new_tre <- ape::multi2di(treefile)

biomMapTree <- merge_phyloseq(biomMapFile, new_tre)

phyloseq <- merge_phyloseq(biomMapTree, repsetFile)
colnames(tax_table(phyloseq)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

marmic = prune_samples(sample_sums(phyloseq) > 1000, phyloseq)
marmictrans <- transform_sample_counts(marmic, function(x) x/sum(x))

otudf <- t(otu_table(marmictrans))
sampledf <- data.frame(sample_data(marmic))

pdf("Variable_correlation.pdf")
pairs.panels(sampledf[,5:9], scale=T)
dev.off()

pred <- subset(sampledf, select=-c(X.SampleID, BarcodeSequence, LinkerPrimerSequence, Description, Group))

rda <- capscale(otudf ~ Species + Diet + Migration + Depth, data = pred, method = "bray", na.action = na.exclude, add = "cailliez")

RsquareAdj(rda)
vif.cca(rda) #Should be below 10

signif.full <- anova.cca(rda, parallel=getOption("mc.cores"))
signif.full
signif.axis <- anova.cca(rda, by="axis", parallel=getOption("mc.cores"))
signif.axis
signif.term <- anova.cca(rda, by="term", parallel=getOption("mc.cores"))
signif.term
signif.margin <- anova.cca(rda, by="margin", parallel=getOption("mc.cores"))
signif.margin

pdf("Variance_partitions.pdf")
var <- varpart(vegdist(otudf), ~ Species, ~ Diet, ~ Migration, ~ Depth, data = pred, add = "cailliez")
plot(var, cutoff = 0, cex = 0.7, bg=2:5, Xnames=c("Species", "Diet", "Migration", "Depth"))
dev.off()

col <- c("Acanthamunnopsis" = "mistyrose", "Cyclothone" = "lightblue1", "Krill" = "red1", "Munneurycope" = "pink2", "Myctophid" = "royalblue1", "Poeobius" = "lightgoldenrod1", "Sergestes" = "maroon4", "Tomopteris" = "darkgoldenrod1", "Vampyroteuthis" = "midnightblue", "Vitreosalpa" = "cyan4")
levels(pred$Taxon) <- c("Acanthamunnopsis", "Cyclothone", "Krill", "Munneurycope", "Myctophid", "Poeobius", "Sergestes", "Tomopteris", "Vampyroteuthis", "Vitreosalpa")

pdf("Microbiome_RDA.pdf")
plot(rda, type="n", scaling=3)
points(rda, display="species", pch=20, cex=0.7, col="gray32", scaling=3)
points(rda, display="sites", pch=21, cex=1.3, col="gray32", scaling=3, bg=col[pred$Taxon])
text(rda, scaling=3, display="bp", col="black", cex=1)
legend("bottomleft", legend=levels(pred$Taxon), bty="n", col="gray32", pch=21, cex=1, pt.bg=col)
dev.off()

load.rda <- scores(rda, choices = c(1), display = "species")
rownames(load.rda) <- colnames(otudf)

hist(load.rda[,1], main="Loadings on RDA1", xlab="RDA1")

outliers <- function(x,z){
  lims <- mean(x) + c(-1, 1) * z * sd(x)    
  x[x < lims[1] | x > lims[2]] 
}
  
cand1 <- outliers(load.rda[,1],2.5)

ncand <- length(cand1)

cand1 <- cbind.data.frame(rep(1,times=length(cand1)), names(cand1), unname(cand1))

colnames(cand1) <- c("Axis","ASV","Loading")

cand <- rbind(cand1)
cand$ASV <- as.character(cand$ASV)

foo <- matrix(nrow=(ncand), ncol=4)
colnames(foo) <- c("Species", "Diet", "Migration", "Depth")

for (i in 1:length(cand$ASV)) {
  nam <- cand[i,2]
  asv.gen <- otudf[,nam]
  foo[i,] <- apply(pred[, 2:5], 2, function(x) cor(x,asv.gen))
}

cand <- cbind.data.frame(cand,foo)  

length(cand$snp[duplicated(cand$ASV)])
cand <- cand[!duplicated(cand$ASV),]

for (i in 1:length(cand$ASV)) {
  bar <- cand[i,]
  cand[i,8] <- names(which.max(abs(bar[4:7])))
  cand[i,9] <- max(abs(bar[4:7]))
}

colnames(cand)[8] <- "Predictor"
colnames(cand)[9] <- "Correlation"

table(cand$Predictor)

write.table(cand, file="ASV_environment_correlation.txt")
