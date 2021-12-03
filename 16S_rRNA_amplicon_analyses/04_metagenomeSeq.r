library(ggplot2)
library(metagenomeSeq)
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

# Remove low abundance ASVs and all samples with less than 1000 reads
marmic = prune_samples(sample_sums(phyloseq) > 1000, phyloseq)
marmic = filter_taxa(marmic, function(x) sum(x) > 10, TRUE)
marmic = filter_taxa(marmic, function (x) {sum(x > 0) > 1}, TRUE)

# CSS normalization for differential abundance testing
metaseq <- phyloseq_to_metagenomeSeq(marmic)
metaseq = metaseq[, -which(pData(metaseq)$Migration == "Unknown")]
# Trim data
metaseq = filterData(metaseq, present = 1, depth = 1)
metaseqcss = cumNorm(metaseq, p=cumNormStatFast(metaseq))

pd <- pData(metaseqcss)
mod1 <- model.matrix(~1 + Migration, data = pd)
metaseqres1 = fitFeatureModel(metaseqcss, mod1, coef = 2, B = 100)
MRcoefs(metaseqres1, number = 100, adjustMethod = "fdr", alpha = 0.1, eff = 0.5, group = 3, coef = 2, file = "MRcoefs_metaseqres1.txt")
MRfulltable(metaseqres1, number = 100, adjustMethod = "fdr", eff = 0.5, group = 3, coef = 2, file = "MRfulltable_metaseqres1.txt")

# Older model functions for assessing impact of covariates - this was not possible with this dataset due to creation of NaNs in the analysis
migration = pData(metaseqcss)$Migration
group = pData(metaseqcss)$Group
genus = pData(metaseqcss)$Genus
diet = pData(metaseqcss)$Diet
depth = pData(metaseqcss)$Depth
mod2 = model.matrix(~ migration + genus)
settings = zigControl(maxit = 100, verbose = TRUE)
metaseqres2 = fitZig(obj = metaseqcss, mod = mod2, useCSSoffset = TRUE, control = settings)
MRcoefs(metaseqres2, number = 100, adjustMethod = "fdr", alpha = 0.1, eff = 0.5, group = 3, coef = 2, file = "MRcoefs_metaseqres2.txt")
MRfulltable(metaseqres2, number = 100, adjustMethod = "fdr", eff = 0.5, group = 3, coef = 2, file = "MRfulltable_metaseqres2.txt")
