library(ggplot2)
library(plyr)
library(dplyr)
library(devtools)
library(phyloseq)
library(ape)
library("biomformat")
library(vegan)

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/Mantel_Test")

microbiome <- import_biom("zotu-table-microbiome_tax.biom", parseFunction = parse_taxonomy_default)
mapFile <- import_qiime_sample_data("Smithsonian_metadata.txt")
biomMapFile <- merge_phyloseq(microbiome, mapFile)
repsetFile <- Biostrings::readDNAStringSet("dna-sequences.fasta")
names(repsetFile) <- gsub("\\s.+$", "", names(repsetFile))
treefile <- read_tree("tree.nwk")
new_tre <- ape::multi2di(treefile)
biomMapTree <- merge_phyloseq(biomMapFile, new_tre)
phyloseq <- merge_phyloseq(biomMapTree, repsetFile)
colnames(tax_table(phyloseq)) = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phyloseq = prune_samples(sample_sums(phyloseq) > 1000, phyloseq)
phyloseq = filter_taxa(phyloseq, function(x) sum(x) > 10, TRUE)
phyloseq = subset_samples(phyloseq, Description != "Acanthamunnopsis1" & Description != "Acanthamunnopsis2" & Description != "Acanthamunnopsis3" & Description != "Cyclothone2" & Description != "Cyclothone4" & Description != "Krill1" & Description != "Krill4" & Description != "Myctophid2" & Description != "Poeobius1" & Description != "Poeobius3" & Description != "Poeobius4" & Description != "Poeobius5" & Description != "Tomopteris1" & Description != "Tomopteris2" & Description != "Vampyroteuthis1" & Description != "Vampyroteuthis2" & Description != "Vitreosalpa2" & Description != "Vitreosalpa3")
phyloseq = filter_taxa(phyloseq, function(x) sum(x) > 0, TRUE)

microbiometrans <- transform_sample_counts(phyloseq, function(x) x/sum(x))

microbiome.dist <- UniFrac(microbiometrans, weighted=TRUE)
microbiome.dist2 <- distance(microbiometrans, method = "bray")

diet <- import_biom("zotu-table-diet_tax.biom", parseFunction = parse_taxonomy_default)
biomMapFile2 <- merge_phyloseq(diet, mapFile)
repsetFile2 <- Biostrings::readDNAStringSet("dna-sequences2.fasta")
names(repsetFile2) <- gsub("\\s.+$", "", names(repsetFile2))
biomMapTree2 <- merge_phyloseq(biomMapFile2, new_tre)
phyloseq2 <- merge_phyloseq(biomMapTree2, repsetFile2)
colnames(tax_table(phyloseq2)) = c("Domain", "Category", "Phylum", "Class", "Order", "Family", "Genus", "Species")

phyloseq2 = prune_samples(sample_sums(phyloseq2) > 500, phyloseq2)
phyloseq2 = filter_taxa(phyloseq2, function(x) sum(x) > 10, TRUE)
diettrans <- transform_sample_counts(phyloseq2, function(x) x/sum(x))

diet.dist <- UniFrac(diettrans, weighted=TRUE)
diet.dist2 <- distance(diettrans, method="bray")

mantel(microbiome.dist, diet.dist, method="spearman", permutations=9999)
mantel(microbiome.dist2, diet.dist2, method="spearman", permutations=9999)
mantel(microbiome.dist, diet.dist2, method="spearman", permutations=9999)

