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

microbiometrans <- transform_sample_counts(phyloseq, function(x) x/sum(x))

microbiome.dist <- UniFrac(microbiometrans, weighted=TRUE)


diet <- import_biom("zotu-table-diet_tax.biom", parseFunction = parse_taxonomy_default)
diettrans <- transform_sample_counts(diet, function(x) x/sum(x))

diet.dist <- distance(diettrans, method="bray")

mantel(microbiome.dist, diet.dist, method="spearman", permutations=9999)
