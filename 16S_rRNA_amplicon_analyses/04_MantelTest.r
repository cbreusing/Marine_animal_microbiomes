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
microbiometrans <- transform_sample_counts(microbiome, function(x) x/sum(x))

microbiome.dist <- distance(microbiometrans, method="bray")


diet <- import_biom("zotu-table-diet_tax.biom", parseFunction = parse_taxonomy_default)
diettrans <- transform_sample_counts(diet, function(x) x/sum(x))

diet.dist <- distance(diettrans, method="bray")

mantel(microbiome.dist, diet.dist, method="spearman", permutations=9999)
