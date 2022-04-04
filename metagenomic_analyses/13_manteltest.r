library(vegan)
library(ggplot2)
library(tidyverse)
library(ape)
library(poppr)
library(pegas)
library(adegenet)
library(ade4)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Marine_animal_microbiomes_UCSD/metagenomics/SHOGUN/utree")

amplicon <- read.table("zotu-table-MT1.txt", header=T, row.names=1)
amplicon.mat <- t(apply(amplicon, 2, function(i) i/sum(i, na.rm=TRUE)))
amplicon.dist <- vegdist(amplicon.mat, method="bray", binary=F, na.rm=TRUE)

shotgun <- read.table("otu-table-MT1.txt", header=T, row.names=1)
shotgun.mat <- t(apply(shotgun, 2, function(i) i/sum(i, na.rm=TRUE)))
shotgun.dist <- vegdist(shotgun.mat, method="bray", binary=F, na.rm=TRUE)

sink("MantelTest1.txt")
mantel(amplicon.dist, shotgun.dist, permutations = 999, na.rm = TRUE, method="spearman", parallel = getOption("mc.cores"))
sink()
