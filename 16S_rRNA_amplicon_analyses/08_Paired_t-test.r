setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/microbiome")

data <- read.table("Actinobacteria_mean_proportions.txt", header=T)

d <- with(data, Value[Group == "Metagenomics"] - Value[Group == "Amplicons"])
shapiro.test(d)

res <- t.test(Value ~ Group, data = data, paired = TRUE)
res