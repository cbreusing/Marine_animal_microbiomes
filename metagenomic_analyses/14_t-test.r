setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Marine_animal_microbiomes_UCSD/metagenomics/SHOGUN/utree")

data <- read.table("ratio-aerobic-anaerobic.txt", header=T)

shapiro.test(data$Anaerobic)

wilcox.test(Anaerobic ~ Group, data = data, paired = FALSE)
