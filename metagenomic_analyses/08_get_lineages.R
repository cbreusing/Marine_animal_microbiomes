library(taxize)
library(myTAI)
library(plyr)
library(dplyr)

# Extract lineage information from TaxIDs
setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Marine_animal_microbiomes_UCSD/metagenomics/Acanthamunnopsis/eukaryotic_fraction") # Acanthamunnopsis, Cyclothone, Krill, Munneurycope, Myctophid, Poeobius, Sergestes

Sys.setenv(ENTREZ_KEY = "063e59da58f82e1a608dc8e31789e4da8c09")

taxids <- read.table("Acanthamunnopsis_all.taxIDs.txt", header=F)

taxid_list <- list()
for (i in 1:nrow(taxids)){
tax  <- classification(as.uid(taxids[i,]), db = "ncbi", return_id = FALSE, max_tries = 10)
unname(tax)
df <- data.frame(lapply(tax[[1]]$name, function(x) data.frame(x)))
colnames(df) <- tax[[1]]$rank
df$taxid <- taxids[i,]
taxid_list[[i]] <- df
Sys.sleep(0.5)
}

taxid_df <- do.call(rbind.fill, taxid_list)

taxdata <- taxid_df %>% select(c('taxid', 'superkingdom', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species'))

write.table(taxdata, file="Acanthamunnopsis.taxLineages.txt")



