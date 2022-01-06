library(ggplot2)
library(vegan)
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

setwd("/Users/Corinna/Documents/PostDoc/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/microbiome")
data1 <- read.table("microbiome-functions-reduced.txt", header=T, row.names=1)
mat1 <- as.matrix(data1)
ann <- read.table("Smithsonian_metadata_reduced.txt", header=T)
colRamp = colorRamp2(c(0, 0.01), c("lightsteelblue1", "lightsteelblue4"))

ann$Depth <- factor(ann$Depth, levels = c("239-244", "278-337", "448-455", "500-900", "814-909", "540-1700", "500-1800"))
ann$Migration <- factor(ann$Migration, levels = c("Yes", "No", "Unknown"))

ha <- HeatmapAnnotation(Migration = ann$Migration, Diet = ann$Diet, Depth = ann$Depth, annotation_height = unit(4, "mm"), col = list(Migration = c("Yes" = "paleturquoise4", "No" = "palevioletred4", "Unknown" = "darkgrey"), Diet = c("Detritus" = "sandybrown", "Phytoplankton" = "darkgreen", "Zooplankton" = "saddlebrown"), Depth = c("239-244" = "#EFF3FF", "278-337" = "#C6DBEF", "448-455" = "#9ECAE1", "500-900" = "#6BAED6", "814-909" = "#4292C6", "540-1700" = "#2171B5", "500-1800" = "#084594")))

pdf("Function_potential_heatmap_full.pdf", width=9, height=8)
map1 <- Heatmap(mat1, width = unit(12, "cm"), height = unit(15, "cm"), top_annotation = ha, km = 1, name = "Abundance", col = colRamp, show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 7), cluster_rows = FALSE, clustering_distance_columns = "euclidean", clustering_method_columns = "complete", show_column_names = TRUE, column_names_gp = gpar(fontsize = 7, fontface = "bold.italic", col = c(rep("mistyrose", 3), rep("lightblue1", 3), "lightsteelblue2", rep("red1", 4), rep("pink2", 3), rep("royalblue1", 4), rep("lightgoldenrod1", 6), rep("maroon4", 4), rep("darkgoldenrod1", 2), rep("midnightblue", 2), rep("cyan4", 5))))
draw(map1, annotation_legend_side = "right", merge_legend = TRUE)
dev.off()

mat1bray <- vegdist(t(mat1), method = "bray")

# Traditional PERMANOVA
sink("PERMANOVA_functional_potential.txt")
# Adonis test
adonis(mat1bray ~ Species, data = ann, add = "cailliez", permutations = 999)
adonis(mat1bray ~ Diet, data = ann, add = "cailliez", permutations = 999)
adonis(mat1bray ~ Migration, data = ann, add = "cailliez", permutations = 999)
adonis(mat1bray ~ Depth, data = ann, add = "cailliez", permutations = 999)

# Dispersion test
beta1 <- betadisper(mat1bray, ann$Species)
permutest(beta1)

beta2 <- betadisper(mat1bray, ann$Diet)
permutest(beta2)

beta3 <- betadisper(mat1bray, ann$Migration)
permutest(beta3)

beta4 <- betadisper(mat1bray, ann$Depth)
permutest(beta4)
sink()
