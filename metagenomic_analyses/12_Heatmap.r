library(ggplot2)
library(vegan)
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Marine_animal_microbiomes_UCSD/metagenomics/SHOGUN/utree")
data1 <- read.table("shogun.metabolic.functions.consolidated.txt", header=T, row.names=1)
com <- apply(data1, 2, function(i) i/sum(i, na.rm=TRUE))
mat1 <- as.matrix(com)
ann <- read.table("Smithsonian_metadata_reduced.txt", header=T)
cat <- read.table("categories.txt", header=T)
colRamp = colorRamp2(c(0, 0.01), c("lightsteelblue1", "lightsteelblue4"))
pal <- c("gainsboro", "snow4")

ann$Depth <- factor(ann$Depth, levels = c("242", "308", "817-862", "1120-1150"))
ann$Migration <- factor(ann$Migration, levels = c("Yes", "No", "Unknown"))

ha <- HeatmapAnnotation(Migration = ann$Migration, Diet = ann$Diet, Depth = ann$Depth, Host = ann$Species, annotation_height = unit(4, "mm"), col = list(Migration = c("Yes" = "paleturquoise4", "No" = "palevioletred4", "Unknown" = "darkgrey"), Diet = c("Detritus" = "sandybrown", "Phytoplankton" = "darkgreen", "Zooplankton" = "saddlebrown"), Depth = c("242" = "#eff3ff", "308" = "#bdd7e7", "817-862" = "#2171b5", "1120-1150" = "#013270"), Host = c("Acanthamunnopsis milleri" = "mistyrose", "Cyclothone atraria" = "lightblue1", "Cyclothone signata" = "lightsteelblue2", "Euphausia pacifica" = "red1", "Munneurycope murrayi" = "pink2", "Stenobrachius leucopsarus" = "royalblue1", "Poeobius meseres" = "lightgoldenrod1", "Eusergestes similis" = "maroon4", "Tomopteris sp." = "darkgoldenrod1", "Vampyroteuthis infernalis" = "midnightblue")))
hb <- rowAnnotation(df = cat, df2 = anno_mark(at = c(11,30,39,44,53,65,72,76,79,84), labels = c("Carbohydrate and Derivative Metabolism", "Various Organic Compound Degradation", "Methane, Methanol and Other Carbon Metabolism", "Protein and Amino Acid Metabolism", "Nitrogen Metabolism", "Sulfur Metabolism", "Hydrogen Metabolism", "Other Respiratory and Energy Metabolism", "Fatty Acid, Lipid and Steroid Metabolism", "Cofactor, Secondary and Other Metabolism"), labels_gp = gpar(fontsize = 9)), show_legend = FALSE, annotation_height = unit(4, "mm"), col = list(Category = c("Carbohydrate and Derivative Metabolism" = pal[1], "Various Organic Compound Degradation" = pal[2], "Methane, Methanol and Other Carbon Metabolism" = pal[1], "Protein and Amino Acid Metabolism" = pal[2], "Nitrogen Metabolism" = pal[1], "Sulfur Metabolism" = pal[2], "Hydrogen Metabolism" = pal[1], "Other Respiratory and Energy Metabolism" = pal[2], "Fatty Acid, Lipid and Steroid Metabolism" = pal[1], "Cofactor, Secondary and Other Metabolism" = pal[2])))

pdf("Function_potential_heatmap_full.pdf", width=13, height=12)
map1 <- Heatmap(mat1, width = unit(12, "cm"), height = unit(24, "cm"), right_annotation = hb, top_annotation = ha, km = 1, name = "Abundance", col = colRamp, show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 7), cluster_rows = FALSE, clustering_distance_columns = "euclidean", clustering_method_columns = "complete", show_column_names = TRUE, column_names_gp = gpar(fontsize = 7, fontface = "bold.italic", col = c(rep("mistyrose", 4), rep("lightblue1", 3), "lightsteelblue2", rep("red1", 4), rep("pink2", 3), rep("royalblue1", 4), rep("lightgoldenrod1", 5), rep("maroon4", 4), rep("darkgoldenrod1", 2), rep("midnightblue", 2))))
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
