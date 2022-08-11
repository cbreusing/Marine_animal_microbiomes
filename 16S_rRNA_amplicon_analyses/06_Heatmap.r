library(ggplot2)
library(vegan)
library(gridExtra)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

setwd("/Users/Corinna/Documents/Work/Beinart_Lab/Marine_animal_microbiomes_UCSD/16S_amplicons/microbiome")
data1 <- read.table("microbiome-functions-reduced.txt", header=T, row.names=1)
mat1 <- as.matrix(data1)
ann <- read.table("Smithsonian_metadata_reduced.txt", header=T)
cat <- read.table("categories.txt", header=T)
colRamp = colorRamp2(c(0, 0.01), c("lightsteelblue1", "lightsteelblue4"))
pal <- c("gainsboro", "snow4")

ann$Depth <- factor(ann$Depth, levels = c("242", "308", "452", "817-862", "1120-1150"))
ann$Migration <- factor(ann$Migration, levels = c("Yes", "No", "Unknown"))

ha <- HeatmapAnnotation(Migration = ann$Migration, Diet = ann$Diet, Depth = ann$Depth, Host = ann$Species, annotation_height = unit(4, "mm"), col = list(Migration = c("Yes" = "paleturquoise4", "No" = "palevioletred4", "Unknown" = "darkgrey"), Diet = c("Detritus" = "sandybrown", "Phytoplankton" = "darkgreen", "Zooplankton" = "saddlebrown"), Depth = c("242" = "#eff3ff", "308" = "#bdd7e7", "452" = "#6baed6", "817-862" = "#2171b5", "1120-1150" = "#013270"), Host = c("Acanthamunnopsis milleri" = "mistyrose", "Cyclothone atraria" = "lightblue1", "Cyclothone signata" = "lightsteelblue2", "Euphausia pacifica" = "red1", "Munneurycope murrayi" = "pink2", "Stenobrachius leucopsarus" = "royalblue1", "Poeobius meseres" = "lightgoldenrod1", "Eusergestes similis" = "maroon4", "Tomopteris sp." = "darkgoldenrod1", "Vampyroteuthis infernalis" = "midnightblue", "Vitreosalpa sp." = "cyan4")))
hb <- rowAnnotation(df = cat, df2 = anno_mark(at = c(2,5,9,13,18,24,28,38,47,51), labels = c("General Chemoheterotrophy", "Carbohydrate and Derivative Metabolism", "Aromatic Compound, Hydrocarbon and Plastic Degradation", "Methane and Methanol Metabolism", "Phototrophy", "Hydrogen Metabolism", "Sulfur Metabolism", "Nitrogen Metabolism", "Symbiosis", "Other Metabolism"), labels_gp = gpar(fontsize = 9)), show_legend = FALSE, annotation_height = unit(4, "mm"), col = list(Category = c("General Chemoheterotrophy" = pal[1], "Carbohydrate and Derivative Metabolism" = pal[2], "Aromatic Compound, Hydrocarbon and Plastic Degradation" = pal[1], "Methane and Methanol Metabolism" = pal[2], "Phototrophy" = pal[1], "Hydrogen Metabolism" = pal[2], "Sulfur Metabolism" = pal[1], "Nitrogen Metabolism" = pal[2], "Symbiosis" = pal[1], "Other Metabolism" = pal[2])))

pdf("Function_potential_heatmap_full.pdf", width=13, height=9)
map1 <- Heatmap(mat1, width = unit(12, "cm"), height = unit(15, "cm"), right_annotation = hb, top_annotation = ha, km = 1, name = "Abundance", col = colRamp, show_row_names = TRUE, row_names_side = "left", row_names_gp = gpar(fontsize = 7), cluster_rows = FALSE, clustering_distance_columns = "euclidean", clustering_method_columns = "complete", show_column_names = TRUE, column_names_gp = gpar(fontsize = 7, fontface = "bold.italic", col = c(rep("mistyrose", 3), rep("lightblue1", 3), "lightsteelblue2", rep("red1", 4), rep("pink2", 3), rep("royalblue1", 4), rep("lightgoldenrod1", 6), rep("maroon4", 4), rep("darkgoldenrod1", 2), rep("midnightblue", 2), rep("cyan4", 5))))
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
