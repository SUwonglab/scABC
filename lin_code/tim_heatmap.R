library(gplots)
library(Heatplus)
library(vegan)
library(RColorBrewer)
HL60GeneSpecificPromoters5K = scan(file = "HL60GeneSpecificPromoters.txt", what = character())
K562GeneSpecificPromoters5K = scan(file = "K562GeneSpecificPromoters.txt", what = character())
gene_specific_promoters = rownames(normed_log10_combined_expression_matrix)[c(which(rownames(normed_log10_combined_expression_matrix) %in% K562GeneSpecificPromoters5K), which(rownames(normed_log10_combined_expression_matrix) %in% HL60GeneSpecificPromoters5K))]
normed_log10_combined_expression_matrix.gene_specific_promoters = normed_log10_combined_expression_matrix[which(rownames(normed_log10_combined_expression_matrix) %in% gene_specific_promoters), ]
normed_log10_combined_expression_matrix.gene_specific_promoters = normed_log10_combined_expression_matrix.gene_specific_promoters[which(rowSums(normed_log10_combined_expression_matrix.gene_specific_promoters) > 0), ]
colcols = c(rep("red", times = length(grep("K562", colnames(normed_log10_combined_expression_matrix.gene_specific_promoters)))), rep("chartreuse", times = length(grep("HL60", colnames(normed_log10_combined_expression_matrix.gene_specific_promoters)))))
rowcols = rep("black", times = length(rownames(normed_log10_combined_expression_matrix.gene_specific_promoters)))
for(i in 1:length(rowcols)){
  if(rownames(normed_log10_combined_expression_matrix.gene_specific_promoters)[i] %in% K562GeneSpecificPromoters5K){rowcols[i] = "red"}
  else if(rownames(normed_log10_combined_expression_matrix.gene_specific_promoters)[i] %in% HL60GeneSpecificPromoters5K){rowcols[i] = "chartreuse"}
}
scaleyellowred <- colorRampPalette(c("lightyellow", "red"), space = "rgb")(100)
x = normed_log10_combined_expression_matrix.gene_specific_promoters
row.clus = hclust(as.dist(1-cor(t(x), method="spearman")), method="complete")
col.clus = hclust(as.dist(1-cor(x, method="spearman")), method="complete")
cols = c(rep("red", times = length(grep("K562", colnames(normed_log10_combined_expression_matrix)))), rep("chartreuse", times = length(grep("HL60", colnames(normed_log10_combined_expression_matrix)))))
heatmap.2(as.matrix(x), Rowv = as.dendrogram(row.clus), Colv = as.dendrogram(col.clus), col = scaleyellowred, ColSideColors = colcols, RowSideColors = rowcols, margins = c(10, 3),  trace = "none", density.info = "none")