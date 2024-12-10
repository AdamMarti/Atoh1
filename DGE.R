# ######################################
# Load Libraries
# ######################################
library(GenomicRanges)
library(dplyr)
library(edgeR)
library(ggplot2)
library(matrixStats)
library(reshape2)
library(ggpubr)
library(pheatmap)
library(org.Mm.eg.db)

# ######################################
# Load and Prepare Data
# ######################################

# Set working directory
setwd("/home/annie/Documents/PGB/Project/RNAseq/03_counts")

# List and read count files
count_files <- list.files(pattern = "*.counts.txt")
counts_list <- lapply(count_files, function(file) {
  counts_data <- read.table(file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  sample_name <- gsub("\\.SRR[0-9]+\\.counts\\.txt$", "", file)  # Extract sample name
  colnames(counts_data)[2] <- sample_name
  counts_data
})

# Merge counts by Geneid
combined_counts <- Reduce(function(x, y) merge(x, y, by = "Geneid", all = TRUE), counts_list)
counts <- combined_counts  # Save merged counts
counts.m <- as.matrix(counts[, -1])  # Remove Geneid for matrix conversion
rownames(counts.m) <- counts$Geneid  # Set Geneid as row names

# ######################################
# Normalization
# ######################################

# Create DGEList object
d <- DGEList(counts = counts.m)

# Normalize counts using TMM
d <- calcNormFactors(d, method = "TMM")

# Extract CPM values
cpm.matrx <- cpm(d, log = FALSE)

# Save normalized CPM matrix
write.csv(cpm.matrx, "cpm_matrix.csv", row.names = TRUE)

# ######################################
# PCA Plot
# ######################################

# Define conditions and reorder columns
condition <- factor(c("WT", "WT", "KO", "KO"), levels = c("WT", "KO"))
current_order <- colnames(cpm.matrx)
desired_order <- c("WT1", "WT2", "KO1", "KO2")
cpm.matrx <- cpm.matrx[, match(desired_order, current_order)]

# Compute PCA
top100 <- head(order(rowVars(cpm.matrx), decreasing = TRUE), 100)
pca <- prcomp(t(cpm.matrx[top100, ]), center = TRUE, scale = TRUE)

# Generate PCA plot
ggplot(data = as.data.frame(pca$x), aes(x = PC1, y = PC2, color = factor(condition))) +
  geom_point(alpha = 0.8, size = 4) +
  geom_label_repel(aes(label = colnames(cpm.matrx)), box.padding = 0.2) +
  xlab(paste0("PC1 (", round(summary(pca)$importance[2, "PC1"] * 100), "% variance)")) +
  ylab(paste0("PC2 (", round(summary(pca)$importance[2, "PC2"] * 100), "% variance)")) +
  theme_minimal() +
  scale_color_manual(values = c("WT" = "blue", "KO" = "red")) +
  ggtitle("PCA Plot of RNA-Seq Data")

# ######################################
# Differential Expression Analysis
# ######################################

# Create DGEList and estimate dispersion
d <- DGEList(counts = counts.m, group = condition)
d <- calcNormFactors(d)
d <- estimateDisp(d)

# Fit GLM model and perform differential expression test
fit <- glmQLFit(d)
result <- glmQLFTest(fit)

# Process results
result_table <- as.data.frame(result$table)
result_table$FDR <- p.adjust(result_table$PValue, method = "BH")
write.table(result_table[order(result_table$FDR), ], file = "DEG_results_sorted_by_FDR.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# Map Ensembl IDs to Gene Symbols
result_table_sorted <- result_table[order(result_table$FDR), ]
ensembl_ids <- rownames(result_table_sorted)
gene_symbols <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL")
result_table_sorted$GeneSymbol <- gene_symbols
result_table_sorted <- result_table_sorted[, c("GeneSymbol", setdiff(names(result_table_sorted), "GeneSymbol"))]

# Save results
write.table(result_table_sorted, "significant_DEGs_with_GeneSymbols_FDR_0.05_logFC_1.txt", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)

# ######################################
# Heatmap for Top 50 Genes
# ######################################

# Select top 50 genes
top_50_genes <- head(result_table_sorted, 50)
selected_matrix <- cpm.matrx[rownames(top_50_genes), ]
rownames(selected_matrix) <- top_50_genes$GeneSymbol
ha_col <- data.frame(Cond = condition, row.names = colnames(selected_matrix))

# Plot heatmap
pheatmap(selected_matrix, 
         col = colorRampPalette(c("blue", "white", "red"))(256),
         cluster_rows = TRUE,
         clustering_method = "ward.D2",
         clustering_distance_cols = "euclidean",
         scale = "row",
         main = "Top 50 Genes Expression Heatmap",
         annotation_col = ha_col)

# ######################################
# Volcano Plot for DEGs
# ######################################

# Create volcano plot data frame
volcano_df <- data.frame(
  logFC = result$table$logFC,
  negLogP = -log10(result$table$PValue),
  FDR = p.adjust(result$table$PValue, method = "BH"),
  Geneid = rownames(result$table)
)
volcano_df$Significance <- ifelse(volcano_df$FDR < 0.05 & abs(volcano_df$logFC) > 1, "Significant", "Not Significant")

# Generate volcano plot
ggplot(volcano_df, aes(x = logFC, y = negLogP)) +
  geom_point(aes(color = Significance), alpha = 0.8, size = 2) +
  scale_color_manual(values = c("Significant" = "red", "Not Significant" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "blue") +
  labs(title = "Volcano Plot of Differential Expression",
       x = "Log2 Fold Change",
       y = "-Log10 P-Value",
       color = "Significance") +
  theme_minimal() +
  theme(legend.position = "top")
