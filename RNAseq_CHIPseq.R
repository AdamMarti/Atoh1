# ######################################
#  Comparing CHIPseq x RNAseq genes
# ############ ############ ############ 
library(dplyr)
library(ggplot2)
library(reshape2)
library(pheatmap)

setwd('/home/annie/Documents/PGB/Project/CHIPseq_RNAseq')
# Read the CSV files
enriched_genes <- read.csv("enriched_genes_with_ensembl.csv")
deg_results <- read.csv("DEG_results_filtered.csv")
cpm.matrx <- read.csv("cpm_matrix.csv", row.names = 1)

# Desired order of columns (WT first, then KO)
desired_order <- c("WT1", "WT2", "KO1", "KO2")

# Reorder columns in the CPM matrix to match the desired order
cpm.matrx <- cpm.matrx[, match(desired_order, colnames(cpm.matrx))]

# Remove duplicates by keeping only the first occurrence of each Ensembl_ID
enriched_genes_unique <- enriched_genes %>%
  distinct(Ensembl_ID, .keep_all = TRUE)

deg_results <- deg_results %>%
  mutate(Ensembl_ID = Geneid) %>%  # Create a new column 'Ensembl_ID' with the same values as 'Geneid'
  select(-Geneid)  # Optionally, remove the old 'Geneid' column if no longer needed

# First, filter the DEG results to keep only the Ensembl_IDs that are present in the enriched_genes dataset
deg_results_matched <- deg_results %>%
  filter(Ensembl_ID %in% enriched_genes_unique$Ensembl_ID)

# Merge the filtered DEG results with the enriched_genes dataset (matching by Ensembl_ID)
final_merged_data <- deg_results_matched %>%
  left_join(enriched_genes_unique %>% select(Ensembl_ID, Description, geneID), by = "Ensembl_ID")

# Check the result
head(final_merged_data)

# Save the final merged dataset (DEGs with matches) to a new CSV file
write.csv(final_merged_data, "final_DEG_with_enriched_genes_matched.csv", row.names = FALSE)

cat("Final DEG file with enriched genes matches saved as 'final_DEG_with_enriched_genes_matched.csv'\n")

# Assuming final_merged_data has the 'Ensembl_ID' and 'GeneSymbol' columns
# Extract the gene symbols from final_merged_data for the 45 selected genes
gene_symbols <- final_merged_data$GeneSymbol
selected_genes <- final_merged_data$Ensembl_ID
selected_cpm <- cpm.matrx[rownames(cpm.matrx) %in% selected_genes, ]

# Replace row names in selected_cpm with GeneSymbol
rownames(selected_cpm) <- gene_symbols[match(rownames(selected_cpm), selected_genes)]

# Create the condition vector (e.g., WT and KO for the samples)
condition <- c("WT", "WT", "KO", "KO")  # Adjust this based on your actual sample conditions

# Create annotation for the columns (condition for each sample)
ha_col <- data.frame(Condition = condition, row.names = colnames(selected_cpm))

# Define color scale for the heatmap
heatcol <- colorRampPalette(c("blue", "white", "red"), space = "rgb")

# Create the heatmap
pheatmap(
  selected_cpm,                        # The matrix of CPM values for the selected genes
  col = heatcol(256),                   # Color palette from blue to red
  cluster_rows = TRUE,                  # Cluster genes based on expression similarity
  cluster_cols = TRUE,                  # Cluster samples (conditions)
  scale = "row",                        # Scale by row (gene), so each gene is centered
  main = "Heatmap of Gene Expression (CPM) for Top 45 Genes",  # Title of the heatmap
  cexRow = 0.8,                         # Size of row labels (gene names)
  cexCol = 0.8,                         # Size of column labels (sample names)
  show_rownames = TRUE,                 # Display gene names
  show_colnames = TRUE,                 # Display sample names
  treeheight_row = 0,                   # Remove row dendrogram height
  annotation_col = ha_col               # Add column annotations (conditions)
)

#Find common hits
# 1. Sort the deg_results by FDR and get the top 50 genes
top_45_deg <- deg_results %>%
  arrange(FDR) %>%  # Sort by FDR (ascending)
  head(45)           # Get top 50

# 2. Extract the gene symbols from the top 50 DEGs and final_merged_data
top_45_deg_symbols <- top_45_deg$GeneSymbol
final_merged_symbols <- final_merged_data$GeneSymbol

# 3. Find matching genes between the top 50 DEGs and final_merged_data
matching_genes <- intersect(top_45_deg_symbols, final_merged_symbols)

# 4. Get the data for the matching genes from final_merged_data
matching_genes_data <- final_merged_data %>%
  filter(GeneSymbol %in% matching_genes)

# 5. Print the matching genes data
print(matching_genes_data)

# Optionally, save the matching genes data to a file
write.csv(matching_genes_data, "matching_genes_from_top_45_DEGs.csv", row.names = FALSE)

cat("Matching genes have been saved to 'matching_genes_from_top_45_DEGs.csv'\n")

# Ensure the 'final_merged_data' dataframe contains 'GeneSymbol' and 'logFC'
# Extract GeneSymbol and logFC columns for the plot
logFC_final_data <- final_merged_data %>%
  select(GeneSymbol, logFC) %>%
  arrange(desc(logFC))  # Optionally, order by logFC values

# Create a bar plot of logFC
ggplot(logFC_final_data, aes(x = reorder(GeneSymbol, logFC), y = logFC, fill = logFC)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip coordinates for better visibility
  labs(title = "Log Fold Change (logFC) of Final Merged Genes",
       x = "Gene Symbol",
       y = "Log Fold Change") +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
