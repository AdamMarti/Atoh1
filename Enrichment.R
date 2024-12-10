
# ######################################
#  Setting up environment
# ############ ############ ############ 

#Set working directory
setwd('/home/annie/Documents/PGB/Project/Pathway_Enrichment')

out_path <- "/home/annie/Documents/PGB/Project/Pathway_Enrichment/Results" 

#Load libraries
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)
library(DOSE)
library(enrichplot)
library(ggplot2)

######################################
#  Load DEG Data
# ############ ############ ############

#Read DGE results
df <- read.delim('significant_DEGs_with_GeneSymbols_FDR_0.05_logFC_1.txt')
colnames(df)[1] <- 'Geneid' #Rename as Gene Id the first column
head(df)

# Sort the data frame by FDR (ascending)
df <- df %>% arrange(FDR)

# View the sorted data
head(df)

#Add a column for up-regulation or down-regulation
# Add the 'diffexpressed' column based on logFC and FDR
df <- df %>% mutate(diffexpressed = case_when(
  logFC > 0 & FDR < 0.05 ~ 'UP',
  logFC < 0 & FDR < 0.05 ~ 'DOWN',
  FDR > 0.05 ~ 'NO'
))

# Remove rows where diffexpressed is "NO", FDR<0.05
df <- df %>% filter(diffexpressed != "NO")


# Save to a new CSV file
write.csv(df, "DEG_results_filtered.csv", row.names = FALSE)

##At this point we have all DEG that are up or down-regulated with p<0.05

# ######################################
#  Prepare Background Genes
# ############ ############ ############

# Get a list of all mouse gene ENSEMBL IDs
background_genes <- keys(org.Mm.eg.db, keytype = "ENSEMBL")

# Convert the ENSEMBL IDs to gene symbols
background_gene_symbols <- mapIds(org.Mm.eg.db, keys = background_genes, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# Remove any NA values (genes that don't have a corresponding symbol)
background_gene_symbols <- background_gene_symbols[!is.na(background_gene_symbols)]

# Now, 'background_gene_symbols' is your list of all mouse gene symbols for use as background
head(background_gene_symbols)

# Save the background gene symbols to a CSV file
write.csv(background_gene_symbols, "background_gene_symbols.csv", row.names = FALSE)

# ######################################
#  Prepare data for clusterProfiler
# ######################################

# Load DEG data (Filtered) and Background Genes
deg_data <- read.csv("DEG_results_filtered.csv")

# Get the Ensembl IDs for your DEG (those that are up- or down-regulated)
# Assuming your DEG data has Ensembl IDs in a column named 'Ensembl_ID'
genes_in_data <- deg_data$Geneid

# Check if your Ensembl IDs are in the background
# Assuming background_genes is a vector of Ensembl IDs
genes_in_background <- genes_in_data[genes_in_data %in% background_genes]

# Optional: Check the number of matching genes
length(genes_in_background)  # This gives you how many of your DEGs are in the background gene list

# Convert DEG Ensembl IDs to Entrez IDs
deg_entrez <- bitr(genes_in_background, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Convert background Ensembl IDs to Entrez IDs
background_entrez <- bitr(background_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# Extract Entrez IDs for DEGs (the relevant genes for enrichment)
deg_entrez_ids <- deg_entrez$ENTREZID

# Extract Entrez IDs for background genes
background_entrez_ids <- background_entrez$ENTREZID


# ######################################
#  clusterProfiler
# ############ ############ ############

#Gene Ontology Enrichment
go_enrichment <- enrichGO(
  gene = deg_entrez_ids,
  OrgDb = org.Mm.eg.db,
  keyType = "ENTREZID",
  ont = "BP",  # Biological Process
  pAdjustMethod = "BH",  # Benjamini-Hochberg correction
  universe = background_entrez_ids,  # Background gene set
  qvalueCutoff = 0.05,  # Adjusted p-value threshold
  minGSSize = 5
)

# ######################################
#  Add diffexpressed column and save
# ############ ############ ############

#Add diffexpressed column in the GO results
go_enrichment_df <- go_enrichment@result  # Convert the result to a data frame

# Convert ENSEMBL IDs from deg_data to ENTREZ IDs
deg_data$ENTREZID <- mapIds(org.Mm.eg.db, keys = deg_data$Geneid, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")

# Add 'diffexpressed' column based on matched ENTREZ IDs
go_enrichment_df$diffexpressed <- sapply(go_enrichment_df$geneID, function(gene_ids) {
  # Split gene IDs (they may be separated by "/")
  genes <- unlist(strsplit(gene_ids, split = "/"))
  
  # Get the diffexpressed values from deg_data for the matched ENTREZ IDs
  diff_status <- deg_data$diffexpressed[deg_data$ENTREZID %in% genes]
  
  # If there are matched genes, return the most frequent diffexpressed status
  if (length(diff_status) > 0) {
    return(names(sort(table(diff_status), decreasing = TRUE))[1])  # Most common diffexpressed status
  } else {
    return(NA)  # If no matching genes, assign NA
  }
})

# Save the GO enrichment results as CSV
write.csv(go_enrichment_df, file = file.path(out_path, "GO_enrichment_results_with_diffexpressed.csv"), row.names = FALSE)

# Save the GO enrichment results as RDS
saveRDS(go_enrichment, file = file.path(out_path, "GO_enrichment_results_with_diffexpressed.rds"))


# ######################################
#  Visualization
# ############ ############ ############


# ######################################
#  Barplot 
# ############ ############ ############

#Go enrichment
# Calculate -log10(qvalue) for better visualization
go_enrichment_df <- go_enrichment_df %>%
  mutate(log_qvalue = -log10(qvalue))

# Filter the top 50 GO terms based on log_qvalue
top_50_go <- go_enrichment_df %>%
  arrange(desc(log_qvalue)) %>%
  head(50)

# Create the barplot with color gradient based on p.adjust and show the legend
ggplot(top_50_go, aes(x = reorder(Description, log_qvalue), y = log_qvalue, fill = p.adjust)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # Flip the barplot to make it horizontal
  labs(x = "GO Term", y = "-log10(q-value)", title = "Top 50 GO Enrichment Terms") +
  scale_fill_gradient(low = "blue", high = "red", name = "p.adjust") +  # Adding color gradient legend
  theme_minimal() +
  theme(legend.position = "right")  # Place the legend to the right


# ######################################
#  Separately
# ############ ############ ############

# Filter GO terms for upregulated and downregulated genes
go_up <- go_enrichment_df %>% filter(diffexpressed == "UP")
go_down <- go_enrichment_df %>% filter(diffexpressed == "DOWN")

# Select top 10 GO terms for each group
go_up_top <- go_up %>% arrange(p.adjust) %>% head(10)
go_down_top <- go_down %>% arrange(p.adjust) %>% head(10)

# Combine upregulated and downregulated data
go_combined <- bind_rows(
  go_up_top %>% mutate(Regulation = "Upregulated"),
  go_down_top %>% mutate(Regulation = "Downregulated")
)

# Create barplot
ggplot(go_combined, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    title = "Top GO Enrichment Terms by Regulation"
  ) +
  scale_fill_manual(values = c("Upregulated" = "red", "Downregulated" = "blue")) +
  theme_minimal()

# ######################################
#  Dotplot
# ############ ############ ############

# Dotplot for GO
dotplot(go_enrichment, showCategory = 20)  # Direct use of go_enrichment object

# ######################################
#  Cnetplot
# ############ ############ ############
# Cnetplot for GO
cnetplot(go_enrichment, 
         showCategory = 10,  # Show the top 10 GO terms
         foldChange = go_enrichment@result$FoldEnrichment,  # Use fold enrichment for node coloring (or logFC if preferred)
         circular = TRUE,    # Circular layout for clarity
         node_label = "category") +  # Display GO process names as labels
  ggtitle("Cnetplot: Top 10 GO Processes") +
  theme_minimal()

# ######################################
#  Treeplots
# ############ ############ ############

#Treeplot Go
# Calculate pairwise similarities of the enriched GO terms
go_enrichment_sim <- pairwise_termsim(go_enrichment)

# Create the treeplot for the GO enrichment results
treeplot(go_enrichment_sim)


# ######################################
#  Enrichment maps
# ############ ############ ############

# Create the Enrichment Map for GO Enrichment
emapplot(go_enrichment_sim)

# ######################################
#  Upset plot
# ############ ############ ############

#Upset plot for GO enrichment
upsetplot(go_enrichment)

