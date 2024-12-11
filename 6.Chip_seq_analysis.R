setwd("~/Desktop/PGB/PGB_FINAL_PROJECT/CHIP")

#This code is pretty similar to the one available in ChIP-Seeker, only using the blocks of code needed to get 
#the grephs of interest

#Loading necessary libraries

library(ChIPseeker)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
library(clusterProfiler)
library(ggupset)
library(rtracklayer)
library(ggimage)

#Choosing our files
files <- import("~/Desktop/PGB/PGB_FINAL_PROJECT/CHIP/IP_summits.bed")

#Representing the Chip Peaks over different chromosomes 
peak <- readPeakFile("~/Desktop/PGB/PGB_FINAL_PROJECT/CHIP/IP_summits.bed")
peak
covplot(peak, weightCol="V5")

#Peak Annotation

peakAnno <- annotatePeak(files, 
                         tssRegion = c(-1000, 1000), 
                         TxDb = txdb,
                         annoDb = "org.Mm.eg.db")

plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
upsetplot(peakAnno)
upsetplot(peakAnno, vennpie = TRUE)

plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")


#Functional enrichment analysis

# Load required packages
library(org.Mm.eg.db)  # Mouse genome annotation
library(clusterProfiler)
library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId, organism = "mouse")
head(pathway1, 2)

gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene, organism = "mouse")
head(pathway2, 2)

dotplot(pathway2)

#Extract Gene IDs from Functional Enrichment Analysis

# Extracting gene IDs involved in enriched pathways
library(dplyr)
library(tidyr)  # Load tidyr for the unnest() function

# Extract genes related to each pathway from the enrichment result
enriched_genes <- pathway2@result
head(enriched_genes)

# Separate the gene IDs for each enriched pathway
enriched_genes_df <- enriched_genes %>%
  select(ID, Description, geneID) %>%
  mutate(geneID = strsplit(geneID, "/")) %>%  # Split gene IDs by "/"
  unnest(geneID)  # Flatten the list into a dataframe

# Display the resulting dataframe
head(enriched_genes_df)

# Save the gene IDs to a CSV file
write.csv(enriched_genes_df, "enriched_genes_mouse.csv", row.names = FALSE)


#Get Ensembl IDs

#This final part was for getting the Ensembl IDs of the genes found so that they could be used in other parts of the project
#such as RNA-Seq comparison with ChIP-Seq

# Load the libraries
library(org.Mm.eg.db)
library(dplyr)

# Extract Entrez IDs from the dataframe
entrez_ids <- enriched_genes_df$geneID

# Map Entrez IDs to Ensembl IDs
ensembl_ids <- mapIds(
  org.Mm.eg.db,
  keys = entrez_ids,
  column = "ENSEMBL",
  keytype = "ENTREZID",
  multiVals = "first"  # If multiple mappings, take the first one
)

# Add the Ensembl IDs as a new column to your dataframe
enriched_genes_df$Ensembl_ID <- ensembl_ids

# Display the updated dataframe
head(enriched_genes_df)

write.csv(enriched_genes_df, "enriched_genes_with_ensembl.csv", row.names = FALSE)

