# Load necessary libraries
library(ggplot2)
library(dplyr)
library(gridExtra)  # For arranging plots side by side

# Set working directory (adjust path as necessary)
setwd("~/PGB/Project/dNdS/")

# Load your files
main_data <- read.table("mart_export.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
second_file <- read.table("HumanTFs_DBD.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

head(main_data)
head(second_file)

# Calculate dN/dS for all genes
main_data$dN_dS <- main_data$dN / main_data$dS

# Gene of interest
gene_of_interest <- "ENSG00000172238"
dnds_interest <- main_data$dN_dS[main_data$Ensembl.Gene.ID == gene_of_interest]

# Plot function with specified customizations
plot_violin <- function(dnds_data, label, dnds_interest) {
  dnds_df <- data.frame(dnds = dnds_data, group = label)
  
  ggplot(dnds_df, aes(x = group, y = dnds)) +
    geom_violin(fill = "lightblue", color = "blue", alpha = 0.5) +  # Blue violin plot
    geom_boxplot(width = 0.1, color = "black", alpha = 0.3, outlier.shape = NA, fatten = 2, na.rm = FALSE) +  # Box plot inside violin
    geom_hline(yintercept = dnds_interest, color = "darkred", linetype = "dashed", size = 0.5) +
    geom_text(aes(x = 1, y = dnds_interest, label = paste("Atoh1", round(dnds_interest, 3))), 
              vjust = -0.5, hjust = -0.25, size = 4, color = "darkred") +# Red dashed line
    labs(title = label, x = "" , y = "dN/dS") +  # Title and y-axis label
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      plot.title = element_text(size = 10, face = "bold"),  # Beautiful text theme
      legend.position = "none",
    ) +
    ylim(0, 1)  # Limit y-axis to 1
}

# Filter subsets for each condition
matching_genes <- main_data[main_data$Ensembl.Gene.ID %in% second_file$Ensembl_ID, ]
matching_bHLH <- matching_genes[grepl("bHLH", second_file$DBD[match(matching_genes$Ensembl.Gene.ID, second_file$Ensembl_ID)]), ]

# Generate individual plots
plot_all <- plot_violin(main_data$dN_dS, "dN/dS of H.Sapiens genes", dnds_interest)
plot_TFs <- plot_violin(matching_genes$dN_dS, "dN/dS of H.Sapiens TFs", dnds_interest)
plot_bHLH <- plot_violin(matching_bHLH$dN_dS, "dN/dS of H.Sapiens bHLH TFs", dnds_interest)

print(plot_all)
print(plot_TFs)
print(plot_bHLH)

# Arrange the three plots in a single row
grid.arrange(plot_all, plot_TFs, plot_bHLH, ncol = 3)


# Combina los datos de main_data y second_file usando una clave común (por ejemplo, Ensembl.Gene.ID)
merged_data <- merge(main_data, second_file, by.x = "Ensembl.Gene.ID", by.y = "Ensembl_ID")


# Step 1: Count the number of components in each TF family
family_counts <- merged_data %>%
  group_by(DBD) %>%
  summarise(count = n())

# Step 2: Filter families to include only those with more than 100 components
large_families <- family_counts %>%
  filter(count > 100) %>% #HERE YOU FILTER THE FAMILIES THAT HAVE LESS THAN N MEMBERS 
  pull(DBD)

# Step 3: Filter the merged data to include only large families
filtered_data <- merged_data %>%
  filter(DBD %in% large_families)

# Step 4: Generate dN/dS summary statistics for large families
tf_summary <- filtered_data %>%
  group_by(DBD) %>%
  summarise(
    mean_dN_dS = mean(dN_dS),
    sd_dN_dS = sd(dN_dS),
    n = n()
  ) %>%
  arrange(desc(mean_dN_dS))

# Display the summary statistics
print(tf_summary)

# Step 5: Generate violin plots for large families
ggplot(filtered_data, aes(x = DBD, y = dN_dS)) +
  geom_violin(fill = "lightblue", color = "blue", alpha = 0.5) +
  geom_boxplot(width = 0.1, color = "black", alpha = 0.3, outlier.shape = NA, fatten = 2) +
  labs(title = "dN/dS Distributions for TF Families with +100 members", x = "TF Family", y = "dN/dS") +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    strip.text = element_text(size = 10, face = "bold")
  ) +
  facet_wrap(~ DBD, scales = "free_x") +
  ylim(0, 1)

ggsave("Large_TF_Families_dN_dS_Violin.png", width = 14, height = 12, dpi = 300)

# Save the plot with a white background using ggsave()
ggsave("Large_TF_Families_dN_dS_Violin_White.png", width = 14, height = 12, dpi = 300, bg = "white")

# Step 6: Save the summary statistics to a CSV file
write.csv(tf_summary, "Large_TF_families_dN_dS_summary.csv", row.names = FALSE)


# Overall test for differences between groups
kruskal_result <- kruskal.test(dN_dS ~ DBD, data = filtered_data)
print(kruskal_result)

# If there are significant differences, perform pairwise comparisons with Wilcoxon
if (kruskal_result$p.value < 0.05) {
  # Pairwise Wilcoxon test with Bonferroni adjustment
  pairwise_results <- pairwise.wilcox.test(
    x = filtered_data$dN_dS, 
    g = filtered_data$DBD, 
    p.adjust.method = "bonferroni"
  )
  
  print(pairwise_results)
  
  # Export the entire comparisons matrix to a CSV file
  pairwise_matrix <- as.data.frame(pairwise_results$p.value)
  
  # Add the names of the rows and columns to the matrix
  rownames(pairwise_matrix) <- rownames(pairwise_results$p.value)
  colnames(pairwise_matrix) <- colnames(pairwise_results$p.value)
  
  # Export the entire matrix to a CSV file
  write.csv(pairwise_matrix, 
            file = "pairwise_wilcoxon_full_results.csv", 
            row.names = TRUE)
  
}
