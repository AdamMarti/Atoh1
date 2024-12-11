setwd('/home/annie/Documents/PGB/Project/BLAST results')
blast_results <- read.csv('species.csv')

# Load necessary libraries
library(ggplot2)
library(pheatmap)
library(dplyr)
library(tidyr)

# Create the 'Group' column based on the species
blast_results$Group <- ifelse(blast_results$Scientific.Name %in% c("Pongo pygmaeus", "Gorilla gorilla gorilla", "Pan paniscus"), "Primates", 
                              ifelse(blast_results$Scientific.Name %in% c("Polyodon spathula", "Albula goreensis", "Chelydra serpentina"), "Fish & Reptiles", 
                                     "Non-primate Mammals"))

# Ensure the 'Group' column is a factor to control the order of groups
blast_results$Group <- factor(blast_results$Group, levels = c("Primates", "Non-primate Mammals", "Fish & Reptiles"))

# Check the structure of the data to confirm everything looks correct
str(blast_results)

# Plot the barplot with species and percentage identity, grouped by the 'Group' column
ggplot(blast_results, aes(x = Scientific.Name, y = Per..ident, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +  # Group bars side by side within each species category
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate species names for readability
  labs(title = "Percentage Identity by Species", x = "Species", y = "Percentage Identity") +  # Add title and axis labels
  scale_fill_manual(values = c("Primates" = "#ADD8E6",  # Light blue for primates
                               "Non-primate Mammals" = "#98FB98",  # Light green for non-primate mammals
                               "Fish & Reptiles" = "#FFA07A")) +  # Light orange for fish & reptiles
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +  # Rotate x-axis labels for better readability
  facet_wrap(~Group, scales = "free_x")  # Create separate panels for each group, grouped together
