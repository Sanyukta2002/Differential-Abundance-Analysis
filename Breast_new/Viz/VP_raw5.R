setwd("C:/Users/sanyu/OneDrive/Desktop/Breast_new/Viz/DESEQ")
# Volcano Plot for DESeq2 Results
# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)

# Read your data
# Replace 'your_file.csv' with your actual file path
data <- read.csv("raw5_all.csv", header = TRUE)

# Add a column for significance and direction
data <- data %>%
  mutate(
    # Determine significance based on padj < 0.05
    significant = ifelse(padj < 0.05, "Significant", "Not Significant"),
    
    # Create a combined column for direction labeling
    direction_label = case_when(
      padj >= 0.05 ~ "Not Significant",
      log2FoldChange > 0 ~ "Tumor_enriched",
      log2FoldChange < 0 ~ "NAT_enriched",
      TRUE ~ "Not Significant"
    )
  )

# Create the volcano plot
volcano_plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
  
  # Add points colored by direction
  geom_point(aes(color = direction_label), size = 2, alpha = 0.6) +
  
  # Set custom colors
  scale_color_manual(
    values = c(
      "Tumor_enriched" = "#E64B35",  # Red
      "NAT_enriched" = "#4DBBD5",     # Blue
      "Not Significant" = "grey70"
    ),
    name = "Direction"
  ) +
  
  # Add horizontal line for significance threshold (padj = 0.05)
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Add vertical line at log2FC = 0
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  
  # Add genus labels for significant points only
  geom_text_repel(
    data = subset(data, padj < 0.05),
    aes(label = genus, color = direction_label),
    size = 3,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size = 0.2,
    show.legend = FALSE
  ) +
  
  # Set axis labels
  labs(
    title = "Volcano Plot: Tumor vs Normal Adjacent Tissue  , Raw data w 5 pct Filter",
    x = expression(log[2]~"Fold Change"),
    y = expression(-log[10]~"(adjusted p-value)")
  ) +
  
  # Use a clean theme
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    legend.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# Display the plot
print(volcano_plot)

# Optional: Save the plot
ggsave("volcano_plot.png", plot = volcano_plot, width = 10, height = 8, dpi = 300)
ggsave("volcano_plot.pdf", plot = volcano_plot, width = 10, height = 8)
