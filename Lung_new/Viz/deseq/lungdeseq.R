setwd("C:/Users/sanyu/OneDrive/Desktop/Lung_new/Viz/deseq")

library(ggplot2)
library(ggrepel)
library(dplyr)

# Load the data (full DESeq2 result with all taxa)
# Expecting columns:
# baseMean, log2FoldChange, lfcSE, pvalue, padj, genus, significant, direction
data <- read.csv("imp10.csv", header = TRUE)

# We will use padj as the adjusted p-value
padj_cutoff <- 0.05

# Add -log10(padj) and recompute direction
data <- data %>%
  mutate(
    padj = padj,  # just to mirror MaAsLin syntax
    
    # Y-axis for volcano: -log10(adjusted p-value)
    neglog10padj = -log10(padj),
    
    # Direction based on log2FoldChange sign + padj (significant ones)
    direction_label = case_when(
      padj < padj_cutoff & log2FoldChange > 0 ~ "Tumor_enriched",
      padj < padj_cutoff & log2FoldChange < 0 ~ "NAT_enriched",
      TRUE                                   ~ "Not Significant"
    )
  )

# Determine which genera to label
# If there are significant taxa, label those
# Otherwise, label the top 10 by lowest padj
n_significant <- sum(data$padj < padj_cutoff, na.rm = TRUE)

if (n_significant > 0) {
  label_data <- data %>%
    filter(padj < padj_cutoff)
  cat("Found", n_significant, "significant genera to label\n")
} else {
  label_data <- data %>%
    arrange(padj) %>%
    slice_head(n = 10)
  cat("No significant genera found. Labeling top 10 by padj\n")
}

# Keep only finite values for plotting
plot_data <- data %>%
  filter(is.finite(log2FoldChange), is.finite(neglog10padj))

# Make direction_label an ordered factor for consistent legend
plot_data$direction_label <- factor(
  plot_data$direction_label,
  levels = c("NAT_enriched", "Tumor_enriched", "Not Significant")
)

# Color palette
dir_cols <- c(
  "Tumor_enriched"  = "#E64B35",  # Red
  "NAT_enriched"    = "#4DBBD5",  # Blue
  "Not Significant" = "grey70"    # Grey
)

# Create the volcano plot
volcano_deseq <- ggplot(plot_data,
                        aes(x = log2FoldChange,
                            y = neglog10padj)) +
  
  # Points colored by direction
  geom_point(aes(color = direction_label),
             size = 2,
             alpha = 0.6) +
  
  # Manual colors
  scale_color_manual(
    values = dir_cols,
    name   = "Direction",
    drop   = FALSE
  ) +
  
  # Horizontal line at significance threshold (padj = 0.05)
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  # Vertical line at log2FC = 0
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  # Labels for significant genera (or top 10)
  geom_text_repel(
    data = label_data,
    aes(label = genus, color = direction_label),
    size          = 3,
    max.overlaps  = 20,
    box.padding   = 0.5,
    point.padding = 0.3,
    segment.color = "grey50",
    segment.size  = 0.2,
    show.legend   = FALSE
  ) +
  
  # Titles and axis labels
  labs(
    title = "DESeq2 Volcano Plot, Lung Imputed: Tumor vs NAT (10% filter)",
    x     = "log2 fold change (Tumor − NAT)",
    y     = expression(-log[10]~"(adjusted p-value, padj)")
  ) +
  
  # Clean theme
  theme_minimal() +
  theme(
    plot.title       = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title       = element_text(size = 12, face = "bold"),
    axis.text        = element_text(size = 10),
    legend.title     = element_text(size = 11, face = "bold"),
    legend.text      = element_text(size = 10),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

# Show plot
print(volcano_deseq)

# Save plot
ggsave("Vpimp10pct.png", plot = volcano_deseq,
       width = 10, height = 8, dpi = 300)
ggsave("Vpimp10pct.pdf", plot = volcano_deseq,
       width = 10, height = 8)

