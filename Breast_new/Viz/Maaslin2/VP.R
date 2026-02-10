setwd("C:/Users/sanyu/OneDrive/Desktop/Breast_new/Viz/Maaslin2")

library(ggplot2)
library(ggrepel)
library(dplyr)

# Load the data (this is the full MaAsLin2 table with all taxa)
# Columns: feature, metadata, value, coef, stderr, N, N.not.0, pval, qval, significant, direction
data <- read.csv("imp5pct.csv", header = TRUE)

# We will use qval as the adjusted p-value
padj_cutoff <- 0.05

# Add padj and -log10 column and direction
data <- data %>%
  mutate(
    # Choose which adjusted p-value to use (here: qval from MaAsLin2)
    padj = qval,
    
    # Y-axis for volcano: -log10(adjusted p-value)
    neglog10padj = -log10(padj),
    
    # Direction based on coef sign + padj (significant ones)
    direction_label = case_when(
      padj < padj_cutoff & coef > 0 ~ "Tumor_enriched",
      padj < padj_cutoff & coef < 0 ~ "NAT_enriched",
      TRUE                          ~ "Not Significant"
    )
  )

# Determine which features to label
# If there are significant taxa, label those
# Otherwise, label the top 10 by lowest qval
n_significant <- sum(data$padj < padj_cutoff, na.rm = TRUE)

if (n_significant > 0) {
  # Label significant features
  label_data <- data %>% 
    filter(padj < padj_cutoff)
  cat("Found", n_significant, "significant genera to label\n")
} else {
  # No significant features, so label top 10 by qval
  label_data <- data %>%
    arrange(padj) %>%
    slice_head(n = 10)
  cat("No significant genera found. Labeling top 10 by qval\n")
}

# Keep only finite values for plotting
plot_data <- data %>%
  filter(is.finite(coef), is.finite(neglog10padj))


# Make direction_label an ordered factor for consistent legend
plot_data$direction_label <- factor(
  plot_data$direction_label,
  levels = c("NAT_enriched", "Tumor_enriched", "Not Significant")
)

# Color palette
dir_cols <- c(
  "Tumor_enriched"  = "#E64B35",   # Red
  "NAT_enriched"    = "#4DBBD5",   # Blue
  "Not Significant" = "grey70"      # Grey
)

# Create the volcano plot
volcano_maaslin <- ggplot(plot_data,
                          aes(x = coef,
                              y = neglog10padj)) +
  
  # Add data points colored by direction
  geom_point(aes(color = direction_label),
             size = 2,
             alpha = 0.6) +
  
  # Manually set the dot colors for each direction
  scale_color_manual(
    values = dir_cols,
    name   = "Direction",
    drop   = FALSE
  ) +
  
  # Horizontal line at significance threshold (qval = 0.05)
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  # Vertical line at coef = 0 (no difference)
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  # Add labels for significant genera (or top 10 if none significant)
  geom_text_repel(
    data = label_data,
    aes(label = feature, color = direction_label),
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
    title = "MaAsLin2 Volcano Plot, Imputed: Tumor vs NAT (5% filter)",
    x     = "MaAsLin2 coefficient (Tumor − NAT)",
    y     = expression(-log[10]~"(adjusted p-value, qval)")
  ) +
  
  # Clean plot theme
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

# Display the plot
print(volcano_maaslin)

# Save plot to files
ggsave("vpimp5pct.png", plot = volcano_maaslin,
       width = 10, height = 8, dpi = 300)
ggsave("vpimp5pct.pdf", plot = volcano_maaslin,
       width = 10, height = 8)

