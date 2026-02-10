setwd("C:/Users/sanyu/OneDrive/Desktop/Ovary new/Viz/maaslin")
library(ggplot2)
library(ggrepel)
library(dplyr)

# Load the data (full MaAsLin2 result table)
# Columns: feature, metadata, value, coef, stderr, N, N.not.0,
#          pval, qval, significant, direction
data <- read.csv("imp10pct.csv", header = TRUE)

# Use qval as adjusted p-value
padj_cutoff <- 0.05

data <- data %>%
  mutate(
    # adjusted p-value
    padj = qval,
    
    # y-axis: -log10(padj)
    neglog10padj = -log10(padj),
    
    # direction based on coef sign + padj
    direction_label = case_when(
      padj < padj_cutoff & coef > 0 ~ "Tumor_enriched",
      padj < padj_cutoff & coef < 0 ~ "NAT_enriched",
      TRUE                          ~ "Not Significant"
    )
  )

# Decide which genera to label:
# if any significant, label those; otherwise label top 10 by lowest padj
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

# Keep only finite values
plot_data <- data %>%
  filter(is.finite(coef), is.finite(neglog10padj))

# Ordered factor for legend
plot_data$direction_label <- factor(
  plot_data$direction_label,
  levels = c("NAT_enriched", "Tumor_enriched", "Not Significant")
)

# Colors
dir_cols <- c(
  "Tumor_enriched"  = "#E64B35",  # red
  "NAT_enriched"    = "#4DBBD5",  # blue
  "Not Significant" = "grey70"    # grey
)

# Volcano plot
volcano_maaslin <- ggplot(plot_data,
                          aes(x = coef,
                              y = neglog10padj)) +
  
  geom_point(aes(color = direction_label),
             size = 2,
             alpha = 0.6) +
  
  scale_color_manual(
    values = dir_cols,
    name   = "Direction",
    drop   = FALSE
  ) +
  
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
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
  
  labs(
    title = "MaAsLin2 Volcano Plot, Ovary Imputed : Tumor vs NAT (10% filter)",
    x     = "MaAsLin2 coefficient (Tumor − NAT)",
    y     = expression(-log[10]~"(adjusted p-value, qval)")
  ) +
  
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
print(volcano_maaslin)

# Save plot
ggsave("vpimp10pct_maaslin.png", plot = volcano_maaslin,
       width = 10, height = 8, dpi = 300)
ggsave("vpimp10pct_maaslin.pdf", plot = volcano_maaslin,
       width = 10, height = 8)

