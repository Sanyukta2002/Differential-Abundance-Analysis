setwd("C:/Users/sanyu/OneDrive/Desktop/Ovary new/Viz/aldex")
library(ggplot2)
library(ggrepel)
library(dplyr)


# Load the data (full ALDEx2 table with all taxa)
# Columns: taxon, we.ep, we.eBH, wi.ep, wi.eBH, rab.*, diff.*, effect,
#          overlap, p_we, q_we, p_wi, q_wi, p_min, significant, direction
data <- read.csv("imp10pct.csv", header = TRUE)

# Use we.eBH as the adjusted p-value
padj_cutoff <- 0.05

# Add padj, -log10(padj) and direction label
data <- data %>%
  mutate(
    # adjusted p-value (Welch test BH)
    padj = we.eBH,
    
    # Y-axis: -log10(adjusted p-value)
    neglog10padj = -log10(padj),
    
    # Direction based on effect sign + padj
    direction_label = case_when(
      padj < padj_cutoff & effect > 0 ~ "Tumor_enriched",
      padj < padj_cutoff & effect < 0 ~ "NAT_enriched",
      TRUE                            ~ "Not Significant"
    )
  )

# Decide which taxa to label:
# If there are significant taxa, label those,
# otherwise label top 10 by lowest padj
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
  filter(is.finite(effect), is.finite(neglog10padj))

# Ordered factor for legend
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

# Volcano plot
volcano_aldex <- ggplot(plot_data,
                        aes(x = effect,
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
  
  # Horizontal line at padj = 0.05
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  # Vertical line at effect = 0
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  # Labels for significant genera (or top 10)
  geom_text_repel(
    data = label_data,
    aes(label = taxon, color = direction_label),
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
    title = "ALDEx2 Volcano Plot,Ovary Impute: Tumor vs NAT (10% filter)",
    x     = "ALDEx2 effect size (Tumor − NAT)",
    y     = expression(-log[10]~"(adjusted p-value, we.eBH)")
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
print(volcano_aldex)

# Save plot
ggsave("vpimp10pct_aldex.png", plot = volcano_aldex,
       width = 10, height = 8, dpi = 300)
ggsave("vpimp10pct_aldex.pdf", plot = volcano_aldex,
       width = 10, height = 8)




























































































