setwd("C:/Users/sanyu/OneDrive/Desktop/Breast_new/Viz/ANCOM")

library(ggplot2)
library(ggrepel)
library(dplyr)

# load the data (this is the full ALDEx table with all taxa)
data <- read.csv("raw5pct.csv", header = TRUE)

# In ALDEx we use qval as the q val
padj_cutoff <- 0.05

# add padj and -ve log col and direction 
data <- data %>%
  mutate(
    # choose which adjusted p-value to use (here: qval from ALDEx)
    padj = qval,
    
    # y-axis for volcano: -log10(adjusted p-value)
    neglog10padj = -log10(padj),
    
    # direction based on effect sign + padj
    direction_label = case_when(
      padj >= padj_cutoff        ~ "Not Significant",
      lfc > 0                    ~ "Tumor_enriched",
      lfc < 0                    ~ "NAT_enriched",
      TRUE                       ~ "Not Significant"
    ),
    
    # make it an ordered factor so all three levels exist,
    # even if one of them has 0 points (e.g. NAT_enriched)
    direction_label = factor(
      direction_label,
      levels = c("NAT_enriched", "Tumor_enriched", "Not Significant")
    )
  )

## ========= CHOOSE TAXA TO LABEL =====================================

# 1) all significant taxa (padj < cutoff)
sig_data <- data %>%
  filter(padj < padj_cutoff)

# 2) extra non-significant taxa: top 10 lowest padj (highest -log10)
extra_data <- data %>%
  filter(padj >= padj_cutoff) %>%      # only non-significant
  arrange(desc(neglog10padj)) %>%      # strongest p-values first
  slice_head(n = 10)                   # change 10 if you want more / fewer

# 3) combine sig + extra and remove duplicates / NA taxon
label_data <- bind_rows(sig_data, extra_data) %>%
  filter(!is.na(taxon)) %>%
  distinct(taxon, .keep_all = TRUE)

## ====================================================================

# making of plot
# 1) initialize a ggplot object
volcano_aldex <- ggplot(data,
                        aes(x = lfc,
                            y = neglog10padj)) +
  
  # data points
  geom_point(aes(color = direction_label),
             size = 2,
             alpha = 0.6) +
  
  # manually set the dot colors for each direction
  scale_color_manual(
    values = c(
      "Tumor_enriched"  = "#E64B35",
      "NAT_enriched"    = "#4DBBD5",
      "Not Significant" = "grey70"
    ),
    name   = "Direction",
    # drop = FALSE keeps levels in legend even if they have 0 points
    drop   = FALSE
  ) +
  
  # horizontal line at significance threshold
  geom_hline(
    yintercept = -log10(padj_cutoff),
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  # vertical line at effect = 0
  geom_vline(
    xintercept = 0,
    linetype   = "dashed",
    color      = "black",
    linewidth  = 0.5
  ) +
  
  # titles and axis labels
  labs(
    title = "ANCOMBC Volcano Plot Raw: Tumor vs NAT (5% filter)",
    x     = "ANCOMBC lfc (Tumor − NAT, CLR scale)",
    y     = expression(-log[10]~"(adjusted p-value, qval)")
  ) +
  
  # clean plot
  theme_minimal() +
  theme(
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.title      = element_text(size = 12, face = "bold"),
    axis.text       = element_text(size = 10),
    legend.title    = element_text(size = 11, face = "bold"),
    legend.text     = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank()
  )

# add labels for significant + extra selected taxa
if (nrow(label_data) > 0) {
  volcano_aldex <- volcano_aldex +
    geom_text_repel(
      data = label_data,
      aes(label = taxon, color = direction_label),
      size          = 3,
      max.overlaps  = 20,
      box.padding   = 0.5,   # space around the text box
      point.padding = 0.3,   # space around the point
      segment.color = "grey50",
      segment.size  = 0.2,
      show.legend   = FALSE  # don't add text labels to the legend
    )
}

# Display the plot
print(volcano_aldex)

# Save plot to files
ggsave("vp_raw5pct.png", plot = volcano_aldex,
       width = 10, height = 8, dpi = 300)
ggsave("vp_raw5pct.pdf", plot = volcano_aldex,
       width = 10, height = 8)

