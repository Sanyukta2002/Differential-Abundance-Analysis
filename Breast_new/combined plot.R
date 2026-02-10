## ------------------------------------------------------------
## Combined script: Dot plot + UpSet plot
##  - Uses same four DA result files
##  - Saves dot plot (PNG + PDF)
##  - Saves UpSet plot (PDF)
## ------------------------------------------------------------

## 0. Setup ----------------------------------------------------

setwd("C:/Users/sanyu/OneDrive/Desktop/Impute check")

suppressPackageStartupMessages({
  library(tidyverse)  # dplyr, tidyr, ggplot2, etc.
  library(UpSetR)
})

alpha <- 0.05  # significance threshold


## 1. Read input files -----------------------------------------

aldex   <- read.csv("ALDEx2_prev10_imputed_Tumor_vs_NAT_all.csv",
                    check.names = FALSE, stringsAsFactors = FALSE)
ancom   <- read.csv("ANCOMBC2_prev10_imputed_Tumor_vs_NAT_full.csv",
                    check.names = FALSE, stringsAsFactors = FALSE)
deseq   <- read.csv("DESeq2_prev10_imputed_Tumor_vs_NAT_all.csv",
                    check.names = FALSE, stringsAsFactors = FALSE)
maaslin <- read.csv("Maaslin2_prev10_imputed_Tumor_vs_NAT_all.csv",
                    check.names = FALSE, stringsAsFactors = FALSE)


## 2. Extract significant taxa per tool ------------------------
## (These are used both for the overlap matrix and for the UpSet)

# ALDEx2: use the logical 'significant' column + 'direction'
aldex_sig <- aldex %>%
  filter(significant) %>%
  transmute(
    taxon,
    Direction_ALDEx2 = direction,
    ALDEx2 = 1L
  )

# ANCOMBC2: qval < alpha
ancom_sig <- ancom %>%
  filter(qval < alpha) %>%
  transmute(
    taxon,
    Direction_ANCOMBC2 = direction,
    ANCOMBC2 = 1L
  )

# DESeq2: padj < alpha
deseq_sig <- deseq %>%
  filter(padj < alpha) %>%
  transmute(
    taxon = genus,
    Direction_DESeq2 = direction,
    DESeq2 = 1L
  )

# MaAsLin2: qval < alpha, and if metadata/value present keep only Group = Tumor
maaslin_sig <- maaslin %>%
  {
    if (all(c("metadata", "value") %in% names(.))) {
      filter(., metadata == "Group", value == "Tumor")
    } else {
      .
    }
  } %>%
  filter(qval < alpha) %>%
  transmute(
    taxon = feature,
    Direction_Maaslin2 = direction,
    Maaslin2 = 1L
  )


## 3. Build overlap matrix -------------------------------------

overlap_mat <- aldex_sig %>%
  full_join(ancom_sig,   by = "taxon") %>%
  full_join(deseq_sig,   by = "taxon") %>%
  full_join(maaslin_sig, by = "taxon") %>%
  mutate(
    ALDEx2   = ifelse(is.na(ALDEx2),   0L, ALDEx2),
    ANCOMBC2 = ifelse(is.na(ANCOMBC2), 0L, ANCOMBC2),
    DESeq2   = ifelse(is.na(DESeq2),   0L, DESeq2),
    Maaslin2 = ifelse(is.na(Maaslin2), 0L, Maaslin2)
  ) %>%
  rowwise() %>%
  mutate(
    Direction = {
      dirs <- c(
        Direction_ALDEx2,
        Direction_ANCOMBC2,
        Direction_DESeq2,
        Direction_Maaslin2
      )
      dirs <- unique(na.omit(dirs))
      if (length(dirs) == 1) dirs else "NA"
    },
    n_tools = sum(c_across(c(ALDEx2, ANCOMBC2, DESeq2, Maaslin2)))
  ) %>%
  ungroup() %>%
  arrange(desc(n_tools), taxon)

# Filter only taxa significant in > 1 tool for the dot plot
overlap_gt1 <- overlap_mat %>%
  filter(n_tools > 1) %>%
  select(taxon, ALDEx2, ANCOMBC2, DESeq2, Maaslin2, n_tools, Direction)

cat("\nTaxa significant in > 1 tool:\n")
print(overlap_gt1)


## 4. DOT PLOT (your compact overlap plot) ---------------------

dot_df <- overlap_gt1 %>%
  pivot_longer(
    cols = c(ALDEx2, ANCOMBC2, DESeq2, Maaslin2),
    names_to = "Tool",
    values_to = "Sig"
  )

# Order taxa (top = most overlapped)
taxa_order <- overlap_gt1 %>%
  arrange(desc(n_tools), taxon) %>%
  pull(taxon)

dot_df <- dot_df %>%
  mutate(
    taxon = factor(taxon, levels = rev(taxa_order)),
    Tool  = factor(Tool, levels = c("ALDEx2", "ANCOMBC2", "DESeq2", "Maaslin2"))
  )

p_dot <- ggplot(dot_df, aes(x = Tool, y = taxon)) +
  # light grey background dots (non-significant)
  geom_point(
    data = subset(dot_df, Sig == 0),
    color = "grey85",
    size  = 2.2,
    alpha = 0.6
  ) +
  # colored dots (significant)
  geom_point(
    data = subset(dot_df, Sig == 1),
    aes(color = Direction),
    size = 4.2
  ) +
  scale_color_manual(
    values = c(
      "Tumor_enriched" = "firebrick",
      "NAT_enriched"   = "steelblue",
      "NA"             = "grey40"
    )
  ) +
  theme_minimal(base_size = 12) +
  scale_x_discrete(expand = c(0.02, 0.02)) +
  scale_y_discrete(expand = c(0.05, 0.05)) +
  theme(
    panel.grid      = element_blank(),
    axis.text.x     = element_text(size = 10, angle = 45, hjust = 1, vjust = 1),
    axis.text.y     = element_text(size = 11),
    axis.title.x    = element_text(size = 12),
    axis.title.y    = element_text(size = 12),
    plot.title      = element_text(size = 15, face = "bold", hjust = 0.5),
    legend.title    = element_text(size = 11),
    legend.text     = element_text(size = 10),
    legend.key.size = unit(0.55, "lines")
  ) +
  labs(
    title = "Overlap of Significant Genera Across Tools Breast (prev10 Imputed)",
    x     = "DA Tool",
    y     = "Genus",
    color = "Direction"
  )

# Print to screen
print(p_dot)

# Save dot plot (PNG + PDF)
ggsave(
  filename = "dotplot_imp_prev10_imp.png",
  plot     = p_dot,
  width    = 8,
  height   = 6,
  dpi      = 300
)

ggsave(
  filename = "dotplot_prev10_imp.pdf",
  plot     = p_dot,
  width    = 8,
  height   = 6
)

cat("Dot plot saved as: dotplot_imp_prev10_imp.png and dotplot_pre10_imp.pdf\n")


## 5. UPSET PLOT (using same overlap_mat, same files) ----------

# Build lists of taxa per tool from overlap_mat
sig_list <- list(
  ALDEx2   = overlap_mat$taxon[overlap_mat$ALDEx2   == 1L] %>% unique(),
  ANCOMBC2 = overlap_mat$taxon[overlap_mat$ANCOMBC2 == 1L] %>% unique(),
  DESeq2   = overlap_mat$taxon[overlap_mat$DESeq2   == 1L] %>% unique(),
  MaAsLin2 = overlap_mat$taxon[overlap_mat$Maaslin2 == 1L] %>% unique()
)

# Report counts
cat("\nNumber of significant genera per tool (alpha =", alpha, "):\n")
cat("  ALDEx2   :", length(sig_list$ALDEx2),   "\n")
cat("  ANCOMBC2 :", length(sig_list$ANCOMBC2), "\n")
cat("  DESeq2   :", length(sig_list$DESeq2),   "\n")
cat("  MaAsLin2 :", length(sig_list$MaAsLin2), "\n\n")

# Plot UpSet on screen
upset(
  fromList(sig_list),
  order.by        = "freq",
  nsets           = 4,
  
  # COLORS MATCHING DOT PLOT:
  main.bar.color  = "steelblue",  # same blue as NAT_enriched
  sets.bar.color  = "firebrick",  # same red as Tumor_enriched
  
  point.size      = 3.5,
  line.size       = 1.5,
  text.scale      = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.2),
  mb.ratio        = c(0.6, 0.4),
  mainbar.y.label = "Intersection Size",
  sets.x.label    = "Number of significant genera (q < 0.05)"
)

# Save UpSet to PDF
pdf("upset_ovary_prev10_imp.pdf", width = 12, height = 7)

upset(
  fromList(sig_list),
  order.by        = "freq",
  nsets           = 4,
  main.bar.color  = "steelblue",
  sets.bar.color  = "firebrick",
  point.size      = 3.5,
  line.size       = 1.5,
  text.scale      = c(1.5, 1.3, 1.2, 1.2, 1.5, 1.2),
  mb.ratio        = c(0.6, 0.4),
  mainbar.y.label = "Intersection Size",
  sets.x.label    = "Number of significant genera (q < 0.05)"
)

dev.off()

cat("UpSet plot saved as: upset_lung_prev10_imp.pdf\n")
