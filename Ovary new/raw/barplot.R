setwd("C:/Users/sanyu/OneDrive/Desktop/Ovary new/raw")

library(tidyverse)
X11()

# 1. Read McNemar results for ovary cohort
df <- read.csv("mcnemar_test_results.csv")

# 2. Take top 10 by absolute prevalence difference
top10_ovary <- df %>%
  arrange(desc(abs(Prevalence_Diff))) %>%
  slice(1:10) %>%
  mutate(
    direction = factor(direction,
                       levels = c("Tumor_enriched", "NAT_enriched"))
  )

# 3. Plot: bars = |Tumor − NAT|, color = which side is enriched
p_ovary <- ggplot(
  top10_ovary,
  aes(
    x = reorder(Taxon, abs(Prevalence_Diff)),
    y = abs(Prevalence_Diff),
    fill = direction
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title    = "Top 10 Bacterial Genera by Prevalence Difference in Ovary Cohort",
    subtitle = "Bars show absolute prevalence difference |Tumor − NAT|; color indicates direction of enrichment",
    x        = "Genus",
    y        = "Absolute prevalence difference |Tumor − NAT|",
    fill     = "Direction of enrichment"
  ) +
  scale_fill_manual(
    values = c(
      "Tumor_enriched" = "#d62728",  # red
      "NAT_enriched"   = "#1f77b4"   # blue
    )
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.y   = element_text(size = 12)
  )

# 4. Print the plot
print(p_ovary)

# 5. (Optional) Save it
ggsave("ovary_top10_prevalence_diff.png",
       plot = p_ovary,
       width = 8, height = 6, dpi = 300)
