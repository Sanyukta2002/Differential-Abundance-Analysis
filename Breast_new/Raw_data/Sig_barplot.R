setwd("C:/Users/sanyu/OneDrive/Desktop/Breast_new/Raw_data")
library(tidyverse)

df <- read.csv("mcnemar_test_results.csv")

# Top 10 taxa by absolute prevalence difference
top10 <- df %>%
  arrange(desc(abs(Prevalence_Diff))) %>%
  slice(1:10) %>%
  mutate(
    direction = factor(direction,
                       levels = c("NAT_enriched", "Tumor_enriched"))
  )

ggplot(
  top10,
  aes(
    x = reorder(Taxon, Prevalence_Diff),
    y = Prevalence_Diff,
    fill = direction
  )
) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_col() +
  coord_flip() +
  labs(
    title = "Top 10 Bacterial Genera by Prevalence Difference (Tumor − NAT)",
    subtitle = "Breast cancer cohort • Bars right of 0: Tumor-enriched • Bars left of 0: NAT-enriched",
    x = "Genus",
    y = "Prevalence Difference (Tumor − NAT)"
  ) +
  scale_fill_manual(
    values = c(
      "Tumor_enriched" = "#d62728",  # red
      "NAT_enriched"   = "#1f77b4"   # blue
    ),
    name = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    axis.text.y   = element_text(size = 12)
  )

X11()