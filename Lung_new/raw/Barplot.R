setwd("C:/Users/sanyu/OneDrive/Desktop/Lung_new/raw")
X11()
library(tidyverse)

df <- read.csv("mcnemar_test_results.csv")

# Top 10 taxa by absolute prevalence difference
top10 <- df %>%
  arrange(desc(abs(Prevalence_Diff))) %>%
  slice(1:10) %>%
  mutate(
    direction = factor(direction,
                       levels = c("Tumor_enriched", "NAT_enriched"))
  )

ggplot(
  top10,
  aes(
    x = reorder(Taxon, abs(Prevalence_Diff)),
    y = abs(Prevalence_Diff),
    fill = direction
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title    = "Top 10 Bacterial Genera by Prevalence Difference in Lung Cohort",
    subtitle = "Bars show absolute prevalence difference; color indicates whether Tumor or NAT is enriched",
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
