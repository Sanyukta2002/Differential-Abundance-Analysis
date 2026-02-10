setwd("C:/Users/sanyu/OneDrive/Desktop/Ovary new/raw")


library(tidyverse)   # dplyr + tidyr + readr + purrr etc.
library(exact2x2)    # for mcnemar.exact()

## 1. Load data -----------------------------------------------------------------

abundance <- read.csv(
  "abundanceint.csv",
  header = TRUE,
  check.names = FALSE   # keep sample IDs & prevalence column names exactly
)

clinical <- read.csv(
  "clinical.csv",
  header = TRUE,
  check.names = FALSE
)

cat("Abundance dimensions:", dim(abundance), "\n")
cat("Clinical dimensions :", dim(clinical),  "\n\n")

## 2. Identify taxonomy and sample columns -------------------------------------
## These must match EXACTLY what is in your abundance CSV

tax_cols <- c(
  "genus",
  "Prevalence in DNA extraction/NTC negative controls",
  "Prevalence in Paraf. Controls",
  "domain", "phylum", "class", "order", "family"
)

# Everything else in abundance is a sample column (Sample_ID)
sample_cols <- setdiff(names(abundance), tax_cols)

cat("Number of sample columns (should ≈ #rows in clinical):",
    length(sample_cols), "\n\n")

## 3. Reshape abundance: Sample_ID rows, genus columns -------------------------

abundance_long <- abundance %>%
  select(all_of(c("genus", sample_cols))) %>%
  pivot_longer(
    cols      = all_of(sample_cols),
    names_to  = "Sample_ID",
    values_to = "abundance"
  )

abundance_by_sample <- abundance_long %>%
  pivot_wider(
    names_from  = genus,
    values_from = abundance
  )
# Now: one row per Sample_ID, one column per genus

## 4. Merge with clinical metadata ---------------------------------------------

# clinical has:
#   "Sample_ID", "Sample_ID (WIS)", "Group", "Material", "Center",
#   "DNA Extraction batch", "PCR batch", "Sequencing lane batch",
#   "Patient_ID", "Age", "Gender", "Sample_ID_original"

merged <- clinical %>%
  inner_join(abundance_by_sample, by = "Sample_ID")

cat("Merged dimensions (samples with both metadata + abundance):",
    dim(merged), "\n")
cat("Unique groups in merged data:", unique(merged$Group), "\n")
cat("Number of unique patients:", length(unique(merged$Patient_ID)), "\n\n")

## 5. Identify paired patients (Tumor + NAT present) ---------------------------

tumor_group  <- "Tumor"
normal_group <- "NAT"

paired_patients <- merged %>%
  group_by(Patient_ID) %>%
  summarise(
    has_tumor  = any(Group == tumor_group),
    has_normal = any(Group == normal_group),
    .groups    = "drop"
  ) %>%
  filter(has_tumor & has_normal) %>%
  pull(Patient_ID)

cat("Number of paired patients (Tumor + NAT):", length(paired_patients), "\n")

if (length(paired_patients) == 0) {
  stop("No paired patients found! Check Group labels / data.")
}

## 6. Filter to paired samples and verify exact pairing ------------------------

merged_paired <- merged %>%
  filter(Patient_ID %in% paired_patients)

tumor_paired <- merged_paired %>%
  filter(Group == tumor_group) %>%
  arrange(Patient_ID)

normal_paired <- merged_paired %>%
  filter(Group == normal_group) %>%
  arrange(Patient_ID)

cat("Tumor samples :", nrow(tumor_paired), "\n")
cat("Normal samples:", nrow(normal_paired), "\n")

if (nrow(tumor_paired) != nrow(normal_paired)) {
  stop("Unequal number of tumor and normal samples!")
}

if (!all(tumor_paired$Patient_ID == normal_paired$Patient_ID)) {
  stop("Patient IDs don't match between tumor and normal!")
}

cat("✓ Pairing verified! Testing", nrow(tumor_paired), "pairs\n\n")

## 7. Identify taxa (genus) columns -------------------------------------------

# All columns in merged_paired that are NOT in clinical are taxa
taxa_columns <- setdiff(names(merged_paired), names(clinical))

cat("Number of taxa (genus columns):", length(taxa_columns), "\n\n")

## 8. Convert abundance to presence/absence (0/1) ------------------------------

threshold <- 0  # > 0 means "present"

tumor_binary  <- tumor_paired[, taxa_columns, drop = FALSE]
normal_binary <- normal_paired[, taxa_columns, drop = FALSE]

tumor_binary[]  <- (tumor_binary  > threshold) * 1
normal_binary[] <- (normal_binary > threshold) * 1

## 9. McNemar test function for one taxon -------------------------------------

mcnemar_test_taxon <- function(normal_vals, tumor_vals, taxon_name) {
  # normal_vals, tumor_vals: numeric 0/1 vectors, same length, paired by patient
  
  n00 <- sum(normal_vals == 0 & tumor_vals == 0)  # absent in both
  n01 <- sum(normal_vals == 0 & tumor_vals == 1)  # absent in NAT, present in Tumor
  n10 <- sum(normal_vals == 1 & tumor_vals == 0)  # present in NAT, absent in Tumor
  n11 <- sum(normal_vals == 1 & tumor_vals == 1)  # present in both
  
  # If no discordant pairs (n01 + n10 == 0), nothing to test
  if ((n01 + n10) == 0) return(NULL)
  
  # 2 x 2 contingency table: rows = Normal (0,1), cols = Tumor (0,1)
  contingency_table <- matrix(
    c(n00, n01,
      n10, n11),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Normal = c("0", "1"),
      Tumor  = c("0", "1")
    )
  )
  
  # Use exact McNemar for small discordant counts, chi-square otherwise
  if (n01 + n10 < 25) {
    test_result <- mcnemar.exact(contingency_table)
    statistic   <- NA_real_
  } else {
    test_result <- mcnemar.test(contingency_table, correct = TRUE)
    statistic   <- as.numeric(test_result$statistic)
  }
  
  prev_normal <- sum(normal_vals) / length(normal_vals)
  prev_tumor  <- sum(tumor_vals)  / length(tumor_vals)
  
  data.frame(
    Taxon                           = taxon_name,
    n00_absent_both                 = n00,
    n01_absent_normal_present_tumor = n01,  # gain in Tumor
    n10_present_normal_absent_tumor = n10,  # loss in Tumor
    n11_present_both                = n11,
    Prevalence_Normal               = round(prev_normal, 4),
    Prevalence_Tumor                = round(prev_tumor, 4),
    Prevalence_Diff                 = round(prev_tumor - prev_normal, 4),
    statistic                       = statistic,
    pvalue                          = test_result$p.value,
    stringsAsFactors                = FALSE
  )
}

## 10. Run McNemar tests for all taxa -----------------------------------------

cat("Performing McNemar tests for all taxa...\n")

results_list <- list()

for (i in seq_along(taxa_columns)) {
  taxon <- taxa_columns[i]
  
  normal_vals <- as.numeric(normal_binary[[taxon]])
  tumor_vals  <- as.numeric(tumor_binary[[taxon]])
  
  res <- mcnemar_test_taxon(
    normal_vals = normal_vals,
    tumor_vals  = tumor_vals,
    taxon_name  = taxon
  )
  
  if (!is.null(res)) {
    results_list[[length(results_list) + 1L]] <- res
  }
  
  if (i %% 100 == 0) {
    cat("  Tested", i, "of", length(taxa_columns), "taxa\n")
  }
}

cat("✓ Tests complete! Taxa with ≥1 discordant pair:",
    length(results_list), "\n\n")

if (length(results_list) == 0) {
  
  cat("No taxa with discordant pairs found. Nothing to test.\n")
  
} else {
  
  ## 11. Combine results, add FDR & direction ---------------------------------
  
  results_df <- bind_rows(results_list)
  
  results_df <- results_df %>%
    mutate(
      pvalue_adjusted = p.adjust(pvalue, method = "fdr"),
      direction = case_when(
        Prevalence_Diff >  0 ~ "Tumor_enriched",
        Prevalence_Diff <  0 ~ "NAT_enriched",
        TRUE                 ~ "No_prevalence_change"
      )
    ) %>%
    arrange(pvalue)
  
  ## 12. Summary + save to CSV -----------------------------------------------
  
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("MCNEMAR TEST RESULTS FOR PAIRED TUMOR vs NAT (presence/absence)\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  cat("Total taxa tested (with discordant pairs):", nrow(results_df), "\n")
  cat("Significant taxa (p < 0.05):",
      sum(results_df$pvalue < 0.05), "\n")
  cat("Significant taxa (FDR < 0.05):",
      sum(results_df$pvalue_adjusted < 0.05), "\n\n")
  
  cat(paste(rep("-", 80), collapse = ""), "\n")
  cat("TOP 20 MOST SIGNIFICANT TAXA (by raw p-value):\n")
  cat(paste(rep("-", 80), collapse = ""), "\n")
  print(head(results_df, 20), row.names = FALSE)
  
  out_file <- "mcnemar_test_results.csv"
  write.csv(results_df, out_file, row.names = FALSE)
  
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat("✓ Full results saved to:", out_file, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n\n")
  
  cat("INTERPRETATION NOTES:\n")
  cat("- Positive Prevalence_Diff: more prevalent in Tumor\n")
  cat("- Negative Prevalence_Diff: more prevalent in NAT\n")
  cat("- n01 (0→1): absent in NAT, present in Tumor (Tumor gain)\n")
  cat("- n10 (1→0): present in NAT, absent in Tumor (Tumor loss)\n")
}

