setwd("C:/Users/sanyu/OneDrive/Desktop/Lung_new/imputed")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(DESeq2)
  library(ALDEx2)
  library(ANCOMBC)
  library(Maaslin2)
})

## --------------------- SETUP -----------------------------

alpha    <- 0.05        # FDR / q / padj threshold
prev_cut <- 0.10        # 10% prevalence
base_dir <- "differential_abundance_prev10_imputed"

dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# Column name assumptions in metadata
sample_id_col <- "Sample_ID_clean"   # matches OTU column names
group_col     <- "Group"             # NAT / Tumor
patient_col   <- "Patient_ID"        # pairing / subject ID

## ------------------ LOAD OTU & META ----------------------

# 1) OTU table (taxa x samples, first column = genus)
otu <- read.csv(
  "otu_imputed_origlibsize_genus_lung_from_paired456.csv",
  check.names      = FALSE,  # keep original sample IDs
  stringsAsFactors = FALSE
)

# First column is the taxon name (e.g., "genus")
taxa_col <- colnames(otu)[1]       # should be "genus" or similar
rownames(otu) <- otu[[taxa_col]]   # set rownames = taxon names
otu[[taxa_col]] <- NULL            # remove taxon column from data matrix

# Ensure numeric (mbImpute outputs can be non-integer; we keep numeric here)
otu[] <- lapply(otu, as.numeric)

counts_mat <- as.matrix(otu)       # taxa (rows) x samples (columns)

cat("OTU table dims (taxa x samples):",
    nrow(counts_mat), "x", ncol(counts_mat), "\n")

# 2) Metadata
meta_full <- read.csv(
  "meta_mbimpute_genus_lung_paired.csv",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

## --------------- ALIGN + RESTRICT TO PAIRED --------------

# Sanity check: columns exist
stopifnot(sample_id_col %in% colnames(meta_full))
stopifnot(group_col     %in% colnames(meta_full))
stopifnot(patient_col   %in% colnames(meta_full))

# Keep only samples that appear in the imputed count matrix
meta <- meta_full %>%
  filter(.data[[sample_id_col]] %in% colnames(counts_mat)) %>%
  distinct(.data[[sample_id_col]], .keep_all = TRUE)

# Keep only NAT / Tumor samples
meta <- meta %>%
  filter(.data[[group_col]] %in% c("NAT", "Tumor"))

# Identify patients with BOTH NAT and Tumor
tab_pair   <- table(meta[[patient_col]], meta[[group_col]])
paired_ids <- names(which(rowSums(tab_pair > 0) == 2))

cat("Patients with both NAT & Tumor:", length(paired_ids), "\n")

# Restrict metadata to paired patients
meta_paired <- meta %>%
  filter(.data[[patient_col]] %in% paired_ids)

# Align metadata rows to columns of the count matrix
common_samples <- intersect(colnames(counts_mat), meta_paired[[sample_id_col]])

meta_paired   <- meta_paired[match(common_samples, meta_paired[[sample_id_col]]), ]
counts_paired <- counts_mat[, common_samples, drop = FALSE]

stopifnot(all(colnames(counts_paired) == meta_paired[[sample_id_col]]))

# Factors for modeling
meta_paired[[group_col]]   <- factor(meta_paired[[group_col]], levels = c("NAT", "Tumor"))
meta_paired[[patient_col]] <- factor(meta_paired[[patient_col]])

cat("After pairing & alignment: taxa =",
    nrow(counts_paired), ", samples =", ncol(counts_paired), "\n")

## ------------------- 10% PREVALENCE FILTER ---------------

n_samples <- ncol(counts_paired)
min_nonzero_samples <- ceiling(prev_cut * n_samples)   # taxa must be >0 in at least this many samples

keep_taxa <- rowSums(counts_paired > 0) >= min_nonzero_samples
counts_prev10 <- counts_paired[keep_taxa, , drop = FALSE]

cat("10% prevalence filter:",
    "keeping", sum(keep_taxa), "of", nrow(counts_paired), "taxa\n")

## =========================================================
## 1) DESeq2 (paired, 10% prev, imputed)
## =========================================================

out_dir_deseq <- file.path(base_dir, "DESeq2")
dir.create(out_dir_deseq, showWarnings = FALSE, recursive = TRUE)

# DESeq2 requires integer counts; round imputed values
counts_deseq <- round(counts_prev10)

dds <- DESeqDataSetFromMatrix(
  countData = counts_deseq,
  colData   = meta_paired,
  design    = as.formula(paste("~", patient_col, "+", group_col))  # ~ Patient_ID + Group
)

# Size factor estimation (poscounts recommended for sparse data)
dds <- estimateSizeFactors(dds, type = "poscounts")

# Run DESeq2 (dispersion + model fit)
dds <- DESeq(dds)

cat("\nDESeq2 result names:\n")
print(resultsNames(dds))

# Tumor vs NAT (log2FC = Tumor / NAT)
res <- results(dds, contrast = c(group_col, "Tumor", "NAT"))
summary(res)

# LFC shrinkage (prefer apeglm if available)
if (requireNamespace("apeglm", quietly = TRUE)) {
  coef_name <- grep("Tumor_vs_NAT", resultsNames(dds), value = TRUE)
  if (length(coef_name) == 1) {
    res_shrunk <- lfcShrink(
      dds,
      coef = coef_name,
      type = "apeglm"
    )
  } else {
    message("Could not uniquely identify Tumor_vs_NAT coef; using 'normal' shrinkage.")
    res_shrunk <- lfcShrink(
      dds,
      contrast = c(group_col, "Tumor", "NAT"),
      type = "normal"
    )
  }
} else {
  message("Package 'apeglm' not available, using 'normal' shrinkage.")
  res_shrunk <- lfcShrink(
    dds,
    contrast = c(group_col, "Tumor", "NAT"),
    type = "normal"
  )
}

res_deseq <- as.data.frame(res_shrunk)
res_deseq[[taxa_col]] <- rownames(res_deseq)   # add genus column

res_deseq <- res_deseq %>%
  mutate(
    significant = !is.na(padj) & padj < alpha,
    direction   = case_when(
      !significant       ~ "ns",
      log2FoldChange > 0 ~ "Tumor_enriched",
      log2FoldChange < 0 ~ "NAT_enriched",
      TRUE               ~ "ns"
    )
  ) %>%
  arrange(padj)

n_sig_deseq  <- sum(res_deseq$significant)
n_up_deseq   <- sum(res_deseq$direction == "Tumor_enriched")
n_down_deseq <- sum(res_deseq$direction == "NAT_enriched")

cat("\nDESeq2 (imputed, prev 10%)\n")
cat("  Significant (padj <", alpha, "):", n_sig_deseq, "\n")
cat("    Tumor_enriched:", n_up_deseq, "\n")
cat("    NAT_enriched  :", n_down_deseq, "\n")

write_csv(
  res_deseq,
  file.path(out_dir_deseq, "DESeq2_prev10_imputed_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_deseq %>% filter(significant),
  file.path(out_dir_deseq, "DESeq2_prev10_imputed_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_deseq,
    n_Tumor_enriched = n_up_deseq,
    n_NAT_enriched   = n_down_deseq
  ),
  file.path(out_dir_deseq, "DESeq2_prev10_imputed_Tumor_vs_NAT_summary.csv")
)

## =========================================================
## 2) ALDEx2 (paired, 10% prev, imputed)
## =========================================================

out_dir_aldex <- file.path(base_dir, "ALDEx2")
dir.create(out_dir_aldex, showWarnings = FALSE, recursive = TRUE)

# Order samples by Patient_ID, then Group (NAT first, then Tumor)
meta_ord <- meta_paired %>%
  arrange(.data[[patient_col]], .data[[group_col]])

# Reorder count matrix to match metadata order
counts_aldex <- counts_prev10[, meta_ord[[sample_id_col]], drop = FALSE]
stopifnot(all(colnames(counts_aldex) == meta_ord[[sample_id_col]]))

# Condition vector (NAT / Tumor) in correct order
conds <- as.character(meta_ord[[group_col]])
cat("\nALDEx2 group counts (imputed, prev 10%):\n")
print(table(conds))

set.seed(123)
# CLR transformation + Monte Carlo sampling; denom = "all" → no internal filtering
x_clr <- aldex.clr(
  counts_aldex,
  conds,
  mc.samples = 128,
  denom      = "all",
  verbose    = TRUE
)

# Paired test (within-patient differences)
x_tt <- aldex.ttest(
  x_clr,
  paired.test = TRUE,
  verbose     = TRUE
)

# Effect sizes
x_eff <- aldex.effect(
  x_clr,
  verbose = TRUE
)

res_aldex <- data.frame(
  taxon = rownames(x_tt),
  x_tt,
  x_eff,
  check.names = FALSE
)

res_aldex <- res_aldex %>%
  mutate(
    p_we   = we.ep,
    q_we   = we.eBH,
    p_wi   = wi.ep,
    q_wi   = wi.eBH,
    p_min  = pmin(q_we, q_wi, na.rm = TRUE),   # combined minimum FDR
    significant = !is.na(p_min) & p_min < alpha,
    direction   = case_when(
      !significant ~ "ns",
      diff.btw > 0 ~ "Tumor_enriched",
      diff.btw < 0 ~ "NAT_enriched",
      TRUE         ~ "ns"
    )
  ) %>%
  arrange(p_min)

n_sig_aldex  <- sum(res_aldex$significant)
n_up_aldex   <- sum(res_aldex$direction == "Tumor_enriched")
n_down_aldex <- sum(res_aldex$direction == "NAT_enriched")

cat("\nALDEx2 (imputed, prev 10%, paired)\n")
cat("  Significant (min FDR <", alpha, "):", n_sig_aldex, "\n")
cat("    Tumor_enriched:", n_up_aldex, "\n")
cat("    NAT_enriched  :", n_down_aldex, "\n")

write_csv(
  res_aldex,
  file.path(out_dir_aldex, "ALDEx2_prev10_imputed_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_aldex %>% filter(significant),
  file.path(out_dir_aldex, "ALDEx2_prev10_imputed_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_aldex,
    n_Tumor_enriched = n_up_aldex,
    n_NAT_enriched   = n_down_aldex
  ),
  file.path(out_dir_aldex, "ALDEx2_prev10_imputed_Tumor_vs_NAT_summary.csv")
)

## =========================================================
## 3) ANCOM-BC2 (paired via random effect, 10% prev, imputed)
## =========================================================

out_dir_ancom <- file.path(base_dir, "ANCOMBC2")
dir.create(out_dir_ancom, showWarnings = FALSE, recursive = TRUE)

meta_ancom <- meta_paired
rownames(meta_ancom) <- meta_ancom[[sample_id_col]]

cat("\nANCOM-BC2 (imputed): taxa =", nrow(counts_prev10),
    ", samples =", ncol(counts_prev10), "\n")

set.seed(123)
out_ancom <- ancombc2(
  data          = counts_prev10,
  taxa_are_rows = TRUE,
  meta_data     = meta_ancom,
  fix_formula   = group_col,                        # "Group"
  rand_formula  = paste0("(1 | ", patient_col, ")"),# "(1 | Patient_ID)"
  p_adj_method  = "BH",
  pseudo        = 0,
  pseudo_sens   = TRUE,   # sensitivity analysis, not a prevalence filter
  prv_cut       = 0,      # no internal prevalence filter
  lib_cut       = 0,      # no internal library-size filter
  s0_perc       = 0.05,
  group         = NULL,
  struc_zero    = FALSE,
  neg_lb        = FALSE,
  alpha         = alpha,
  n_cl          = 1,
  verbose       = TRUE,
  global        = FALSE,
  pairwise      = FALSE,
  dunnet        = FALSE,
  trend         = FALSE
)

res_ancom <- out_ancom$res

res_ancom_tbl <- res_ancom %>%
  as.data.frame() %>%
  mutate(
    lfc         = lfc_GroupTumor,          # log2FC Tumor vs NAT
    qval        = q_GroupTumor,            # BH-adjusted q-value
    diff        = diff_GroupTumor,         # ANCOM-BC2 "significant" flag
    diff_robust = diff_robust_GroupTumor   # robust flag from sensitivity analysis
  ) %>%
  mutate(
    significant        = !is.na(qval) & qval < alpha,
    significant_robust = diff_robust,
    direction          = case_when(
      !significant ~ "ns",
      lfc > 0      ~ "Tumor_enriched",
      lfc < 0      ~ "NAT_enriched",
      TRUE         ~ "ns"
    ),
    direction_robust   = case_when(
      !significant_robust ~ "ns",
      lfc > 0             ~ "Tumor_enriched",
      lfc < 0             ~ "NAT_enriched",
      TRUE                ~ "ns"
    )
  ) %>%
  arrange(qval)

n_sig_ancom         <- sum(res_ancom_tbl$significant)
n_sig_robust_ancom  <- sum(res_ancom_tbl$significant_robust)
n_up_ancom          <- sum(res_ancom_tbl$direction == "Tumor_enriched")
n_down_ancom        <- sum(res_ancom_tbl$direction == "NAT_enriched")
n_up_robust_ancom   <- sum(res_ancom_tbl$direction_robust == "Tumor_enriched")
n_down_robust_ancom <- sum(res_ancom_tbl$direction_robust == "NAT_enriched")

cat("\nANCOM-BC2 (imputed, prev 10%, paired)\n")
cat("  Significant (q <", alpha, "):", n_sig_ancom, "\n")
cat("    Tumor_enriched:", n_up_ancom, "\n")
cat("    NAT_enriched  :", n_down_ancom, "\n")
cat("  Robust significant:", n_sig_robust_ancom, "\n")
cat("    Tumor_enriched (robust):", n_up_robust_ancom, "\n")
cat("    NAT_enriched   (robust):", n_down_robust_ancom, "\n")

write_csv(
  res_ancom_tbl,
  file.path(out_dir_ancom, "ANCOMBC2_prev10_imputed_Tumor_vs_NAT_full.csv")
)

write_csv(
  res_ancom_tbl %>% filter(significant),
  file.path(out_dir_ancom, "ANCOMBC2_prev10_imputed_Tumor_vs_NAT_significant_q05.csv")
)

write_csv(
  res_ancom_tbl %>% filter(significant_robust),
  file.path(out_dir_ancom, "ANCOMBC2_prev10_imputed_Tumor_vs_NAT_significant_robust_q05.csv")
)

write_csv(
  tibble(
    alpha                   = alpha,
    n_significant           = n_sig_ancom,
    n_Tumor_enriched        = n_up_ancom,
    n_NAT_enriched          = n_down_ancom,
    n_significant_robust    = n_sig_robust_ancom,
    n_Tumor_enriched_robust = n_up_robust_ancom,
    n_NAT_enriched_robust   = n_down_robust_ancom
  ),
  file.path(out_dir_ancom, "ANCOMBC2_prev10_imputed_Tumor_vs_NAT_summary.csv")
)

## =========================================================
## 4) MaAsLin2 (paired via random effect, 10% prev, imputed)
## =========================================================

out_dir_maaslin <- file.path(base_dir, "Maaslin2")
dir.create(out_dir_maaslin, showWarnings = FALSE, recursive = TRUE)

# Minimal metadata for MaAsLin2
meta_maaslin <- meta_paired[, c(sample_id_col, group_col, patient_col)]
rownames(meta_maaslin) <- meta_maaslin[[sample_id_col]]

meta_maaslin[[group_col]]   <- factor(meta_maaslin[[group_col]], levels = c("NAT", "Tumor"))
meta_maaslin[[patient_col]] <- factor(meta_maaslin[[patient_col]])

stopifnot(all(colnames(counts_prev10) == rownames(meta_maaslin)))

fit_maaslin <- Maaslin2(
  input_data      = counts_prev10,
  input_metadata  = meta_maaslin,
  output          = out_dir_maaslin,
  fixed_effects   = group_col,      # Group (NAT vs Tumor)
  random_effects  = patient_col,    # random intercept for Patient_ID
  normalization   = "TSS",
  transform       = "LOG",
  analysis_method = "LM",
  correction      = "BH",
  standardize     = TRUE,
  min_abundance   = 0,
  min_prevalence  = 0,              # no internal prevalence filter
  min_variance    = 0,
  max_significance = 1
)

all_res <- read_tsv(
  file.path(out_dir_maaslin, "all_results.tsv"),
  show_col_types = FALSE
)

res_group <- all_res %>%
  filter(metadata == group_col) %>%  # only Group effect
  mutate(
    significant = !is.na(qval) & qval < alpha,
    direction   = case_when(
      !significant ~ "ns",
      coef > 0     ~ "Tumor_enriched",
      coef < 0     ~ "NAT_enriched",
      TRUE         ~ "ns"
    )
  ) %>%
  arrange(qval)

n_sig_maaslin  <- sum(res_group$significant)
n_up_maaslin   <- sum(res_group$direction == "Tumor_enriched")
n_down_maaslin <- sum(res_group$direction == "NAT_enriched")

cat("\nMaAsLin2 (imputed, prev 10%, paired)\n")
cat("  Significant (q <", alpha, "):", n_sig_maaslin, "\n")
cat("    Tumor_enriched:", n_up_maaslin, "\n")
cat("    NAT_enriched  :", n_down_maaslin, "\n")

write_csv(
  res_group,
  file.path(out_dir_maaslin, "Maaslin2_prev10_imputed_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_group %>% filter(significant),
  file.path(out_dir_maaslin, "Maaslin2_prev10_imputed_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_maaslin,
    n_Tumor_enriched = n_up_maaslin,
    n_NAT_enriched   = n_down_maaslin
  ),
  file.path(out_dir_maaslin, "Maaslin2_prev10_imputed_Tumor_vs_NAT_summary.csv")
)

cat("\nDONE: all four DA tools run on mbImpute-imputed data with 10% prevalence filter and paired design.\n")

