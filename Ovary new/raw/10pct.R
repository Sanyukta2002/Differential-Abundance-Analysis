setwd("C:/Users/sanyu/OneDrive/Desktop/Ovary new/raw")
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(DESeq2)
  library(ALDEx2)
  library(ANCOMBC)
  library(Maaslin2)
})

## -------------------- SETUP ------------------------------

# FDR threshold and prevalence cutoff
alpha    <- 0.05
prev_cut <- 0.10   # 10% prevalence

# Create base output directory
base_dir <- "DA_10pct"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

# Reusable column names
sampleid <- "Sample_ID"
groupcol <- "Group"
patienid <- "Patient_ID"

## -------------------- LOAD ABUNDANCE ---------------------

# Load genus-level abundance table (integer counts)
abund <- read.csv(
  "abundanceint.csv",
  check.names      = FALSE,  # do not alter column names
  stringsAsFactors = FALSE   # keep strings as character, not factor
)

# Non-count columns (taxonomy + prevalence info)
taxa <- c(
  "genus",
  "Prevalence in DNA extraction/NTC negative controls",
  "Prevalence in Paraf. Controls",
  "domain",
  "phylum",
  "class",
  "order",
  "family" 
)

# Set genus column as rownames (rows = taxa / genus)
genus_col <- "genus"
rownames(abund) <- abund[[genus_col]]

# Keep only sample columns (remove taxonomy + prevalence columns)
# setdiff(A, B) = everything in A that is not in B
otu <- abund[, setdiff(colnames(abund), taxa), drop = FALSE]

# Ensure numeric counts (defensive, in case some columns were read as character)
otu[] <- lapply(otu, as.numeric)

# Convert the OTU data frame to a numeric matrix: taxa x samples
counts_mat <- as.matrix(otu)


## -------------------- LOAD METADATA ----------------------

meta_full <- read.csv(
  "clinical.csv",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

# Align metadata rows to columns of counts_mat using Sample_ID
# match() gives row indices of meta_full in the order of colnames(counts_mat)
meta_aligned <- meta_full[match(colnames(counts_mat), meta_full[[sampleid]]), ]

# Check alignment: Sample_ID in metadata must exactly match column order of counts_mat
stopifnot(all(meta_aligned[[sampleid]] == colnames(counts_mat)))

# Factor Group and Patient_ID for modeling
meta_aligned[[groupcol]] <- factor(meta_aligned[[groupcol]], levels = c("NAT", "Tumor"))
meta_aligned[[patienid]] <- factor(meta_aligned[[patienid]])


## ---------------- PREVALENCE FILTERING (10%) -------------

# Number of samples
n_samples <- ncol(counts_mat)

# Minimum number of non-zero samples required for a taxon to be kept
# e.g. 10% of 308 samples = 30.8 → ceiling = 31 samples
min_non_zero_samples <- ceiling(prev_cut * n_samples)

# For each taxon (row), count the number of samples where count > 0
keep_taxa <- rowSums(counts_mat > 0) >= min_non_zero_samples

# Filter count matrix to keep only taxa passing 10% prevalence
count_prev10 <- counts_mat[keep_taxa, , drop = FALSE]

cat("\nPrevalence filter (10%)\n")
cat("  Taxa before filter:", nrow(counts_mat), "\n")
cat("  Taxa after  10%   :", nrow(count_prev10), "\n")


## =========================================================
## 1) DESeq2 (paired NAT vs Tumor, 10% prev)
## =========================================================

# Output directory for DESeq2
out_dir_deseq <- file.path(base_dir, "DESEQ2")
dir.create(out_dir_deseq, showWarnings = FALSE, recursive = TRUE)

# Create DESeq2 dataset: design = ~ Patient_ID + Group (paired design)
dds <- DESeqDataSetFromMatrix(
  countData = count_prev10,
  colData   = meta_aligned,
  design    = ~ Patient_ID + Group
)

# Estimate size factors (poscounts recommended for sparse microbiome data)
dds <- estimateSizeFactors(dds, type = "poscounts")

# Run DESeq2 pipeline (dispersion estimation + model fitting)
dds <- DESeq(dds)

# Inspect result names
cat("\nDESeq2 result names:\n")
print(resultsNames(dds))

# Extract results for Tumor vs NAT (log2FC = Tumor / NAT)
res <- results(dds, contrast = c("Group", "Tumor", "NAT"))
summary(res)

# LFC shrinkage to stabilize log2FC, especially for low counts
if (requireNamespace("apeglm", quietly = TRUE)) {
  # Try to find the Tumor_vs_NAT coefficient name
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
      contrast = c("Group", "Tumor", "NAT"),
      type = "normal"
    )
  }
} else {
  message("Package 'apeglm' not available, using 'normal' shrinkage.")
  res_shrunk <- lfcShrink(
    dds,
    contrast = c("Group", "Tumor", "NAT"),
    type = "normal"
  )
}

# Convert to data frame and add genus column
res_deseq <- as.data.frame(res_shrunk)
res_deseq[[genus_col]] <- rownames(res_deseq)

# Add significance (padj < alpha) and direction (Tumor_enriched / NAT_enriched)
res_deseq <- res_deseq %>%
  mutate(
    significant = !is.na(padj) & padj < alpha,
    direction   = case_when(
      !significant       ~ "ns",
      log2FoldChange > 0 ~ "Tumor_enriched",  # higher in Tumor
      log2FoldChange < 0 ~ "NAT_enriched",    # higher in NAT
      TRUE               ~ "ns"
    )
  ) %>%
  arrange(padj)

# Summary counts
n_sig_deseq  <- sum(res_deseq$significant)
n_up_deseq   <- sum(res_deseq$direction == "Tumor_enriched")
n_down_deseq <- sum(res_deseq$direction == "NAT_enriched")

cat("\nDESeq2 (prev 10%)\n")
cat("  Significant (padj <", alpha, "):", n_sig_deseq, "\n")
cat("    Tumor_enriched:", n_up_deseq, "\n")
cat("    NAT_enriched  :", n_down_deseq, "\n")

# Write outputs
write_csv(
  res_deseq,
  file.path(out_dir_deseq, "DESeq2_prev10_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_deseq %>% filter(significant),
  file.path(out_dir_deseq, "DESeq2_prev10_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_deseq,
    n_Tumor_enriched = n_up_deseq,
    n_NAT_enriched   = n_down_deseq
  ),
  file.path(out_dir_deseq, "DESeq2_prev10_Tumor_vs_NAT_summary.csv")
)


## =========================================================
## 2) ALDEx2 (paired NAT vs Tumor, 10% prev, no internal filter)
## =========================================================

out_dir_aldex <- file.path(base_dir, "ALDEx2")
dir.create(out_dir_aldex, showWarnings = FALSE, recursive = TRUE)

# Order metadata by Patient_ID, then Group, so NAT/Tumor pairs are adjacent
meta_ord <- meta_aligned %>%
  arrange(.data[[patienid]], .data[[groupcol]])

# Reorder 10%-filtered count matrix columns to match this metadata order
counts_aldex <- count_prev10[, meta_ord[[sampleid]], drop = FALSE]
stopifnot(all(colnames(counts_aldex) == meta_ord[[sampleid]]))

# Condition vector (NAT / Tumor) in the matching order
conds <- as.character(meta_ord[[groupcol]])
cat("\nALDEx2 group counts (10% prev):\n")
print(table(conds))

# CLR transformation + Monte Carlo sampling
# denom = "all" → no internal feature filtering
set.seed(123)
x_clr <- aldex.clr(
  counts_aldex,
  conds,
  mc.samples = 128,
  denom      = "all",
  verbose    = TRUE
)

# Paired within-patient NAT vs Tumor test
x_tt <- aldex.ttest(
  x_clr,
  paired.test = TRUE,
  verbose     = TRUE
)

# Effect sizes (difference, dispersion)
x_eff <- aldex.effect(
  x_clr,
  verbose = TRUE
)

# Combine into one result table
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
    # Use minimum of the two FDR values as combined evidence
    p_min  = pmin(q_we, q_wi, na.rm = TRUE),
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

cat("\nALDEx2 (prev 10%, paired)\n")
cat("  Significant (min FDR <", alpha, "):", n_sig_aldex, "\n")
cat("    Tumor_enriched:", n_up_aldex, "\n")
cat("    NAT_enriched  :", n_down_aldex, "\n")

write_csv(
  res_aldex,
  file.path(out_dir_aldex, "ALDEx2_prev10_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_aldex %>% filter(significant),
  file.path(out_dir_aldex, "ALDEx2_prev10_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_aldex,
    n_Tumor_enriched = n_up_aldex,
    n_NAT_enriched   = n_down_aldex
  ),
  file.path(out_dir_aldex, "ALDEx2_prev10_Tumor_vs_NAT_summary.csv")
)


## =========================================================
## 3) ANCOM-BC2 (paired via random effect, 10% prev, no internal filter)
## =========================================================

out_dir_ancom <- file.path(base_dir, "ANCOMBC2")
dir.create(out_dir_ancom, showWarnings = FALSE, recursive = TRUE)

meta_ancom <- meta_aligned
rownames(meta_ancom) <- meta_ancom[[sampleid]]

cat("\nANCOM-BC2: taxa =", nrow(count_prev10),
    ", samples =", ncol(count_prev10), "\n")

set.seed(123)
out_ancom <- ancombc2(
  data          = count_prev10,
  taxa_are_rows = TRUE,
  meta_data     = meta_ancom,
  fix_formula   = groupcol,                 # Group (NAT vs Tumor)
  rand_formula  = paste0("(1 | ", patienid, ")"),  # random intercept per patient
  p_adj_method  = "BH",
  pseudo        = 0,
  pseudo_sens   = TRUE,
  prv_cut       = 0,                        # no internal prevalence filter
  lib_cut       = 0,                        # no internal library-size filter
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
    lfc         = lfc_GroupTumor,
    qval        = q_GroupTumor,
    diff        = diff_GroupTumor,
    diff_robust = diff_robust_GroupTumor
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

cat("\nANCOM-BC2 (prev 10%, paired)\n")
cat("  Significant (q <", alpha, "):", n_sig_ancom, "\n")
cat("    Tumor_enriched:", n_up_ancom, "\n")
cat("    NAT_enriched  :", n_down_ancom, "\n")
cat("  Robust significant:", n_sig_robust_ancom, "\n")
cat("    Tumor_enriched (robust):", n_up_robust_ancom, "\n")
cat("    NAT_enriched   (robust):", n_down_robust_ancom, "\n")

write_csv(
  res_ancom_tbl,
  file.path(out_dir_ancom, "ANCOMBC2_prev10_Tumor_vs_NAT_full.csv")
)

write_csv(
  res_ancom_tbl %>% filter(significant),
  file.path(out_dir_ancom, "ANCOMBC2_prev10_Tumor_vs_NAT_significant_q05.csv")
)

write_csv(
  res_ancom_tbl %>% filter(significant_robust),
  file.path(out_dir_ancom, "ANCOMBC2_prev10_Tumor_vs_NAT_significant_robust_q05.csv")
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
  file.path(out_dir_ancom, "ANCOMBC2_prev10_Tumor_vs_NAT_summary.csv")
)


## =========================================================
## 4) MaAsLin2 (paired via random effect, 10% prev, no internal filter)
## =========================================================

out_dir_maaslin <- file.path(base_dir, "Maaslin2")
dir.create(out_dir_maaslin, showWarnings = FALSE, recursive = TRUE)

meta_maaslin <- meta_aligned[, c(sampleid, groupcol, patienid)]
rownames(meta_maaslin) <- meta_maaslin[[sampleid]]

meta_maaslin[[groupcol]] <- factor(meta_maaslin[[groupcol]], levels = c("NAT", "Tumor"))
meta_maaslin[[patienid]] <- factor(meta_maaslin[[patienid]])

stopifnot(all(colnames(count_prev10) == rownames(meta_maaslin)))

fit_maaslin <- Maaslin2(
  input_data      = count_prev10,
  input_metadata  = meta_maaslin,
  output          = out_dir_maaslin,
  fixed_effects   = groupcol,
  random_effects  = patienid,
  normalization   = "TSS",
  transform       = "LOG",
  analysis_method = "LM",
  correction      = "BH",
  standardize     = TRUE,
  min_abundance   = 0,
  min_prevalence  = 0,
  min_variance    = 0,
  max_significance = 1
)

all_res <- read_tsv(
  file.path(out_dir_maaslin, "all_results.tsv"),
  show_col_types = FALSE
)

res_group <- all_res %>%
  filter(metadata == groupcol) %>%
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

cat("\nMaAsLin2 (prev 10%, paired)\n")
cat("  Significant (q <", alpha, "):", n_sig_maaslin, "\n")
cat("    Tumor_enriched:", n_up_maaslin, "\n")
cat("    NAT_enriched  :", n_down_maaslin, "\n")

write_csv(
  res_group,
  file.path(out_dir_maaslin, "Maaslin2_prev10_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_group %>% filter(significant),
  file.path(out_dir_maaslin, "Maaslin2_prev10_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_maaslin,
    n_Tumor_enriched = n_up_maaslin,
    n_NAT_enriched   = n_down_maaslin
  ),
  file.path(out_dir_maaslin, "Maaslin2_prev10_Tumor_vs_NAT_summary.csv")
)

cat("\nDONE: DESeq2, ALDEx2, ANCOM-BC2, and MaAsLin2 run with 10% prevalence filter on RAW integer data and paired design.\n")

