setwd("C:/Users/sanyu/OneDrive/Desktop/Lung_new/raw")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(DESeq2)
  library(ALDEx2)
  library(ANCOMBC)
  library(Maaslin2)
})


#Setup fdr and filtering values
alpha <- 0.05
prev_cut <- 0.05

#Create directory to save outputs
base_dir <- "DA_5pct"
dir.create( base_dir , showWarnings = FALSE , recursive = TRUE) # create all parent dirs if they dont exist


#Set up reusable column names:
sampleid <- "Sample_ID"
groupcol <- "Group"
patienid <- "Patient_ID"

#Load Abundance table
abund <- read.csv("abundanceint.csv", check.names = FALSE, stringsAsFactors = FALSE) ##checknames false = does not let R alter column name , stringfactor = keep all charatcers as strings, not factors

#Exclude the non-count columns
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


#Set genus column as rownames
#take the genus column
genus_col <- "genus"
rownames(abund) <- abund[[genus_col]]

#Make a table of only abundance data
otu <- abund[ , setdiff(colnames(abund), taxa) , drop =FALSE] #everything from A that is not in B , drop = false forces result to be a df

#convert the otu(dataframe) to a matrix
counts_mat <- as.matrix(otu)

#Load the metadata
meta_full <- read.csv("clinical.csv" , check.names = FALSE , stringsAsFactors = FALSE)

#Align clinic file rows to columns of count_mat
meta_aligned <- meta_full[match(colnames(counts_mat) , meta_full[[sampleid]]) , ] #comma here means take rows as indicated in match logic and take all columns from clincial file

#Check if aligned
stopifnot(all(meta_aligned[[sampleid]] == colnames(counts_mat)))

#Factor Group and PatientID
meta_aligned[[groupcol]] <- factor(meta_aligned[[groupcol]] , levels = c("NAT" , "Tumor"))
meta_aligned[[patienid]] <- factor(meta_aligned[[patienid]])

#--------Prevalence Filtering----------
n_samples <- ncol(counts_mat)
n_samples

min_non_zero_samples <- ceiling(prev_cut * n_samples) #ceiling rounds upto nearest integer

keep_taxa <- rowSums(counts_mat > 0 ) >= min_non_zero_samples
count_prev5 <- counts_mat[keep_taxa , , drop = FALSE]


#-----------DA - ANALYSIS-------

#Create output directory to store output
out_dir_deseq <- file.path(base_dir , "DESEQ2")
dir.create(out_dir_deseq , showWarnings = FALSE, recursive = TRUE )

#Create a dds object
dds <- DESeqDataSetFromMatrix(
  countData = count_prev5,
  colData = meta_aligned,
  design = ~ Patient_ID + Group )

#
dds <- estimateSizeFactors(dds , type = "poscounts")
dds <- DESeq(dds)

# check what DESeq2 thinks the coefficients are
cat("\nDESeq2 result names:\n")
print(resultsNames(dds))

# Extract results for Tumor vs NAT (log2FC = Tumor / NAT)
res <- results(dds, contrast = c("Group", "Tumor", "NAT"))
summary(res)

# LFC shrinkage (stabilize log2FC, especially for low counts)
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

# Annotate significance and direction
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

# Quick summary counts
n_sig_deseq  <- sum(res_deseq$significant)
n_up_deseq   <- sum(res_deseq$direction == "Tumor_enriched")
n_down_deseq <- sum(res_deseq$direction == "NAT_enriched")

cat("\nDESeq2 (prev 5%)\n")
cat("  Significant (padj <", alpha, "):", n_sig_deseq, "\n")
cat("    Tumor_enriched:", n_up_deseq, "\n")
cat("    NAT_enriched  :", n_down_deseq, "\n")

# Write outputs: all, significant only, and a tiny summary
write_csv(
  res_deseq,
  file.path(out_dir_deseq, "DESeq2_prev5_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_deseq %>% filter(significant),
  file.path(out_dir_deseq, "DESeq2_prev5_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_deseq,
    n_Tumor_enriched = n_up_deseq,
    n_NAT_enriched   = n_down_deseq
  ),
  file.path(out_dir_deseq, "DESeq2_prev5_Tumor_vs_NAT_summary.csv")
)

## -ALDEX------
## =========================================================
## 2) ALDEx2 (paired NAT vs Tumor, 5% prev, no internal filter)
## =========================================================

# Create output directory for ALDEx2 results
out_dir_aldex <- file.path(base_dir, "ALDEx2")
dir.create(out_dir_aldex, showWarnings = FALSE, recursive = TRUE)

# Order metadata by Patient_ID, then Group
# (so NAT / Tumor from same patient are adjacent; needed for paired ALDEx2)
meta_ord <- meta_aligned %>%
  arrange(.data[[patienid]], .data[[groupcol]])

# Reorder count matrix columns to match this ordered metadata
counts_aldex <- count_prev5[, meta_ord[[sampleid]], drop = FALSE]
stopifnot(all(colnames(counts_aldex) == meta_ord[[sampleid]]))

# Condition vector for ALDEx2 (NAT / Tumor per sample, in the correct order)
conds <- as.character(meta_ord[[groupcol]])
cat("\nALDEx2 group counts:\n")
print(table(conds))  # sanity check: counts of NAT vs Tumor

# CLR transformation + Monte Carlo sampling
# denom = "all" → use all features as denominator (no internal prevalence filtering)
set.seed(123)
x_clr <- aldex.clr(
  counts_aldex,
  conds,
  mc.samples = 128,   # number of Monte Carlo Dirichlet instances
  denom      = "all",
  verbose    = TRUE
)

# Paired t-test (within-patient NAT vs Tumor)
x_tt <- aldex.ttest(
  x_clr,
  paired.test = TRUE,  # tells ALDEx2 this is a paired design
  verbose     = TRUE
)

# Effect size estimates (difference and dispersion between conditions)
x_eff <- aldex.effect(
  x_clr,
  verbose = TRUE
)

# Combine t-test and effect size results into one table
res_aldex <- data.frame(
  taxon = rownames(x_tt),  # taxon / genus name
  x_tt,
  x_eff,
  check.names = FALSE
)

# Add convenience columns for p/q-values and significance/direction
res_aldex <- res_aldex %>%
  mutate(
    p_we   = we.ep,   # Welch p-value
    q_we   = we.eBH,  # Welch BH-FDR
    p_wi   = wi.ep,   # Wilcoxon p-value
    q_wi   = wi.eBH,  # Wilcoxon BH-FDR
    # Use the minimum FDR across the two tests as a conservative combined measure
    p_min  = pmin(q_we, q_wi, na.rm = TRUE),
    significant = !is.na(p_min) & p_min < alpha,
    direction   = case_when(
      !significant ~ "ns",
      diff.btw > 0 ~ "Tumor_enriched",  # higher in Tumor
      diff.btw < 0 ~ "NAT_enriched",    # higher in NAT
      TRUE         ~ "ns"
    )
  ) %>%
  arrange(p_min)

# Summary counts
n_sig_aldex  <- sum(res_aldex$significant)
n_up_aldex   <- sum(res_aldex$direction == "Tumor_enriched")
n_down_aldex <- sum(res_aldex$direction == "NAT_enriched")

cat("\nALDEx2 (prev 5%, paired)\n")
cat("  Significant (min FDR <", alpha, "):", n_sig_aldex, "\n")
cat("    Tumor_enriched:", n_up_aldex, "\n")
cat("    NAT_enriched  :", n_down_aldex, "\n")

# Save: all results, significant only, and a small summary table
write_csv(
  res_aldex,
  file.path(out_dir_aldex, "ALDEx2_prev5_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_aldex %>% filter(significant),
  file.path(out_dir_aldex, "ALDEx2_prev5_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_aldex,
    n_Tumor_enriched = n_up_aldex,
    n_NAT_enriched   = n_down_aldex
  ),
  file.path(out_dir_aldex, "ALDEx2_prev5_Tumor_vs_NAT_summary.csv")
)


## =========================================================
## 3) ANCOM-BC2 (paired via random effect, 5% prev, no internal filter)
## =========================================================

# Create output directory for ANCOM-BC2
out_dir_ancom <- file.path(base_dir, "ANCOMBC2")
dir.create(out_dir_ancom, showWarnings = FALSE, recursive = TRUE)

# Metadata for ANCOM-BC2 (rownames must be sample IDs)
meta_ancom <- meta_aligned
rownames(meta_ancom) <- meta_ancom[[sampleid]]

cat("\nANCOM-BC2: taxa =", nrow(count_prev5),
    ", samples =", ncol(count_prev5), "\n")

set.seed(123)
out_ancom <- ancombc2(
  data          = count_prev5,     # taxa x samples count matrix
  taxa_are_rows = TRUE,            # rows are taxa, columns are samples
  meta_data     = meta_ancom,      # sample metadata
  fix_formula   = groupcol,        # fixed effect: Group (NAT vs Tumor)
  rand_formula  = paste0("(1 | ", patienid, ")"),  # random effect: Patient_ID (paired)
  p_adj_method  = "BH",
  pseudo        = 0,               # no pseudo-count added by default
  pseudo_sens   = TRUE,            # sensitivity analysis for pseudo-counts (not a prevalence filter)
  prv_cut       = 0,               # no internal prevalence filtering (we already did 5% externally)
  lib_cut       = 0,               # no internal library-size filtering
  s0_perc       = 0.05,            # small constant for variance stabilization
  group         = NULL,
  struc_zero    = FALSE,           # do not force structural zeros
  neg_lb        = FALSE,           # do not use negative lower bound filtering
  alpha         = alpha,
  n_cl          = 1,               # number of CPU cores (can increase if you want parallel)
  verbose       = TRUE,
  global        = FALSE,
  pairwise      = FALSE,
  dunnet        = FALSE,
  trend         = FALSE
)

# Extract main ANCOM-BC2 results table
res_ancom <- out_ancom$res

# Tidy up and add convenience columns
# NOTE: column names like lfc_GroupTumor, q_GroupTumor come from Group factor (NAT vs Tumor)
res_ancom_tbl <- res_ancom %>%
  as.data.frame() %>%
  mutate(
    lfc         = lfc_GroupTumor,         # log2-fold-change Tumor vs NAT
    qval        = q_GroupTumor,           # BH-adjusted q-value
    diff        = diff_GroupTumor,        # difference (presence/absence) flag
    diff_robust = diff_robust_GroupTumor  # robust difference flag from sensitivity analysis
  ) %>%
  mutate(
    significant        = !is.na(qval) & qval < alpha,  # standard FDR-based significance
    significant_robust = diff_robust,                  # robust significance from ANCOM-BC2
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

# Summary counts
n_sig_ancom         <- sum(res_ancom_tbl$significant)
n_sig_robust_ancom  <- sum(res_ancom_tbl$significant_robust)
n_up_ancom          <- sum(res_ancom_tbl$direction == "Tumor_enriched")
n_down_ancom        <- sum(res_ancom_tbl$direction == "NAT_enriched")
n_up_robust_ancom   <- sum(res_ancom_tbl$direction_robust == "Tumor_enriched")
n_down_robust_ancom <- sum(res_ancom_tbl$direction_robust == "NAT_enriched")

cat("\nANCOM-BC2 (prev 5%, paired)\n")
cat("  Significant (q <", alpha, "):", n_sig_ancom, "\n")
cat("    Tumor_enriched:", n_up_ancom, "\n")
cat("    NAT_enriched  :", n_down_ancom, "\n")
cat("  Robust significant:", n_sig_robust_ancom, "\n")
cat("    Tumor_enriched (robust):", n_up_robust_ancom, "\n")
cat("    NAT_enriched   (robust):", n_down_robust_ancom, "\n")

# Save: full table, standard significant, and robust significant
write_csv(
  res_ancom_tbl,
  file.path(out_dir_ancom, "ANCOMBC2_prev5_Tumor_vs_NAT_full.csv")
)

write_csv(
  res_ancom_tbl %>% filter(significant),
  file.path(out_dir_ancom, "ANCOMBC2_prev5_Tumor_vs_NAT_significant_q05.csv")
)

write_csv(
  res_ancom_tbl %>% filter(significant_robust),
  file.path(out_dir_ancom, "ANCOMBC2_prev5_Tumor_vs_NAT_significant_robust_q05.csv")
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
  file.path(out_dir_ancom, "ANCOMBC2_prev5_Tumor_vs_NAT_summary.csv")
)


## =========================================================
## 4) MaAsLin2 (paired via random effect, 5% prev, no internal filter)
## =========================================================

# Create output directory for MaAsLin2
out_dir_maaslin <- file.path(base_dir, "Maaslin2")
dir.create(out_dir_maaslin, showWarnings = FALSE, recursive = TRUE)

# Minimal metadata for MaAsLin2 (Sample_ID, Group, Patient_ID)
meta_maaslin <- meta_aligned[, c(sampleid, groupcol, patienid)]
rownames(meta_maaslin) <- meta_maaslin[[sampleid]]

# Ensure factors and reference level for Group
meta_maaslin[[groupcol]]   <- factor(meta_maaslin[[groupcol]], levels = c("NAT", "Tumor"))
meta_maaslin[[patienid]]   <- factor(meta_maaslin[[patienid]])

# Sanity check: metadata rows and count matrix columns must align exactly
stopifnot(all(colnames(count_prev5) == rownames(meta_maaslin)))

# Run MaAsLin2
# fixed_effects = Group (NAT vs Tumor)
# random_effects = Patient_ID (paired / subject random effect)
# min_prevalence = 0 → no internal prevalence filter (we already did 5%)
fit_maaslin <- Maaslin2(
  input_data      = count_prev5,          # taxa x samples matrix
  input_metadata  = meta_maaslin,         # sample metadata
  output          = out_dir_maaslin,      # MaAsLin2 will write its own files here
  fixed_effects   = groupcol,             # Group effect (NAT vs Tumor)
  random_effects  = patienid,             # random intercept for Patient_ID
  normalization   = "TSS",                # total sum scaling
  transform       = "LOG",                # log transform after TSS
  analysis_method = "LM",                 # linear model
  correction      = "BH",
  standardize     = TRUE,
  min_abundance   = 0,                    # no internal abundance filter
  min_prevalence  = 0,                    # no internal prevalence filter
  min_variance    = 0,                    # no internal variance filter
  max_significance = 1                    # keep all features in output
)

# Read MaAsLin2 main results table
all_res <- read_tsv(
  file.path(out_dir_maaslin, "all_results.tsv"),
  show_col_types = FALSE
)

# Keep only results for Group (NAT vs Tumor)
res_group <- all_res %>%
  filter(metadata == groupcol) %>%      # effect of Group on each taxon
  mutate(
    significant = !is.na(qval) & qval < alpha,
    direction   = case_when(
      !significant ~ "ns",
      coef > 0     ~ "Tumor_enriched",   # positive coefficient → higher in Tumor
      coef < 0     ~ "NAT_enriched",     # negative coefficient → higher in NAT
      TRUE         ~ "ns"
    )
  ) %>%
  arrange(qval)

# Summary counts
n_sig_maaslin  <- sum(res_group$significant)
n_up_maaslin   <- sum(res_group$direction == "Tumor_enriched")
n_down_maaslin <- sum(res_group$direction == "NAT_enriched")

cat("\nMaAsLin2 (prev 5%, paired)\n")
cat("  Significant (q <", alpha, "):", n_sig_maaslin, "\n")
cat("    Tumor_enriched:", n_up_maaslin, "\n")
cat("    NAT_enriched  :", n_down_maaslin, "\n")

# Save: all Group results, significant only, and summary
write_csv(
  res_group,
  file.path(out_dir_maaslin, "Maaslin2_prev5_Tumor_vs_NAT_all.csv")
)

write_csv(
  res_group %>% filter(significant),
  file.path(out_dir_maaslin, "Maaslin2_prev5_Tumor_vs_NAT_significant.csv")
)

write_csv(
  tibble(
    alpha            = alpha,
    n_significant    = n_sig_maaslin,
    n_Tumor_enriched = n_up_maaslin,
    n_NAT_enriched   = n_down_maaslin
  ),
  file.path(out_dir_maaslin, "Maaslin2_prev5_Tumor_vs_NAT_summary.csv")
)

cat("\nDONE: DESeq2, ALDEx2, ANCOM-BC2, and MaAsLin2 run with 5% prevalence filter on RAW integer data and paired design.\n")



