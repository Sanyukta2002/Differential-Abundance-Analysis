## ------------------------------------------------------------
## 0. Setup
## ------------------------------------------------------------

setwd("C:/Users/sanyu/OneDrive/Desktop/Lung_new/raw/bmdd_lung")

suppressPackageStartupMessages({
  library(BMDD)
  library(tidyverse)
})

## ------------------------------------------------------------
## 1. Load data
## ------------------------------------------------------------

abundance_raw <- read.csv(
  "genus_abundance_lung_paired.csv",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

clinical_raw <- read.csv(
  "clinical.csv",
  check.names      = FALSE,
  stringsAsFactors = FALSE
)

cat("Abundance table: ", nrow(abundance_raw), " taxa x ",
    ncol(abundance_raw), " columns\n")
cat("Clinical table:  ", nrow(clinical_raw), " samples x ",
    ncol(clinical_raw), " columns\n")

## ------------------------------------------------------------
## 2. Split taxonomy vs sample columns
##    (names must match your lung abundance file)
## ------------------------------------------------------------

tax_cols <- c(
  "genus","Prevalence in DNA extraction/NTC negative controls",
  "Prevalence in Paraf. Controls", 
  "domain", "phylum", "class", "order", "family"
)



stopifnot(all(tax_cols %in% colnames(abundance_raw)))

tax_info   <- abundance_raw[, tax_cols]
sample_cols <- setdiff(colnames(abundance_raw), tax_cols)

otu_mat <- abundance_raw[, sample_cols]
rownames(otu_mat) <- abundance_raw$genus

cat("OTU matrix: ", nrow(otu_mat), " taxa x ",
    ncol(otu_mat), " samples (raw)\n")

## ------------------------------------------------------------
## 3. Align abundance with clinical metadata
## ------------------------------------------------------------

stopifnot("Sample_ID" %in% colnames(clinical_raw))

common_ids <- intersect(colnames(otu_mat), clinical_raw$Sample_ID)

cat("Common sample IDs between abundance and clinical: ",
    length(common_ids), "\n")

otu_mat  <- otu_mat[, common_ids, drop = FALSE]
clinical <- clinical_raw[match(common_ids, clinical_raw$Sample_ID), ]

stopifnot(all(colnames(otu_mat) == clinical$Sample_ID))

otu_mat <- as.matrix(otu_mat)
mode(otu_mat) <- "numeric"

## ------------------------------------------------------------
## 4. Taxon prevalence filtering ONLY (10%)
## ------------------------------------------------------------

prev <- 0.10  # 10%

keep_taxa <- apply(otu_mat, 1, function(x) {
  sum(x > 0) >= (ncol(otu_mat) * prev)
})

feature.dat <- otu_mat[keep_taxa, , drop = FALSE]

cat("After 10% taxon prevalence filtering: ",
    nrow(feature.dat), " taxa x ",
    ncol(feature.dat), " samples\n")

## ------------------------------------------------------------
## 5. Run BMDD
## ------------------------------------------------------------

m <- nrow(feature.dat)
n <- ncol(feature.dat)

cat("Running BMDD on ", m, " taxa x ", n, " samples ...\n")

bmdd.obj <- bmdd(
  W     = feature.dat,
  type  = "count",
  trace = TRUE
)

## ------------------------------------------------------------
## 6. Posterior mean composition
## ------------------------------------------------------------

beta <- bmdd.obj$beta

post.mean <- t(t(beta) / colSums(beta))

cat("Posterior mean matrix: ",
    nrow(post.mean), " taxa x ",
    ncol(post.mean), " samples\n")
cat("Column sum range (should be ~1): ",
    paste(range(colSums(post.mean)), collapse = " - "), "\n")

## ------------------------------------------------------------
## 7. Save BMDD-imputed abundance
##    (name this to match what your DA script will read)
## ------------------------------------------------------------

imputed_abundance <- as.data.frame(post.mean)

tax_info_ordered <- tax_info[match(rownames(imputed_abundance),
                                   tax_info$genus), ]

bmdd_imputed_final <- cbind(tax_info_ordered, imputed_abundance)

out_file <- "genus_abundance_lung_paired_BMDD_imputed_posterior_mean.csv"

write.csv(
  bmdd_imputed_final,
  out_file,
  row.names = FALSE
)

cat("BMDD-imputed table saved to: ", out_file, "\n")
