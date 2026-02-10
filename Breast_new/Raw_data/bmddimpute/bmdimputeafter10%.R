## =========================================================
## 0) Setup
## =========================================================


setwd("C:/Users/sanyu/OneDrive/Desktop/Lung_new/raw/bmdd_lung")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ALDEx2)
  library(ANCOMBC)
  library(Maaslin2)
})

alpha    <- 0.05
sampleid <- "Sample_ID"
groupcol <- "Group"
patienid <- "Patient_ID"

base_dir <- "DA_BMDD_genus"
dir.create(base_dir, showWarnings = FALSE, recursive = TRUE)

## =========================================================
## 1) Load BMDD-imputed abundance + clinical
## =========================================================

bmdd_abund <- read.csv(
  "genus_abundance_lung_paired_BMDD_imputed_posterior_mean.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

clinical <- read.csv(
  "clinical.csv",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

taxa_cols <- c(
  "genus","family","order","class","phylum","domain",
  "Prevalence in DNA extraction/NTC negative controls",
  "Prevalence in Paraf. Controls"
)

rownames(bmdd_abund) <- bmdd_abund$genus

bmdd_prop <- bmdd_abund[, setdiff(colnames(bmdd_abund), taxa_cols), drop = FALSE]
bmdd_prop <- as.matrix(bmdd_prop)
mode(bmdd_prop) <- "numeric"

cat("BMDD matrix:", nrow(bmdd_prop), "taxa x",
    ncol(bmdd_prop), "samples\n")
cat("Column sum range:",
    paste(range(colSums(bmdd_prop)), collapse=" - "), "\n")

## =========================================================
## 2) Align with clinical
## =========================================================

common_ids <- intersect(colnames(bmdd_prop), clinical[[sampleid]])
bmdd_prop  <- bmdd_prop[, common_ids, drop = FALSE]
clinical   <- clinical[match(common_ids, clinical[[sampleid]]), ]

stopifnot(all(colnames(bmdd_prop) == clinical[[sampleid]]))

clinical[[groupcol]] <- factor(clinical[[groupcol]], levels = c("NAT","Tumor"))
clinical[[patienid]] <- factor(clinical[[patienid]])

## =========================================================
## 3) Create pseudo-counts for ALDEx2 + ANCOM-BC2
## =========================================================

lib_size     <- 1e6
bmdd_pseudo  <- round(bmdd_prop * lib_size)

cat("Pseudo-count range:",
    paste(range(bmdd_pseudo), collapse=" - "), "\n")

## =========================================================
## 4) ALDEx2 (paired)
## =========================================================

out_dir_aldex <- file.path(base_dir, "ALDEx2_BMDD")
dir.create(out_dir_aldex, showWarnings=FALSE, recursive=TRUE)

meta_ord <- clinical %>% arrange(.data[[patienid]], .data[[groupcol]])
counts_aldex <- bmdd_pseudo[, meta_ord[[sampleid]], drop=FALSE]

conds <- as.character(meta_ord[[groupcol]])
print(table(conds))

set.seed(123)
x_clr <- aldex.clr(counts_aldex, conds, mc.samples=128, denom="all")
x_tt  <- aldex.ttest(x_clr, paired.test=TRUE)
x_eff <- aldex.effect(x_clr)

res_aldex <- data.frame(taxon=rownames(x_tt), x_tt, x_eff)

res_aldex <- res_aldex %>%
  mutate(
    p_min = pmin(we.eBH, wi.eBH, na.rm = TRUE),
    significant = p_min < alpha,
    direction = case_when(
      !significant ~ "ns",
      diff.btw > 0 ~ "Tumor_enriched",
      diff.btw < 0 ~ "NAT_enriched"
    )
  ) %>% arrange(p_min)

write_csv(res_aldex,
          file.path(out_dir_aldex,"ALDEx2_BMDD_Tumor_vs_NAT_all.csv"))

write_csv(res_aldex %>% filter(significant),
          file.path(out_dir_aldex,"ALDEx2_BMDD_Tumor_vs_NAT_significant.csv"))

## =========================================================
## 5) ANCOM-BC2 (paired)
## =========================================================

out_dir_ancom <- file.path(base_dir, "ANCOMBC2_BMDD")
dir.create(out_dir_ancom, showWarnings=FALSE, recursive=TRUE)

meta_ancom <- clinical
rownames(meta_ancom) <- meta_ancom[[sampleid]]

set.seed(123)
out_ancom <- ancombc2(
  data = bmdd_pseudo,
  taxa_are_rows = TRUE,
  meta_data = meta_ancom,
  fix_formula = groupcol,
  rand_formula = paste0("(1 | ", patienid, ")"),
  prv_cut = 0,
  lib_cut = 0,
  pseudo = 0,
  pseudo_sens = TRUE,
  alpha = alpha
)

res_ancom <- out_ancom$res

res_ancom_tbl <- res_ancom %>%
  as.data.frame() %>%
  mutate(
    lfc = lfc_GroupTumor,
    qval = q_GroupTumor,
    significant = qval < alpha,
    direction = case_when(
      !significant ~ "ns",
      lfc > 0 ~ "Tumor_enriched",
      lfc < 0 ~ "NAT_enriched"
    )
  ) %>% arrange(qval)

write_csv(res_ancom_tbl,
          file.path(out_dir_ancom,"ANCOMBC2_BMDD_Tumor_vs_NAT_full.csv"))

write_csv(res_ancom_tbl %>% filter(significant),
          file.path(out_dir_ancom,"ANCOMBC2_BMDD_Tumor_vs_NAT_significant.csv"))

## =========================================================
## 6) MaAsLin2 (direct BMDD compositions)
## =========================================================

out_dir_maaslin <- file.path(base_dir, "Maaslin2_BMDD")
dir.create(out_dir_maaslin, showWarnings=FALSE, recursive=TRUE)

meta_maaslin <- clinical[, c(sampleid, groupcol, patienid)]
rownames(meta_maaslin) <- meta_maaslin[[sampleid]]

fit_maaslin <- Maaslin2(
  input_data     = bmdd_prop,
  input_metadata = meta_maaslin,
  output         = out_dir_maaslin,
  fixed_effects  = groupcol,
  random_effects = patienid,
  normalization  = "NONE",
  transform      = "LOG",
  correction     = "BH",
  min_abundance  = 0,
  min_prevalence = 0,
  min_variance   = 0
)

all_res <- read_tsv(file.path(out_dir_maaslin,"all_results.tsv"))
res_group <- all_res %>% filter(metadata == groupcol)

write_csv(res_group,
          file.path(out_dir_maaslin,"Maaslin2_BMDD_Tumor_vs_NAT_all.csv"))

write_csv(res_group %>% filter(qval < alpha),
          file.path(out_dir_maaslin,"Maaslin2_BMDD_Tumor_vs_NAT_significant.csv"))

cat("\n✅ DONE: ALDEx2, ANCOM-BC2, MaAsLin2 completed using BMDD-imputed genus data\n")
