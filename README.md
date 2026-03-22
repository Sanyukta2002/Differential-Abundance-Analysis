# Differential Abundance Analysis of the Human Tumor Microbiome

Benchmarking of four DA tools applied to SMURF-normalized 16S abundance data 
from Nejman et al. (*Science*, 2020), comparing tumor vs. normal adjacent tissue 
(NAT) across three cancer types using paired patient samples.

---

![Analysis Workflow](<img width="1585" height="499" alt="Screenshot 2026-03-22 174510" src="https://github.com/user-attachments/assets/f24df927-3d1e-4104-8f09-1f08b775e7db" />)

*Figure 1. Analysis workflow. Three cancer types (breast n=308, lung n=456, ovary 
n=24) were analyzed across three data conditions (raw integers, mbImpute, BMDD) 
with two prevalence filters (5%, 10%), using four DA tools with a paired patient 
design (q < 0.05, Benjamini-Hochberg correction).*

---

## Data

| Property | Detail |
|---|---|
| Source | Nejman et al. "The human tumor microbiome is composed of tumor type-specific intracellular bacteria." *Science* 2020 |
| Input table | Table S2 — SMURF-normalized 16S abundance data, genus level |
| Comparison | Normal adjacent tissue (NAT) vs. tumor, paired patients only |
| Samples | Breast: 308 (154 NAT + 154 T) · Lung: 456 (228 NAT + 228 T) · Ovary: 24 (12 NAT + 12 T) |
| Note | Raw sequencing data was not publicly available. SMURF-normalized float values were rounded to integers to satisfy count-based tool requirements. |

---

## Methods

### DA tools

| Tool | Model | Normalization | Paired design |
|---|---|---|---|
| DESeq2 | Negative binomial | poscounts size factors | Patient_ID in design formula |
| ALDEx2 | CLR + Monte Carlo Dirichlet | Compositional (internal) | `paired.test = TRUE` |
| ANCOM-BC2 | Linear model + bias correction | Sampling fraction correction | Patient_ID as random effect |
| MaAsLin2 | Linear model | TSS + log transform | Patient_ID as random effect |

### Analysis conditions

| Condition | Description | Cancers |
|---|---|---|
| Raw integers | No imputation | Breast, Lung, Ovary |
| mbImpute | Regression-based zero imputation | Breast, Lung, Ovary |
| BMDD | Bimodal Dirichlet Distribution imputation | Breast, Lung only (n too small for Ovary) |

Each condition run with 5% and 10% prevalence filters. Significance threshold: q < 0.05 (Benjamini-Hochberg).

---

## Key findings

- **Staphylococcus** is the most robust finding — consistently tumor-enriched in breast across DESeq2, ALDEx2, and MaAsLin2 under both raw and imputed conditions, and in ovary via DESeq2 and ANCOM-BC2. This cross-tool consistency on pre-normalized data corroborates Nejman et al.'s prevalence findings and published literature (Urbaniak et al. 2016; Hieken et al. 2016).
- **Tool sensitivity varies dramatically**: DESeq2 identified 51 significant genera in breast (raw, 5% filter) while ALDEx2 identified 1, consistent with the benchmarking findings of Nearing et al. (*Nature Communications*, 2022). This discordance is likely amplified by the pre-normalized input data, which violates the raw count assumptions of count-based tools (DESeq2, MaAsLin2).
- **ANCOM-BC2 results were significant under the standard test but not the robust test** (`diff_robust = FALSE`) across all cancer types, indicating that findings are sensitive to the reference taxa assumption — expected in low-biomass tumor microbiome settings where many taxa may be simultaneously altered.
- **Ovary (n=24) was severely underpowered**: ALDEx2 and MaAsLin2 found no significant genera under any condition.

---

## Results

Full per-tool, per-condition results are in:
- [`Summary.txt`](Summary.txt) — significant genera per tool per condition for all three cancers
- [`Detailed_Interpretation_Report.txt`](Detailed_Interpretation_Report.txt) — full biological interpretation, tool-by-tool discussion, and limitations

---

## Repository structure
```
Differential-Abundance-Analysis/
├── Breast_new/
│   ├── Raw_data/          # R scripts + outputs for raw integer data
│   └── Imputed_data/      # R scripts + outputs for mbImpute and BMDD
├── Lung_new/              # Same structure as Breast_new
├── Ovary new/             # Same structure (no BMDD)
├── Supplementary_Tables/  # Original supplementary tables from Nejman et al.
├── abundance_tumor.csv    # Integer-rounded genus-level abundance matrix
├── clinical_tumor.csv     # Paired clinical metadata (NAT vs Tumor)
├── Summary.txt            # Significant genera per condition
└── Detailed_Interpretation_Report.txt
```

---

## References

1. Nejman D, et al. The human tumor microbiome is composed of tumor type-specific intracellular bacteria. *Science*. 2020;368(6494):973–980.
2. Nearing JT, et al. Microbiome differential abundance methods produce different results across 38 datasets. *Nature Communications*. 2022;13:342.
3. Love MI, Huber W, Anders S. DESeq2. *Genome Biology*. 2014;15:550.
4. Fernandes AD, et al. ALDEx2. *Microbiome*. 2014;2:15.
5. Lin H, Peddada SD. ANCOM-BC2. *Nature Communications*. 2020;11:3514.
6. Mallick H, et al. MaAsLin2. *PLoS Computational Biology*. 2021;17(11):e1009442.
7. Jiang R, Li WV, Li JJ. mbImpute. *Genome Biology*. 2021;22:192.
8. Zhou H, Chen J, Zhang X. BMDD. 2025. https://github.com/zhouhj1994/BMDD

<img width="1585" height="499" alt="Screenshot 2026-03-22 174510" src="https://github.com/user-attachments/assets/f24df927-3d1e-4104-8f09-1f08b775e7db" />


