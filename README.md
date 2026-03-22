# Differential Abundance Analysis of the Human Tumor Microbiome

<img width="1585" height="499" alt="Screenshot 2026-03-22 174510" src="https://github.com/user-attachments/assets/f24df927-3d1e-4104-8f09-1f08b775e7db" />

# Differential Abundance Analysis of the Human Tumor Microbiome

Benchmarking of four DA tools (DESeq2, ANCOM-BC2, ALDEx2, MaAsLin2) applied to 
SMURF-normalized 16S abundance data from Nejman et al. (*Science*, 2020), comparing 
normal adjacent tissue (NAT) vs. tumor across three cancer types using paired patient 
samples. Data from Table S2 was rounded to integers and analyzed at the genus level 
across three data conditions (raw, mbImpute, BMDD) with 5% and 10% prevalence filters.

---

![Analysis Workflow](workflow.png)

*Figure 1. Analysis pipeline: three cancer types × three data conditions × two 
prevalence filters × four DA tools. Significance threshold: q < 0.05 
(Benjamini-Hochberg). Paired patient design accounted for in all tools.*

---

## Key findings

Tool sensitivity varied dramatically across conditions, consistent with Nearing et al. 
(*Nature Communications*, 2022). DESeq2 identified the most genera (up to 71 in lung, 
raw 5% filter), while ALDEx2 and ANCOM-BC2 were most conservative. Discordance was 
amplified by the pre-normalized input data, which violates the raw count assumptions 
of count-based tools. ANCOM-BC2 findings were significant under the standard test but 
not the robust test across all cancer types, reflecting the low-biomass nature of tumor 
microbiome data.

The table below shows genera found significant by more than one tool across any condition.

| Genus | Cancer | Direction | Tools |
|---|---|---|---|
| *Staphylococcus* | Breast | Tumor-enriched | DESeq2, ALDEx2, MaAsLin2 |
| *Prevotella* | Breast | Tumor-enriched | ALDEx2, MaAsLin2 |
| *Afipia* | Breast | NAT-enriched | DESeq2, MaAsLin2 |
| *Kocuria* | Breast | NAT-enriched | DESeq2, MaAsLin2 |
| *Sphaerobacter* | Breast | NAT-enriched | DESeq2, MaAsLin2 |
| *Bacillus* | Lung | NAT-enriched | DESeq2, MaAsLin2 |
| *Devosia* | Lung | NAT-enriched | DESeq2, MaAsLin2 |
| *Friedmanniella* | Lung | Tumor-enriched | ALDEx2, MaAsLin2 |
| *Staphylococcus* | Ovary | Tumor-enriched | DESeq2, ANCOM-BC2 |
| *Escherichia/Shigella* | Ovary | Tumor-enriched | DESeq2, ANCOM-BC2 |

*Staphylococcus* is the most robust finding — tumor-enriched in both breast and ovary 
across the most conservative tools, corroborating prevalence findings in Nejman et al. 
and published literature (Urbaniak et al. 2016; Hieken et al. 2016). Full per-tool, 
per-condition results are in [`Summary.txt`](Summary.txt) and 
[`Detailed_Interpretation_Report.txt`](Detailed_Interpretation_Report.txt).

---

## Repository structure
```
Differential-Abundance-Analysis/
├── Breast_new/
│   ├── Raw_data/          # R scripts + outputs for raw integer data
│   └── Imputed_data/      # R scripts + outputs for mbImpute and BMDD
├── Lung_new/              # Same structure as Breast_new
├── Ovary new/             # Same structure (no BMDD — sample size too small)
├── Supplementary_Tables/  # Original supplementary tables from Nejman et al.
├── abundance_tumor.csv    # Integer-rounded genus-level abundance matrix
├── clinical_tumor.csv     # Paired clinical metadata (NAT vs Tumor)
├── Summary.txt            # Significant genera per tool per condition
└── Detailed_Interpretation_Report.txt  # Full biological interpretation
```
