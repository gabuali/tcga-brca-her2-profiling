# Project Configuration

- The project focus is TCGA-BRCA HER2 profiling using IHC, FISH, copy number, and RNA-seq RSEM data.

## Completed Work
1. **Clinical HER2 Status Adjudication (`src/adjudicate_her2.py` -> `results/her2_clinical_status.csv`, `results/her2_adjudication_summary.xlsx`):**
   - Reconciled 1,108 clinical pathology records using ASCO/CAP 2023 priority guidelines (IHC score -> reflex FISH).
   - Final cohort counts: `Negative` = 767, `Positive` = 195, `Indeterminate` = 143, `Equivocal` = 3.
2. **RNA-Seq Diagnostic, QC, Filtering, Normalization & PCA Outlier Removal (`src/normalize_rnaseq.R` -> `reports/QC_Report.html`):**
   - Profiled `data/tcga.brca.rsem.csv` (1,212 samples x 20,502 genes); confirmed 75th percentile scaling triggering Gate C (`edgeR` TMM + `limma-voom`).
   - Prevalence filter retained 15,632 genes. Unit-variance PCA screening flagged and removed 19 outliers (|Z_MAD| >= 3), yielding clean 1,193-sample cohort (1,075 Primary Tumors, 111 Normal, 7 Metastasis).
3. **Multimodal Machine Learning & Comparative Evaluation (`src/her2_pipeline_*.py` -> `reports/Classification_Report.html`, `results/her2_classification_results.xlsx`):**
   - Milestone 1: Evaluated ERBB2 RNA (AUC = 0.8317) vs. GISTIC copy number (AUC = 0.8139); point estimates favor RNA across all metrics (higher AUC, sensitivity 65.8% vs. 60.4%, balanced accuracy 81.3% vs. 79.3%); DeLong test ($p = 0.1604$) confirms overlapping 95% CIs. Derived ESTIMATE tumor purity and 17q12 core amplicon co-expression (`Amplicon_Z`).
   - Milestone 2: Trained 5-fold cross-validated cost-sensitive Random Forest (OOF AUC = 0.8011) and Logistic Regression (OOF AUC = 0.8317) models. Documented that multimodal ML does not exceed univariate RNA in global AUC due to 17q12 collinearity ($r = 0.893$). Calibrated threshold shifting ($0.50 \rightarrow 0.5880$) provides a dedicated cardiotoxic rule-in tool (96.70% Precision, 99.56% Specificity, only 3 FP across 680 TN) trading off sensitivity (59.06%).
   - Milestone 3 (Unsupervised Subsetting, `src/her2_unsupervised_subsets.py` -> Plots 9 & 10): Profiled 50 PAM50 markers across 1,075 primary tumors, revealing that clinical HER2+ disease is biologically bifurcated into HER2-Enriched ($58.4\%$, pure 17q12 amplicon driver) and Luminal B / Luminal-HER2 ($27.1\%$, ER+/PR+ co-expression with dual pathway addiction), mechanistically explaining single-gene false-negative rates.
   - Milestone 4/5: Prioritized 16 of 123 Indeterminate primary tumors as HER2-positive reflex candidates (13.0%; 12 Tier 1, 4 Tier 2) with 100% concordance on 7 held-out distant metastatic controls.
4. **Standalone Client Deliverables Package & Public GitHub Repository:**
   - GitHub Repository: `https://github.com/gabuali/tcga-brca-her2-profiling` (MIT License, live clickable HTML reports, S3 bulk download instructions, 67 MB compressed RSEM matrix).
   - Standalone deliverables folder `TCGA_BRCA_HER2_Deliverables/` and zip archive `TCGA_BRCA_HER2_Deliverables.zip` (12.85 MB) with 11-sheet master Excel workbook and all 10 analytical plots.
