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
   - Milestone 1: Evaluated ERBB2 RNA (AUC = 0.8317) vs. GISTIC copy number (AUC = 0.8139); paired DeLong test demonstrated univariate parity ($p = 0.1604$). Derived ESTIMATE tumor purity and 17q12 core amplicon co-expression (`Amplicon_Z`).
   - Milestone 2: Trained 5-fold cross-validated cost-sensitive Random Forest and Logistic Regression models. RF was selected over LR for superior collinearity handling and cardiotoxic false positive suppression. Optimizing probability threshold (0.50 -> 0.5880) via Youden's J achieved 96.70% Precision and 99.56% Specificity across 680 true negatives.
   - Milestone 3 / Step 7: Successfully reclassified 16 of 123 Indeterminate primary tumors as HER2-positive (13.0%; 12 Tier 1, 4 Tier 2) and generalized with 100% concordance on 7 held-out distant metastatic controls.
4. **Standalone Client Deliverables Package (`TCGA_BRCA_HER2_Deliverables/` & `TCGA_BRCA_HER2_Deliverables.zip`):**
   - Self-contained, responsive HTML and Markdown reports (`Classification_Report.html`, `QC_Report.html`) conforming to `plannotator/effective-html` specifications.
   - Complete prompt audit trail (`prompts/`), master audited Excel workbooks (`results/`), serialized Python models (`production_models.pkl`), and clean reproduction instructions (`README.txt`).
   - Raw client `data/` and large re-normalized CSV matrices omitted for a lightweight 10.20 MB release package.
