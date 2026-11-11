# TCGA-BRCA HER2 Profiling

### Multimodal Machine Learning, Clinical Adjudication, Tumor Purity Deconvolution & Phenotypic Reclassification

[![Cohort Size](https://img.shields.io/badge/Cohort_Size-1%2C193_Biospecimens-blue.svg)](#1-master-cohort-flow--pathology-reconciliation)
[![Model](https://img.shields.io/badge/Production_Model-Cost--Sensitive_Random_Forest-emerald.svg)](#4-machine-learning-architecture--probability-shifting)
[![Precision](https://img.shields.io/badge/Precision-96.70%25-brightgreen.svg)](#4-machine-learning-architecture--probability-shifting)
[![Specificity](https://img.shields.io/badge/Specificity-99.56%25-brightgreen.svg)](#4-machine-learning-architecture--probability-shifting)
[![License](https://img.shields.io/badge/License-MIT-purple.svg)](LICENSE)

---

## Executive Overview & Translational Context

Accurate determination of Human Epidermal Growth Factor Receptor 2 (**HER2 / ERBB2**) status is critical for guiding targeted monoclonal antibody (trastuzumab, pertuzumab) and antibody-drug conjugate (trastuzumab deruxtecan) therapies in breast invasive carcinoma. Clinical standard-of-care relies on **Immunohistochemistry (IHC)** with reflex **Fluorescence In Situ Hybridization (FISH)** for equivocal (IHC 2+) cases. However, clinical pathology workflows frequently encounter non-definitive assessments, technical discordance, and **tumor microenvironment (TME) stromal dilution** that obscure true oncogenic driver amplification.

This repository provides an audited, end-to-end translational bioinformatics and machine learning pipeline developed on **The Cancer Genome Atlas Breast Invasive Carcinoma (TCGA-BRCA)** cohort. By integrating diagnostic IHC/FISH pathology records, GISTIC 2.0 copy-number calls, and Level 3 RNA-Seq expression profiles, the framework:
1. **Adjudicates Clinical Discordance:** Harmonizes 1,108 pathology records under ASCO/CAP 2023 priority guidelines, resolving equivocal IHC 2+ cases via reflex FISH and flagging diagnostic discrepancies.
2. **Standardizes High-Throughput RNA-Seq:** Programmatically profiles raw RSEM counts, confirms uniform 75th percentile scaling (~1,000 counts, CV = 0.0456%), and routes data through `edgeR` TMM library scaling and `limma-voom` variance stabilization across 15,632 filtered genes.
3. **Screens Multidimensional PCA Outliers:** Employs unit-variance PCA (`scale. = TRUE`) to prevent single-gene hyper-variable lineage drivers (*ESR1*, $r = 0.80$ with unscaled PC2) from falsely discarding 145 ER-negative/Basal tumors, isolating 19 technical outliers ($|Z_{\text{MAD}}| \ge 3.0$).
4. **Determines RNA Predictive Superiority:** Continuous *ERBB2* RNA-seq outperforms discrete GISTIC copy number across every empirical metric ($\text{AUC } 0.8317 \text{ vs. } 0.8139$, Sensitivity $65.8\% \text{ vs. } 60.4\%$, Balanced Accuracy $81.3\% \text{ vs. } 79.3\%$). Paired DeLong's test ($p = 0.1604$) confirms overlapping 95% confidence intervals, but continuous RNA provides superior resolution for borderline and low-gain (+1) cases.
5. **Deconvolutes Subtypes via Unsupervised Subsetting:** Profiles 50 PAM50 intrinsic markers across $1,075$ primary tumors, revealing that clinical HER2+ disease is biologically bifurcated: **$58.4\%$ are HER2-Enriched** (pure 17q12 amplicon addiction, low ER/PR, high proliferation), while **$27.1\%$ are Luminal B / Luminal-HER2** (co-expressing high estrogen and progesterone receptors). This biological dichotomy directly explains why single-gene univariate cutoffs miss ~34% of clinical positives.
6. **Appraises Multimodal ML & Calibrates Cardiotoxic Safety:** Demonstrates transparently that multimodal ML (RF OOF AUC $0.8011$, LR OOF AUC $0.8317$) does not expand global AUC beyond univariate RNA ($0.8317$) due to heavy 17q12 collinearity ($r = 0.893$). However, shifting the decision threshold ($0.50 \rightarrow 0.5880$) configures the Random Forest as a specialized high-precision rule-in classifier (**$96.70\%$ precision**, **$99.56\%$ specificity**, only 3 false positives across 680 true negatives) to prevent cardiotoxic overtreatment.
7. **Prioritizes Indeterminate Patients:** Identifies **16 of 123 Indeterminate primary tumors as high-priority HER2-positive candidates** ($13.0\%$) exhibiting multi-omic concordance and enrichment in the HER2-Enriched and Luminal-HER2 clusters for targeted reflex testing, while demonstrating **100% concordance** across 7 held-out distant metastatic controls.

---

## Interactive Technical Reports (Clickable Deliverables)

The technical deliverables have been compiled as standalone, responsive HTML documents adhering to publication specifications (featuring sticky sidebar navigation, MathJax formula cards, image lightbox modal zoom, and horizontal table scrolling) alongside GitHub-flavored Markdown versions:

| Report | Description | Direct Repo Link | 1-Click Interactive Web Preview |
| :--- | :--- | :---: | :---: |
| **Multimodal HER2 Status Classification** | Comparative evaluation of RNA vs. DNA predictive potential, ESTIMATE tumor purity, 17q12 amplicon co-expression, cost-sensitive Random Forest vs. Logistic Regression, Youden threshold shifting, and multi-omic reclassification of 123 Indeterminate tumors. | [`Classification_Report.html`](reports/Classification_Report.html) &bull; [`Markdown`](reports/Classification_Report.md) | [**Launch Interactive Report**](https://htmlpreview.github.io/?https://github.com/gabuali/tcga-brca-her2-profiling/blob/main/reports/Classification_Report.html) |
| **RNA-Seq Diagnostic QC & PCA Screening** | Blind profiling logic gates (Gate C), upper-quartile RSEM validation, `edgeR` TMM library scaling, `limma-voom` variance stabilization, prevalence filtering, and unit-variance PCA outlier screening. | [`QC_Report.html`](reports/QC_Report.html) &bull; [`Markdown`](reports/QC_Report.md) | [**Launch Interactive Report**](https://htmlpreview.github.io/?https://github.com/gabuali/tcga-brca-her2-profiling/blob/main/reports/QC_Report.html) |

---

## Data Sourcing & Direct S3 Bulk Package Download

The primary datasets for this pipeline originate from the **TCGA Breast Invasive Carcinoma (TCGA-BRCA)** project (cBioPortal study identifier: `brca_tcga`).

### 1. Direct Bulk Package Download (All 3 Files in One Archive)
To download the complete TCGA-BRCA study archive containing the clinical records (`data_clinical_patient.txt`, `data_clinical_sample.txt`), discrete copy-number calls (`data_cna.txt`), and the RSEM expression matrix (`data_mrna_seq_v2_rsem.txt`) directly from the official cBioPortal S3 bucket:

```bash
# Direct S3 archive for the complete brca_tcga study
wget https://cbioportal-datahub.s3.amazonaws.com/brca_tcga.tar.gz

# Extract contents
tar -xzf brca_tcga.tar.gz
```

### 2. Repository Data Layout
The `data/` directory in this repository contains pre-formatted versions of these files ready for pipeline execution:
* `data/brca_tcga_clinical_data.csv` (1.08 MB): 1,108 clinical pathology records with 140 diagnostic fields.
* `data/brca_tcga_erbb2_copy_number.csv` (30 KB): GISTIC 2.0 discrete copy-number calls for *ERBB2*.
* `data/tcga.brca.rsem.csv.gz` (66.97 MB): Level 3 RSEM Upper-Quartile normalized expression matrix (1,212 samples $\times$ 20,502 genes). 
  * *Note on Compression:* To comply with GitHub's strict 100 MB per-file limit, the raw 182.90 MB CSV matrix is stored as a standard `.gz` archive. `src/normalize_rnaseq.R` automatically detects and streams `tcga.brca.rsem.csv.gz` via `data.table::fread` without requiring manual decompression.

### 3. Note on Omitted Large Re-Normalized Matrices
To keep git clones lightweight and responsive, the two re-normalized Level 3 matrices (~632 MB combined):
* `results/matrix_normalized_linear.csv` (313.6 MB)
* `results/matrix_normalized_log2_stabilized.csv` (319.0 MB)

are excluded from git tracking. All machine learning models, cross-validation runs, and reports run standalone using the serialized binary bundles enclosed in `results/`. Both CSV matrices can be regenerated on-demand in ~45 seconds by running `Rscript src/normalize_rnaseq.R`.

---

## Master Cohort Flow & Pathology Reconciliation

```
[Master Clinical Pathology Records] (N = 1,108)
   ├── Negative: 767
   ├── Positive: 195
   ├── Indeterminate: 143
   └── Equivocal: 3
         │
         ▼
[Level 3 RNA-Seq RSEM Biospecimens] (N = 1,212)
   ├── Primary Solid Tumors: 1,093
   ├── Adjacent Normal Tissue: 112
   └── Distant Metastatic: 7
         │
         ▼ (Unit-Variance PCA Screening |Z_MAD| >= 3.0)
[Flagged Technical Outliers Removed] (N = 19)
   ├── Primary Tumors Flagged: 18
   └── Normal Adjacent Flagged: 1
         │
         ▼
[Post-QC Clean Genomic Cohort] (N = 1,193)
   ├── Clean Primary Tumors: 1,075  <--- Evaluated in ML Pipeline
   ├── Clean Adjacent Normal: 111   <--- Normal Baseline Reference
   └── Distant Metastatic Controls: 7 <--- Held-Out Validation Controls
```

### Complete Accounting of the 33-Sample Primary Exclusion Delta
* Master Clinical Primary Tumors: **1,108**
* Clean Post-QC Primary Tumors with RNA-seq: **1,075**
* **Total Delta: 33 Samples**
  1. **18 Technical PCA Outliers:** Flagged and removed during unit-variance PCA screening (`results/pca_outliers_removed.csv`).
  2. **7 Distant Metastatic Biospecimens:** Preserved as held-out controls for biological generalizability testing.
  3. **8 Clinical Records Lacking Genomic Profiles:** Clinical records lacking corresponding high-throughput RNA-seq profiles in the primary sequencing release.

---

## Machine Learning Architecture & Probability Shifting

```
Multimodal Input Vector:
  ├── ERBB2_log2 (Normalized RNA Expression) [50.89% RF Decision Weight]
  ├── Amplicon_Z (17q12 Core Co-expression: GRB7, PGAP3, STARD3, MIEN1) [31.95% Weight]
  └── ESTIMATE_Purity (Tumor Microenvironment Deconvolution) [17.16% Weight]
         │
         ▼
5-Fold Stratified Cross-Validation (Cost-Sensitive 4.56x Class Balancing)
         │
         ├── Logistic Regression (LR): Precision = 89.52%, Specificity = 98.38% (11 False Positives)
         └── Random Forest (RF):       Precision = 94.00%, Specificity = 99.12% (6 False Positives)
                 │
                 ▼ (Youden's J Threshold Shift: 0.5000 -> 0.5880)
         Calibrated High-Precision Rule-In Model:
                 ├── Precision:   96.70% (+3.08% improvement; only 3 FP across 680 TN)
                 ├── Specificity: 99.56% (Halved false positives: 6 -> 3)
                 ├── Sensitivity: 59.06% (Deliberate trade-off to minimize cardiotoxic risk)
                 └── Cardiotoxic False Positive Prevention: Maximized
```

---

## Repository Directory Architecture

```
.
├── README.md                                 # Primary repository guide & documentation
├── README.txt                                # Plain-text technical release guide
├── .gitignore                                # Excludes raw uncompressed RSEM and large matrices
├── data/                                     # Clinical, copy-number and genomic expression data
│   ├── brca_tcga_clinical_data.csv           # Diagnostic clinical pathology records (N = 1,108)
│   ├── brca_tcga_erbb2_copy_number.csv       # GISTIC 2.0 discrete copy-number calls
│   └── tcga.brca.rsem.csv.gz                 # Compressed Level 3 RSEM matrix (66.97 MB)
├── reports/                                  # Publication-grade interactive deliverables
│   ├── Classification_Report.html            # Standalone responsive HTML classification report (with Section 8 Subtypes)
│   ├── Classification_Report.md              # Standard GitHub-flavored Markdown report
│   ├── QC_Report.html                        # Standalone responsive HTML diagnostic QC report
│   └── QC_Report.md                          # Standard GitHub-flavored Markdown report
├── prompts/                                  # Complete AI governance, audit trail & execution prompts
│   ├── 01_her2_status_adjudication_prompt.txt   # Phase 1: Clinical adjudication prompt
│   ├── 02_qc_agent_prompt.txt                   # Phase 2: RNA-seq diagnostic gates & PCA QC prompt
│   ├── 03_her2_classification_agent_prompt.txt  # Phase 3: Multimodal ML classification prompt
│   ├── 04a_user_reports_request.txt             # Phase 4 User Directive: Strategic reporting prompt
│   └── 04b_client_reports_generation_prompt.txt # Phase 4 Technical Prompt: Technical specification
├── src/                                      # Modular pipeline source code
│   ├── adjudicate_her2.py                    # Phase 1: Clinical adjudication & discordance auditing
│   ├── normalize_rnaseq.R                    # Phase 2: Diagnostic profiling, TMM-voom & PCA screening
│   ├── her2_pipeline_milestone1.py           # Phase 3A: Milestone 1 (Plots 1-6, purity, DeLong)
│   ├── her2_pipeline_milestone2.py           # Phase 3B: Milestone 2 (5-Fold CV, RF/LR threshold shifting)
│   ├── her2_unsupervised_subsets.py          # Phase 3C: Milestone 3 (PAM50 intrinsic clustering, Plots 9 & 10)
│   ├── her2_pipeline_step7.py                # Phase 3D: Milestone 4/5 (Prioritization, Plot 8, metastasis)
│   └── compile_client_reports.py             # Phase 4: Standalone HTML/MD technical report compiler
└── results/                                  # Quantitative outputs, models, workbooks & figures
    ├── her2_classification_results.xlsx      # Master 11-sheet audited Excel classification workbook
    ├── her2_adjudication_summary.xlsx        # 3-sheet clinical contingency and discordance workbook
    ├── her2_clinical_status.csv              # Adjudicated clinical status table (1,108 samples)
    ├── brca_her2_status_adjudicated.tsv      # Tab-delimited status table with discordance flags
    ├── rsem_sample_metadata.csv              # Complete sample QC metadata and library sizes
    ├── pca_outliers_removed.csv              # Clinical and genomic metadata for 19 flagged outliers
    ├── milestone1_bundle.pkl                 # Serialized Step 1-5 data structures and ROC metrics
    ├── production_models.pkl                 # Serialized production RF and LR models (6.37 MB)
    ├── unsupervised_bundle.pkl               # Serialized PAM50 subtype clusters, PCA & contingency
    ├── step7_bundle.pkl                      # Serialized Step 7 prioritization tiers & metastasis
    ├── pca_before_outlier_removal.png        # Figure 1A: PCA before outlier removal
    ├── pca_after_outlier_removal.png         # Figure 1B: PCA clean cohort distribution
    ├── plot1_complete_cases_kde.png          # Figure 1: ERBB2 expression KDE across complete cases
    ├── plot2_rnaseq_roc.png                  # Figure 2: Univariate ERBB2 RNA-seq ROC curve
    ├── plot3_purity_vs_erbb2.png             # Figure 3: ESTIMATE tumor purity vs. ERBB2 expression
    ├── plot4_amplicon_clustermap.png         # Figure 4: 17q12 core amplicon co-expression clustermap
    ├── plot5_gistic_roc.png                  # Figure 5: Discrete GISTIC copy number ROC curve
    ├── plot6_combined_roc.png                # Figure 6: Paired RNA vs. GISTIC ROC overlay with DeLong
    ├── plot7_rf_feature_importance.png       # Figure 7: 5-fold cross-validated RF feature importance
    ├── plot8_rescue_facet_plot.png           # Figure 8: 3-panel Indeterminate prioritization facet plot
    ├── plot9_unsupervised_pca_subtypes.png   # Figure 9: PAM50 unsupervised PCA subtype segregation biplot
    └── plot10_subtype_marker_heatmap.png     # Figure 10: Clustered heatmap of canonical lineage markers
```

---

## Reproduction & Pipeline Execution

To execute the complete pipeline from data processing to client report compilation, run the following sequential commands:

```bash
# Step 1: Clinical Adjudication & Discordance Auditing
python3 src/adjudicate_her2.py

# Step 2: Diagnostic RNA-Seq QC, TMM-Voom Normalization & PCA Screening
# (Reads data/tcga.brca.rsem.csv.gz directly; completes in ~45 seconds)
Rscript src/normalize_rnaseq.R

# Step 3: Milestone 1 Analysis (Plots 1-6, Purity, Amplicon_Z, DeLong Test)
python3 src/her2_pipeline_milestone1.py

# Step 4: Milestone 2 Analysis (5-Fold Stratified CV, RF/LR, Plot 7)
python3 src/her2_pipeline_milestone2.py

# Step 5: Milestone 3 Analysis (Unsupervised PAM50 Subsetting, Plots 9 & 10)
python3 src/her2_unsupervised_subsets.py

# Step 6: Step 7 Analysis (Indeterminate Prioritization, Plot 8, Metastasis)
python3 src/her2_pipeline_step7.py

# Step 7: Compile Standalone Client Technical Reports (HTML & Markdown)
python3 src/compile_client_reports.py
```

---

## Software Dependencies & Environment

* **Python (>= 3.9):** `numpy`, `pandas`, `scipy`, `scikit-learn`, `matplotlib`, `seaborn`, `openpyxl`, `jinja2`, `markdown`.
* **R (>= 4.1):** `data.table`, `edgeR`, `limma`, `ggplot2`.

---

## License

This project is licensed under the [MIT License](LICENSE) - see the LICENSE file for details.
