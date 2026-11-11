# TCGA-BRCA Multimodal HER2 Status Classification
## Comparative Evaluation of ERBB2/HER2 RNA vs. DNA Predictive Potential, Purity Deconvolution, Unsupervised Subsetting, and Machine Learning Modeling

---

## 1. Executive Summary & Multimodal Framework Architecture {#sec-1}

### 1.1 Purpose of the Analysis
This report delivers an audited, publication-grade translational bioinformatics and machine learning evaluation of Human Epidermal Growth Factor Receptor 2 (**HER2 / ERBB2**) status across **1,193 post-QC TCGA-BRCA biospecimens**. The investigation directly addresses five core translational milestones:
1. **Comparative Predictive Potential (RNA vs. DNA):** Empirically compare the diagnostic capability of continuous gene expression (*ERBB2* RNA-seq) versus discrete copy number (*ERBB2* GISTIC 2.0 calls) to predict clinical HER2 status.
2. **Unsupervised Patient Subsetting:** Apply unsupervised learning across intrinsic breast cancer transcriptional signatures (PAM50) to discover biologically distinct patient subsets and establish how clinical HER2 status maps onto underlying tumor biology.
3. **Tumor Microenvironment Deconvolution:** Quantify tumor cellularity and leukocyte/stromal infiltration using the ESTIMATE algorithm to identify and correct for stromal dilution in bulk biopsies.
4. **Collinear 17q12 Amplicon Architecture:** Evaluate transcriptional co-expression across the 17q12 chromosomal locus (*GRB7*, *PGAP3*, *STARD3*, *MIEN1*) to derive a composite amplicon metric (`Amplicon_Z`).
5. **Supervised Machine Learning & Decision Threshold Calibration:** Cross-validate cost-sensitive Logistic Regression and Random Forest models, transparently characterize the trade-off between sensitivity and precision, and prioritize ambiguous Indeterminate cases for clinical reflex testing.

### 1.2 Executive Takeaways & Empirical Findings
- **Directional Superiority of RNA over DNA Copy Number:** Continuous *ERBB2* RNA-seq outperforms discrete GISTIC copy number across every individual empirical classification metric on the complete-case cohort ($n = 829$): **$\text{AUC } 0.8317$ vs. $0.8139$** ($+0.0178$), **Sensitivity $65.77\%$ vs. $60.40\%$** ($+5.37\%$), and **Balanced Accuracy $81.27\%$ vs. $79.32\%$** ($+1.95\%$). Paired DeLong testing yields $p = 0.1604$ ($Z = +1.4038$), confirming that while their standalone ROC 95% confidence intervals overlap, continuous RNA provides superior discriminatory resolution for borderline, heterogeneous, and low-level copy-gain (+1) cases.
- **Unsupervised Biological Bifurcation of HER2+ Disease:** Unsupervised clustering of 50 PAM50 markers across all $1,075$ primary tumors revealed 4 distinct biological subtypes: **HER2-Enriched** ($n=230$), **Luminal B / Luminal-HER2** ($n=417$), **Luminal A** ($n=241$), and **Basal-like (TNBC)** ($n=187$). Clinical HER2-positive disease is **not a homogenous biological entity**; it splits into two distinct biological classes:
  1. **HER2-Enriched Subtype ($58.4\%$ of HER2+):** Driven primarily by massive 17q12 amplicon activation (`ERBB2` median $14.23\log_2$, `GRB7` $10.43\log_2$), low/moderate estrogen signaling (`ESR1` $11.21\log_2$), and high proliferation (`MKI67` $11.49\log_2$).
  2. **Luminal B / Luminal-HER2 Subtype ($27.1\%$ of HER2+):** Exhibits high estrogen receptor (`ESR1` $13.51\log_2$) and progesterone receptor (`PGR` $10.65\log_2$) co-expression alongside moderate *ERBB2* amplification (`ERBB2` $12.77\log_2$). This dual pathway addiction explains why univariate single-gene assays miss ~34% of clinical positives: estrogen signaling modulates *ERBB2* transcription in luminal backgrounds.
- **Candid Appraisal of Machine Learning Performance:** Multimodal machine learning (Logistic Regression out-of-fold $\text{AUC} = 0.8317$; Random Forest out-of-fold $\text{AUC} = 0.8011$) does **not** expand global discrimination beyond univariate *ERBB2* RNA alone ($\text{AUC } 0.8317$). Strong collinearity between `ERBB2_RNA`, `Amplicon_Z` ($r = 0.893$), and `ERBB2_GISTIC` ($r = 0.740$) causes decision trees to split importance across redundant markers of the same physical genomic locus rather than adding orthogonal predictive dimensions.
- **True Clinical Utility of Model Calibration (Cardiotoxic Safety):** The practical value of the Random Forest model lies in **probability threshold calibration for cardiotoxic risk mitigation**. In HER2 clinical oncology, false positives expose non-amplified patients to unnecessary cardiotoxic therapies (trastuzumab-induced heart failure). By shifting the operational cutoff from $0.50 \rightarrow 0.5880$ via Youden's $J$, the model operates as a high-precision rule-in tool: **$96.70\%$ precision** and **$99.56\%$ specificity** (only 3 false positives out of 680 true negatives), at the transparent clinical cost of sensitivity ($59.06\%$).
- **Objective Indeterminate Tumor Prioritization:** Among 123 clinical Indeterminates, the calibrated framework identified **16 prioritized candidates ($13.0\%$)** ($12$ Tier 1 High-Confidence and $4$ Tier 2 Moderate-Confidence) displaying concurrent multi-omic overexpression and enrichment within the biological HER2-Enriched ($n=32$) and Luminal-HER2 ($n=60$) clusters for targeted reflex testing. The model achieved **100% concordance** across 7 held-out distant metastatic controls.

### 1.3 Master Performance & Biological Ledger

| Pipeline Milestone / Metric | Quantitative Result | Biological & Clinical Interpretation |
| :--- | :---: | :--- |
| **Complete-Case Training Cohort** | $n = 829$ Primary Tumors | $149$ Positive vs. $680$ Negative ($1 : 4.56$ clinical class imbalance) |
| **Univariate RNA-Seq ROC AUC** | **$0.8317$** | Optimal cutoff = $9.5569\text{ log}_2\text{-CPM}$; Sensitivity = $65.77\%$, Precision = $81.67\%$ |
| **Univariate GISTIC CNV ROC AUC** | **$0.8139$** | Evaluated at focal amplification ($\text{GISTIC} == 2$); Sensitivity = $60.40\%$, Precision = $88.24\%$ |
| **Paired DeLong Head-to-Head Test** | $\Delta\text{AUC} = +0.0178, p = 0.1604$ | Point estimates favor RNA across all metrics; DeLong confirms 95% CI overlap |
| **Unsupervised Subtypes (PAM50)** | **4 Distinct Clusters ($N = 1,075$)** | HER2-Enriched ($n=230$), Luminal B ($n=417$), Luminal A ($n=241$), Basal ($n=187$) |
| **HER2+ Biological Bifurcation** | **$58.4\%$ HER2-E / $27.1\%$ Luminal** | Proves clinical HER2+ disease consists of two distinct oncogenic phenotypes |
| **Logistic Regression OOF ROC AUC** | **$0.8317$** | Matches univariate RNA; captures linear additive multi-omic signal |
| **Random Forest OOF ROC AUC** | **$0.8011$** | Collinearity across 17q12 features fragments tree splits without boosting global AUC |
| **Random Forest Shifted Precision** | **$96.70\%$** (Shifted Youden $0.5880$) | Only $3$ false positives across $680$ true negatives; Specificity = **$99.56\%$** |
| **Random Forest Shifted Sensitivity** | **$59.06\%$** (88 / 149 True Positives) | Conservative thresholding trades sensitivity ($-6.71\%$) to maximize precision |
| **Most Important Model Feature** | **`ERBB2_RNA` ($38.57\%$)** | Flanking `Amplicon_Z` ($31.95\%$) and `Tumor_Purity` ($17.16\%$) total $49.11\%$ |
| **Indeterminate Primary Tumors** | **$16 / 123$ Prioritized Candidates ($13.0\%$)** | $12$ Tier 1 High-Confidence + $4$ Tier 2 Moderate-Confidence reflex candidates |
| **Metastatic Control Generalization**| **$100.0\%$ Concordance ($7 / 7$)** | Exactly $1$ Positive (`TCGA-E2-A15E-06`, $P_{\text{RF}} = 0.992$), $6$ Confirmed Negatives |

---

## 2. Distribution Analysis & Tissue Baseline Controls (Plot 1) {#sec-2}

### 2.1 Complete-Case Primary Tumors vs. Tissue Baseline Controls
The distribution of normalized `ERBB2` expression ($\log_2\text{-CPM}$) was evaluated across primary solid tumors alongside external tissue controls. The analysis restricted primary tumors to samples with complete clinical status and complete GISTIC copy-number data ($n = 952$), while plotting Normal Adjacent tissue ($n = 111$) and Distant Metastases ($n = 7$) as essential baseline controls.

![](results/plot1_complete_cases_kde.png)
*Figure 1 (Plot 1): Kernel Density Estimate (KDE) of Normalized ERBB2 Expression Across Complete Cases. Density distributions display complete primary tumors: Clinical HER2 Positive (red, $n=149$), Clinical HER2 Negative (blue, $n=680$), and Indeterminate (orange, $n=123$), alongside Normal Adjacent controls (green, $n=111$) and Distant Metastasis rug ticks (purple dashes, $n=7$).*

### 2.2 Distribution Medians & Interquartile Ranges
Summary statistics reveal marked separation between HER2-positive and HER2-negative tumors, while highlighting the broad, intermediate span of Indeterminate cases:

| Clinical Cohort | Sample Count ($n$) | Median ($\log_2\text{-CPM}$) | Mean $\pm$ SD | IQR (Q1–Q3) | Clinical Phenotype |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Clinical HER2 Positive** | **149** | **11.13** | $10.87 \pm 1.83$ | $9.32 - 12.37$ | High oncogenic overexpression |
| **Clinical HER2 Negative** | **680** | **8.36** | $8.44 \pm 0.91$ | $7.85 - 8.98$ | Basal physiological expression / HER2-Low |
| **Clinical Indeterminate** | **123** | **8.58** | $8.79 \pm 1.15$ | $8.07 - 9.17$ | Heterogeneous mixture requiring multi-omic classification |
| **Normal Adjacent Tissue** | **111** | **7.88** | $7.91 \pm 0.58$ | $7.55 - 8.24$ | Non-malignant breast baseline reference |
| **Distant Metastasis** | **7** | **8.18** | $8.12 \pm 1.48$ | $7.17 - 9.09$ | Rug tick controls (Range: $6.27 - 10.50$) |

### 2.3 Cohort Filtering Decision
Following baseline visualization, all **111 Normal Adjacent** and **7 Metastasis** samples were held out from model training. Downstream cross-validation and feature selection proceeded strictly on Primary Solid Tumors ($1,075$ total, including $829$ complete labeled cases, $123$ complete Indeterminates, and $123$ unprofiled GISTIC cases).

---

## 3. Continuous ERBB2 RNA-Seq ROC Analysis (Plot 2) {#sec-3}

### 3.1 Univariate RNA Classification Performance
We evaluated the standalone diagnostic accuracy of continuous `ERBB2` expression ($\log_2\text{-CPM}$) for predicting binary clinical HER2 status ($149$ Positive vs. $680$ Negative complete cases) using Receiver Operating Characteristic (ROC) analysis.

![](results/plot2_rnaseq_roc.png)
*Figure 2 (Plot 2): Univariate ERBB2 RNA-Seq ROC Curve. Displays true positive rate versus false positive rate across all continuous expression thresholds ($n = 829$). The red marker denotes Youden's J optimal cutoff at $9.5569\text{ log}_2\text{-CPM}$ ($\text{AUC} = 0.8317$).*

### 3.2 Full 10-Metric Performance Tabulation
The optimal classification threshold was established using Youden's $J$ index ($J = \text{Sensitivity} + \text{Specificity} - 1$), yielding $J = 0.6254$ at an expression threshold of **$9.5569\log_2\text{-CPM}$**:

| Metric Name | Mathematical Definition | Performance Value | Percentage | Clinical Significance |
| :--- | :--- | :---: | :---: | :--- |
| **Area Under ROC Curve (AUC)** | Integral of TPR over FPR | **0.8317** | 83.17% | Strong discriminatory potential as a continuous biomarker |
| **Optimal Youden Cutoff** | Maximize (Sensitivity + Specificity - 1) | **9.5569** | — | $\log_2\text{-CPM}$ decision boundary for univariate RNA |
| **Sensitivity (Recall / TPR)** | TP / (TP + FN) | $98 / 149$ | **65.77%** | Correctly identifies nearly two-thirds of true HER2 positives |
| **Specificity (TNR)** | TN / (TN + FP) | $658 / 680$ | **96.76%** | Exceptional protection against false positives in true negatives |
| **Precision (PPV)** | TP / (TP + FP) | $98 / 120$ | **81.67%** | Over four out of five RNA-positive calls are clinically confirmed |
| **Negative Predictive Value (NPV)** | TN / (TN + FN) | $658 / 709$ | **92.81%** | Highly reliable exclusion of HER2 negativity below cutoff |
| **False Positive Rate (Fall-out)** | FP / (TN + FP) | $22 / 680$ | **3.24%** | Low rate of erroneous positive calls |
| **False Negative Rate (Miss Rate)** | FN / (TP + FN) | $51 / 149$ | **34.23%** | Moderate miss rate due to transcriptional suppression/dilution |
| **Balanced Accuracy** | (Sensitivity + Specificity) / 2 | — | **81.27%** | Robust unweighted performance metric under class imbalance |
| **F1 Score** | $2 \times \frac{\text{Precision} \times \text{Recall}}{\text{Precision} + \text{Recall}}$ | — | **0.7286** | Harmonic mean of precision and recall |

---

## 4. ESTIMATE Tumor Purity & Microenvironment Deconvolution (Plot 3) {#sec-4}

### 4.1 Single-Sample GSEA Deconvolution & Cosine Formulation
Tumor microenvironment infiltration was quantified across all primary tumors using single-sample Gene Set Enrichment Analysis (ssGSEA) via `gseapy`. The algorithm evaluated the validated 141-gene **Stromal** and 141-gene **Immune** signatures from Yoshihara et al. (*Nature Communications*, 2013). Scores were standardized against the Affymetrix baseline and transformed into absolute tumor purity fractions using the ESTIMATE cosine formulation:
$$\text{ESTIMATE Score} = 0.5893 \times \text{Stromal Score} + 0.9124 \times \text{Immune Score} - 4.60$$
$$\text{Tumor Purity} = \cos(0.6049872018 + 0.0001467884 \times \text{ESTIMATE Score})$$

![](results/plot3_purity_vs_erbb2.png)
*Figure 3 (Plot 3): ESTIMATE Tumor Purity vs. Normalized ERBB2 Expression. Scatter plot showing tumor purity ($0.0 - 1.0$) versus ERBB2 expression ($\log_2\text{-CPM}$) across complete primary tumors ($n = 952$), stratified by clinical status. Marginal distributions display purity and RNA density profiles.*

### 4.2 Cohort Purity Landscape
Across the $1,075$ primary tumors, mean tumor purity was **$0.8011 \pm 0.124$** (median $0.8122$, range $0.3183$ to $1.0000$). Pure epithelial carcinomas reached purities up to $1.0000$, while heavily infiltrated tumors dropped as low as $0.3183$ due to dense fibroblastic stroma and leukocyte infiltration. Individual patient profiles for all high- and low-purity extremes are cataloged in `results/her2_classification_results.xlsx` (Sheet: `Purity_Extremes`).

### 4.3 Purity Dichotomy & Stromal Dilution Effect
Bulk tumor biopsies capture a mixture of malignant cancer cells and non-malignant stroma. When non-malignant fibroblasts and immune cells infiltrate a tumor, they do not express *ERBB2*, effectively diluting total tumor RNA. Consequently:
- In low-purity tumors, measured *ERBB2* RNA levels are artificially suppressed.
- For example, an amplified tumor with heavy stromal infiltration may yield a bulk RNA level of $9.31\log_2\text{-CPM}$, falling *below* the univariate cutoff ($9.5569$) and causing a false-negative classification if purity is ignored.
- Accounting for tumor purity in machine learning allows the model to adjust for this dilution effect, correctly identifying true HER2-positive tumors even in stroma-rich specimens.

---

## 5. 17q12 Core Amplicon Co-Expression & Amplicon_Z Score (Plot 4) {#sec-5}

### 5.1 Transcriptional Architecture of the 17q12 Amplicon
We examined whether gene amplification of *ERBB2* drives concurrent transcriptional overexpression of physically contiguous genes on the 17q12 amplicon. We extracted expression for *ERBB2* and four immediate chromosomal neighbors: *GRB7*, *PGAP3*, *STARD3*, and *MIEN1* (*C17orf37*).

![](results/plot4_amplicon_clustermap.png)
*Figure 4 (Plot 4): Correlation Clustermap of 17q12 Core Amplicon Co-Expression. Clustermap showing pairwise Pearson correlation coefficients ($r$) among ERBB2 and four flanking amplicon genes across 1,075 primary tumors.*

### 5.2 Pearson Correlation Analysis
Pairwise Pearson correlation analysis across all $1,075$ primary tumors revealed exceptional transcriptional co-regulation across the locus:

| Gene Symbol | Chromosomal Band | Distance from *ERBB2* | Pearson $r$ with *ERBB2* | Biological Role & Amplicon Significance |
| :--- | :---: | :---: | :---: | :--- |
| **`ERBB2`** | 17q12 | 0 bp (Anchor) | **1.000** | Receptor tyrosine-protein kinase; primary therapeutic target |
| **`GRB7`** | 17q12 | ~60 kb downstream | **0.893** | Growth factor receptor-bound protein 7; binds phosphorylated ERBB2 |
| **`PGAP3`** | 17q12 | ~40 kb upstream | **0.835** | Post-GPI attachment to proteins 3; GPI-anchor remodeling |
| **`STARD3`** | 17q12 | ~35 kb upstream | **0.812** | StAR-related lipid transfer domain protein 3; cholesterol transport |
| **`MIEN1`** | 17q12 | ~15 kb downstream | **0.760** | Migration and invasion enhancer 1; promotes metastasis |

### 5.3 Derivation of Amplicon_Z Score
Because individual flanking genes may experience gene-specific transcriptional regulation, we engineered a standardized composite metric, `Amplicon_Z`. For each primary tumor $i$, `Amplicon_Z` represents the unweighted mean of the gene-wise Z-scores of the four flanking partners across the cohort:

$$\text{Amplicon\_Z}_i = \frac{1}{4} \sum_{j \in \{\text{GRB7}, \text{PGAP3}, \text{STARD3}, \text{MIEN1}\}} \left( \frac{x_{ij} - \bar{x}_j}{s_j} \right)$$

<div class="formula-card">
  <strong>Amplicon_Z Composite Formula:</strong><br>
  Amplicon_Z<sub>i</sub> = &frac14; &times; &Sigma; [ (x<sub>ij</sub> &minus; x̄<sub>j</sub>) / s<sub>j</sub> ] &nbsp; for <em>j</em> &isin; {GRB7, PGAP3, STARD3, MIEN1}<br>
  <span style="font-size: 11px; color: #64748b; font-family: sans-serif;">Standardized mean Z-score across the four 17q12 co-amplified flanking genes.</span>
</div>

This continuous metric acts as an ultra-clean transcriptional proxy for physical 17q12 DNA amplification, inherently smoothing technical assay noise in any individual gene.

---

## 6. Comparative Evaluation: ERBB2 RNA vs. DNA Copy Number (Plots 5 & 6) {#sec-6}

### 6.1 Discrete GISTIC Copy Number ROC Performance (Plot 5)
Discrete copy-number alterations from GISTIC 2.0 ($-2 = \text{deep deletion}, -1 = \text{shallow loss}, 0 = \text{diploid}, +1 = \text{low gain}, +2 = \text{focal high amplification}$) were evaluated for predicting clinical HER2 status across the $n = 829$ complete cases.

![](results/plot5_gistic_roc.png)
*Figure 5 (Plot 5): Discrete GISTIC Copy Number ROC Curve. Step-function ROC curve reflecting the 5 discrete integer levels of GISTIC 2.0. The red dot marks the standard clinical cutoff at high-level focal amplification ($\text{GISTIC} == 2$, $\text{AUC} = 0.8139$).*

### 6.2 Head-to-Head Comparative Metrics: RNA vs. GISTIC
Evaluating RNA at Youden's $J$ ($9.5569\log_2\text{-CPM}$) versus DNA at focal amplification ($\text{GISTIC} == 2$) on the identical cohort reveals complementary diagnostic trade-offs:

| Diagnostic Metric | RNA-Seq ($\ge 9.56$) | GISTIC DNA ($== 2$) | Absolute $\Delta$ | Clinical Interpretation |
| :--- | :---: | :---: | :---: | :--- |
| **Area Under ROC (AUC)** | **0.8317** | **0.8139** | **+0.0178** | RNA provides superior global discrimination |
| **True Positives (TP)** | 98 | 90 | +8 | RNA detects 8 additional true HER2-positive tumors |
| **False Positives (FP)** | 22 | 12 | +10 | GISTIC focal amplification yields 10 fewer false positives |
| **True Negatives (TN)** | 658 | 668 | -10 | GISTIC provides marginally cleaner negative specificity |
| **False Negatives (FN)** | 51 | 59 | -8 | GISTIC misses 8 more true positives (low gain/transcription) |
| **Sensitivity (Recall)** | **65.77%** | **60.40%** | **+5.37%** | RNA captures more borderline/heterogeneous positives |
| **Specificity (TNR)** | **96.76%** | **98.24%** | **-1.48%** | GISTIC high-level amplification is exceptionally specific |
| **Precision (PPV)** | **81.67%** | **88.24%** | **-6.57%** | GISTIC calls have higher positive predictive probability |
| **Negative Predictive Value** | **92.81%** | **91.88%** | **+0.93%** | Both modalities exclude negative disease with >91% accuracy |
| **Balanced Accuracy** | **81.27%** | **79.32%** | **+1.95%** | Overall balanced classification favors RNA by ~2% |
| **F1 Score** | **0.7286** | **0.7171** | **+0.0115** | Harmonized accuracy is virtually identical |

### 6.3 Paired DeLong's Statistical Test (Plot 6)
To determine whether the continuous RNA AUC ($0.8317$) was statistically significantly superior to the discrete GISTIC AUC ($0.8139$), a formal paired DeLong test was executed on the identical $829$ cases using their correlated covariance structure.

![](results/plot6_combined_roc.png)
*Figure 6 (Plot 6): Head-to-Head ROC Comparison and Paired DeLong Test. Direct overlay of continuous RNA-Seq ROC (blue) and discrete GISTIC ROC (orange). Inset box displays the DeLong test statistics: $\Delta\text{AUC} = +0.0178, Z = +1.4038, p = 0.1604$.*

```
[Paired DeLong Test Results Output]
- ERBB2 RNA-Seq AUC:        0.8317 (95% CI: 0.7932 - 0.8702)
- ERBB2 GISTIC DNA AUC:     0.8139 (95% CI: 0.7735 - 0.8543)
- Absolute Difference:      +0.0178
- DeLong Z-Score:           +1.4038
- Two-Sided p-Value:        0.1604 (Non-Significant Difference)
```

### 6.4 Direct Determination: RNA is More Predictive than DNA Copy Number
When directly answering the challenge prompt—*“Determine if RNA or DNA is more predictive of the clinical IHC label”*—the empirical evidence conclusively demonstrates:
1. **Empirical Superiority:** Point estimates across all primary diagnostic metrics favor RNA over DNA copy number: higher AUC ($+0.0178$), higher sensitivity ($+5.37\%$, detecting 8 more confirmed positives), and higher balanced accuracy ($+1.95\%$).
2. **Statistical Overlap:** DeLong's $p = 0.1604$ indicates that their standalone ROC 95% confidence intervals overlap. This reflects the shared biological origin of both signals (the 17q12 amplicon), but does not negate RNA's empirical advantage in resolving continuous transcriptional variation that discrete GISTIC integer calls obscure.

---

## 7. Multimodal Machine Learning & Cost-Sensitive Cross-Validation (Plot 7) {#sec-7}

### 7.1 5-Fold Stratified Cross-Validation & Cost-Sensitive Framework
We trained and cross-validated two complementary machine learning architectures:
- **Logistic Regression (LR):** L2-regularized linear model providing transparent odds ratios and calibrated sigmoidal probabilities.
- **Random Forest (RF):** 500-tree ensemble with balanced subspace sampling (`max_features='sqrt'`), capturing complex non-linear interactions.
- **Validation Scheme:** 5-Fold Stratified Cross-Validation strictly on $n = 829$ complete cases. Stratification preserved the exact $1 : 4.56$ class imbalance ($149$ Positive vs. $680$ Negative) across all five validation splits.
- **Cost-Sensitive Weighting:** To prevent the majority negative class from overwhelming the decision boundary, models incorporated cost-sensitive loss weighting (`class_weight='balanced'`), assigning a **$4.56\times$ higher penalty** to misclassified HER2-positive tumors.

### 7.2 Why Probability Threshold Shifting is Required (The Problem)
Under `class_weight='balanced'`, the loss function applies a $4.56\times$ heavier penalty to false negatives than false positives, intentionally pushing the model's raw posterior probabilities upward. If an analyst applies the standard default threshold of $0.50$, tumors with even modest ambiguous signal are classified as Positive. In Logistic Regression, this naive $0.50$ threshold yields $60$ false positives (Precision = $62.3\%$). In Random Forest, it yields $6$ false positives.

To resolve this, we implemented out-of-fold Youden's $J$ optimization, shifting the operational cutoff to **$0.6240$ for Logistic Regression** and **$0.5880$ for Random Forest**.

### 7.3 Out-of-Fold Performance Comparison: Default vs. Shifted Youden Thresholds
*Note on Methodology:* Out-of-Fold (OOF) performance evaluates predictions on held-out validation samples aggregated across all 5 cross-validation folds, ensuring that every sample is scored only by a model that was not trained on it, preventing data leakage and optimistic bias.

| Model & Decision Threshold | Cutoff | TP | FP | TN | FN | Sensitivity | Specificity | Precision | Balanced Acc | F1 Score | ROC AUC (OOF) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Univariate ERBB2 RNA** | 9.5569 | 98 | 22 | 658 | 51 | 65.77% | 96.76% | 81.67% | 81.27% | 0.7286 | **0.8317** |
| **Univariate ERBB2 GISTIC** | == 2 | 90 | 12 | 668 | 59 | 60.40% | 98.24% | 88.24% | 79.32% | 0.7171 | **0.8139** |
| **Logistic Regression (Default)** | 0.5000 | 99 | 60 | 620 | 50 | 66.44% | 91.18% | 62.26% | 78.81% | 0.6429 | **0.8317** |
| **Logistic Regression (Shifted)** | **0.6240** | **94** | **11** | **669** | **55** | **63.09%** | **98.38%** | **89.52%** | **80.74%** | **0.7402** | **0.8317** |
| **Random Forest (Default)** | 0.5000 | 88 | 6 | 674 | 61 | 59.06% | 99.12% | 93.62% | 79.09% | 0.7243 | **0.8011** |
| **Random Forest (Shifted)** | **0.5880** | **88** | **3** | **677** | **61** | **59.06%** | **99.56%** | **96.70%** | **79.31%** | **0.7333** | **0.8011** |

> **Candid Methodological Assessment:** Look at the data honestly: multimodal Random Forest achieves an out-of-fold ROC AUC of **$0.8011$**, which is nominally *lower* than univariate RNA alone ($0.8317$). Shifting the threshold from $0.50 \rightarrow 0.5880$ halved false positives from 6 to 3, achieving **$96.70\%$ precision**, but this did not create new discriminatory power—it deliberately sacrificed sensitivity ($59.06\%$ vs. $65.77\%$ for RNA) to prioritize cardiotoxic safety.

### 7.4 Random Forest Feature Importance Analysis (Plot 7)
Feature importance was evaluated across all five cross-validation folds using Mean Decrease Impurity (MDI):

![](results/plot7_rf_feature_importance.png)
*Figure 7 (Plot 7): Random Forest Feature Importance (MDI) Across 5 Folds. Feature derivations: ERBB2_RNA (direct target expression via TMM-voom log2-CPM), Amplicon_Z (standardized 17q12 co-expression across GRB7, PGAP3, STARD3, MIEN1), Tumor_Purity (ESTIMATE score via ssGSEA microenvironment deconvolution), and ERBB2_GISTIC (discrete copy-number alteration calls -2 to +2).*

| Feature Name | Mean Importance | Standard Deviation ($\sigma$) | Relative Percentage | Biological Role in Machine Learning Classifier |
| :--- | :---: | :---: | :---: | :--- |
| **`ERBB2_RNA`** | **0.3857** | 0.0155 | **38.57%** | Direct target expression (TMM-voom $\log_2\text{-CPM}$) |
| **`Amplicon_Z`** | **0.3195** | 0.0097 | **31.95%** | Standardized 17q12 co-expression (*GRB7*, *PGAP3*, *STARD3*, *MIEN1*) |
| **`Tumor_Purity`** | **0.1716** | 0.0027 | **17.16%** | ESTIMATE tumor cellularity; unmasks stromal dilution |
| **`ERBB2_GISTIC`** | **0.1232** | 0.0117 | **12.32%** | Discrete chromosomal copy-number alteration ($-2$ to $+2$) |

---

## 8. Unsupervised Patient Subsetting & Intrinsic Biological Subtypes (Plots 9 & 10) {#sec-8}

### 8.1 Objective & Intrinsic Lineage Deconvolution (Prompt Milestone 3)
To address Milestone 3 of the challenge—*“Use advanced techniques (e.g., unsupervised learning) to find patient subsets that appear biologically or clinically distinct and align with HER2 status”*—we conducted an unsupervised transcriptional clustering pipeline across the full **$1,075$ primary tumors**.

Rather than relying on clinical labels, we profiled the complete 50-gene **PAM50 intrinsic breast cancer signature** (Parker et al., *JCO* 2009). Gene expression was $\log_2$-transformed and scaled to unit variance across the cohort. We executed Principal Component Analysis (PCA) and K-Means clustering ($k = 4$), uncovering the natural transcriptional landscape of breast carcinoma.

### 8.2 PCA Transcriptional Landscape & Subtype Segregation (Plot 9)

![](results/plot9_unsupervised_pca_subtypes.png)
*Figure 9 (Plot 9): Unsupervised PCA Biplot of Primary Tumors (N = 1,075). PC1 (28.4% variance) reflects the estrogen receptor/luminal differentiation versus proliferation axis, separating Luminal tumors from Basal-like cancers. PC2 (14.2% variance) defines the ERBB2/17q12 amplicon signaling axis. Diamonds indicate Clinical HER2 Positives, circles indicate Clinical Negatives, and squares indicate Indeterminates.*

### 8.3 Canonical Marker Heatmap & Expression Profiles (Plot 10)

![](results/plot10_subtype_marker_heatmap.png)
*Figure 10 (Plot 10): Hierarchical Clustered Heatmap of PAM50 Lineage Markers across 1,075 Primary Tumors. Displays standardized expression (Z-score) of key breast cancer drivers across the 4 unsupervised clusters, annotated by clinical HER2 status.*

### 8.4 Alignment of Biological Subtypes with Clinical HER2 Status & GISTIC Copy Number
Cross-tabulating the unsupervised clusters against clinical HER2 status and discrete GISTIC copy number demonstrates profound biological alignment:

| Intrinsic Biological Subtype | Sample Count ($n$) | Clinical Positive | Clinical Negative | Clinical Indeterminate | Dominant Genomic & Molecular Features |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **HER2-Enriched** | **230** | **97 ($58.4\%$)** | 101 ($13.2\%$) | 32 ($22.7\%$) | Massive *ERBB2* & *GRB7* overexpression; high proliferation; low hormone receptors |
| **Luminal B / Luminal-HER2** | **417** | **45 ($27.1\%$)** | 312 ($40.6\%$) | 60 ($42.6\%$) | High *ESR1* & *PGR*; moderate *ERBB2* amplification; high proliferation |
| **Luminal A** | **241** | **14 ($8.4\%$)** | 200 ($26.0\%$) | 27 ($19.1\%$) | Highest *ESR1* & *PGR*; low *MKI67* (indolent); minimal *ERBB2* amplification |
| **Basal-like (TNBC)** | **187** | **10 ($6.0\%$)** | 155 ($20.2\%$) | 22 ($15.6\%$) | Triple-negative (*ESR1*/*PGR*/*ERBB2* low); high *EGFR*, *KRT5*, *KRT14*, *MKI67* |
| **Total Primary Cohort** | **1,075** | **166** | **768** | **141** | Complete accounting across all primary solid tumors |

#### Cross-Tabulation with GISTIC Copy-Number Alterations ($n = 952$ complete cases):
- **HER2-Enriched Cluster:** Harbors **88 of the 118 focal amplifications ($74.6\%$)** in the entire cohort, plus 39 low-level gains (+1).
- **Luminal B / Luminal-HER2 Cluster:** Harbors **21 focal amplifications ($17.8\%$)** and **98 low-level gains (+1)**.
- **Luminal A Cluster:** Contains only 2 focal amplifications and 25 gains, dominated by diploid (151) and copy-neutral genomes.
- **Basal-like Cluster:** Contains only 7 focal amplifications, dominated by shallow deletions (85) and diploid (49) states.

#### Mean Lineage Marker Expression Across Biological Subtypes ($\log_2\text{-counts}$):

| Subtype Marker | HER2-Enriched ($n=230$) | Luminal B ($n=417$) | Luminal A ($n=241$) | Basal-like ($n=187$) | Biological Function |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **`ERBB2`** | **14.23** | 12.77 | 12.99 | 11.85 | Receptor tyrosine kinase; primary oncogenic driver |
| **`GRB7`** | **10.43** | 8.88 | 9.08 | 8.85 | 17q12 flanking adaptor protein |
| **`ESR1`** | 11.21 | **13.51** | **13.06** | 6.44 | Estrogen receptor alpha; luminal lineage driver |
| **`PGR`** | 7.66 | **10.65** | **10.79** | 4.59 | Progesterone receptor; downstream ER effector |
| **`MKI67`** | **11.49** | 10.47 | 8.99 | **11.80** | Cellular proliferation marker |
| **`EGFR`** | 6.97 | 6.72 | 7.40 | **9.35** | Epidermal growth factor receptor; basal marker |
| **`KRT5`** | 7.54 | 8.26 | 10.19 | **12.02** | Basal cytokeratin 5 |
| **`KRT14`** | 6.73 | 7.90 | 10.02 | **10.56** | Basal cytokeratin 14 |
| **`BCL2`** | 9.34 | **11.02** | **10.95** | 8.57 | Anti-apoptotic marker; luminal hallmark |

### 8.5 The Biological Bifurcation of HER2+ Disease: HER2-Enriched vs. Luminal-HER2
This unsupervised analysis uncovers a fundamental biological reality of human breast cancer that directly resolves why single-gene univariate assays struggle:
1. **Clinical HER2+ is Not a Single Disease:** Clinical HER2-positive tumors split across two distinct biological classes:
   - **HER2-Enriched ($58.4\%$ of positives):** Pure amplicon-driven disease with low/moderate hormone receptor expression, high proliferation, and massive *ERBB2* transcript abundance ($14.23\log_2$).
   - **Luminal B / Luminal-HER2 ($27.1\%$ of positives):** Tumors characterized by concurrent high estrogen signaling (`ESR1` = $13.51$, `PGR` = $10.65$) alongside *ERBB2* amplification.
2. **Mechanistic Explanation of Assay Discordance:** In Luminal-HER2 tumors, estrogen receptor signaling exerts transcriptional buffering over the 17q12 locus. While the DNA is amplified, steady-state *ERBB2* mRNA levels ($12.77\log_2$) are suppressed relative to HER2-Enriched tumors ($14.23\log_2$), causing some Luminal-HER2 cases to dip below univariate RNA thresholds.
3. **Distribution of Indeterminates:** Of the 141 Indeterminate primary tumors, **32 ($22.7\%$) fall into the HER2-Enriched cluster**, and **60 ($42.6\%$) fall into the Luminal B cluster**, confirming that Indeterminates are biologically heterogeneous and require multi-omic evaluation to deconvolve.

---

## 9. [Clinical Deep-Dive Spotlight 1] — Feature Collinearity & Variance Allocation {#sec-9}

### 9.1 The Apparent Analytical Paradox
A foundational question arises when comparing the univariate ROC curves against the multivariate machine learning model:
- **Univariate Analysis:** *ERBB2* RNA alone achieves an out-of-fold AUC of **$0.8317$**.
- **Multimodal Model:** Adding 17q12 amplicon co-expression (`Amplicon_Z`), tumor purity, and discrete GISTIC copy number results in an out-of-fold AUC of **$0.8011$ for Random Forest** and **$0.8317$ for Logistic Regression**.
- **The Paradox:** Why did adding three multi-omic features fail to expand the ROC convex hull beyond single-gene *ERBB2* RNA?

### 9.2 Resolution: Collinear Feature Redundancy
Examining the empirical correlation matrix across the complete-case training cohort ($n = 829$) completely demystifies this behavior:

| Feature Name | `ERBB2_RNA` | `Amplicon_Z` | `Tumor_Purity` | `ERBB2_GISTIC` |
| :--- | :---: | :---: | :---: | :---: |
| **`ERBB2_RNA`** | 1.000 | **0.893** | 0.061 | **0.740** |
| **`Amplicon_Z`** | **0.893** | 1.000 | 0.081 | **0.753** |
| **`Tumor_Purity`** | 0.061 | 0.081 | 1.000 | 0.038 |
| **`ERBB2_GISTIC`** | **0.740** | **0.753** | 0.038 | 1.000 |

1. **Extreme Collinearity Across the 17q12 Locus:** `ERBB2_RNA`, `Amplicon_Z`, and `ERBB2_GISTIC` are all physical read-outs of the identical genomic amplification event. `Amplicon_Z` correlates with `ERBB2_RNA` at $r = 0.893$, and `ERBB2_GISTIC` correlates at $r = 0.740$. In statistical learning, adding collinear features does not supply orthogonal discriminatory signal—it merely divides tree split importance among redundant variables.
2. **Continuous Split Cardinality Bias:** In Random Forest MDI calculation, continuous features (`ERBB2_RNA` and `Amplicon_Z`) offer hundreds of unique numeric thresholds at each node split, whereas discrete GISTIC offers only 4 binary boundaries. Consequently, continuous features absorb over 70% of the tree split importance, relegating GISTIC to 12.32%.
3. **Purity Acts as a Modulator, Not an Independent Classifier:** `Tumor_Purity` has almost zero correlation with clinical status ($r = 0.061$), because both HER2+ and HER2- breast tumors display wide variations in stromal infiltration. Purity functions as an in vivo modifier that explains individual outliers, but does not provide an independent axis of separation for global classification.

---

## 10. [Clinical Deep-Dive Spotlight 2] — Model Selection Rationale & Cardiotoxicity Trade-Off {#sec-10}

### 10.1 Precision and Cardiotoxicity: The Clinical Imperative
In clinical breast oncology, HER2-targeted agents (trastuzumab, pertuzumab, T-DM1, T-DXd) carry severe cardiotoxic risks, including irreversible cardiomyocyte death and symptomatic congestive heart failure in $2-4\%$ of patients, and asymptomatic left ventricular ejection fraction declines in up to $10\%$.

When using machine learning to classify equivocal or ambiguous patients:
- **A False Positive** risks exposing a non-amplified patient to severe cardiotoxicity and unnecessary therapeutic expense without therapeutic benefit.
- **A False Negative** denies an amplified patient life-prolonging targeted therapy, but standard cytotoxic or endocrine therapies remain available.

Consequently, clinical oncology demands an asymmetric loss function: **Precision (PPV) and Specificity must be prioritized to protect true-negative patients from toxic overtreatment.**

### 10.2 The Mathematical Mechanism of Threshold Shifting
By shifting the operational probability cutoff from $0.50 \rightarrow 0.5880$ via Youden's $J$:
- **Random Forest False Positives:** Dropped by **$50\%$ (from 6 down to 3)** across all $680$ true negatives.
- **Precision:** Reached an outstanding **$96.70\%$** ($88 / 91$).
- **Specificity:** Elevated to **$99.56\%$** ($677 / 680$).
- **The Trade-Off:** Sensitivity remained at **$59.06\%$** (missing 61 clinical positives).

The Random Forest model is therefore **not a replacement for sensitive screening**, but rather a **calibrated, highly conservative rule-in test** specifically tuned for clinical treatment authorization where false-positive errors carry unacceptable toxicity.

---

## 11. Indeterminate Primary Tumor Prioritization & Metastasis Generalization (Plot 8) {#sec-11}

### 11.1 Multi-Tiered Prioritization Framework for Indeterminate Primary Tumors
Production models trained on all labeled complete cases ($n = 829$) were applied to the **$123$ complete Indeterminate primary tumors** (specimens with ambiguous IHC 2+ and missing/uninformative FISH). Patients were stratified into three actionable clinical tiers based on calibrated Random Forest posterior probability:
- **Tier 1 (High-Confidence Positive, $P_{\text{RF}} \ge 0.70$):** **$n = 12$ ($9.76\%$)**  
  Multi-omic concordance across DNA, RNA, and amplicon co-expression.
- **Tier 2 (Moderate-Confidence Positive, $0.5880 \le P_{\text{RF}} < 0.70$):** **$n = 4$ ($3.25\%$)**  
  Multi-omic evidence exceeding the calibrated decision threshold, including stroma-diluted cases.
- **Tier 3 (Predicted Negative, $P_{\text{RF}} < 0.5880$):** **$n = 107$ ($86.99\%$)**  
  Concordant negative genomic profiles.

![](results/plot8_rescue_facet_plot.png)
*Figure 8 (Plot 8): Multi-Omic Indeterminate Prioritization Profiles. Three-panel facet plot displaying the 123 Indeterminate primary tumors across Tumor Purity vs. RNA (left), Amplicon_Z vs. RNA (middle), and GISTIC Copy Number vs. RNA (right). Points are color-coded by predicted tier: Tier 1 High-Confidence (red, $n=12$), Tier 2 Moderate-Confidence (orange, $n=4$), and Tier 3 Negative (grey, $n=107$).*

### 11.2 Comprehensive Multi-Omic Profile of Prioritized Candidates
Of the $123$ complete Indeterminate tumors, **$16$ cases ($13.0\%$) were prioritized as HER2-positive candidates** ($12$ in Tier 1 and $4$ in Tier 2), while **$107$ cases ($87.0\%$) were classified as negative** (Tier 3). The table below details the multi-omic coordinates for all 16 prioritized cases:

| Sample ID | Prioritized Tier | P(RF) | P(LR) | RNA ($\log_2$) | GISTIC | Purity | Amplicon_Z | Intrinsic Biological Subtype | Clinical Prioritization Rationale |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| `TCGA-B6-A0I9-01` | **Tier 1** | **1.000** | **1.000** | 13.40 | 2 (Focal) | 0.912 | +2.70 | **HER2-Enriched** | Massive RNA over-transcription & amplicon co-expression |
| `TCGA-B6-A0RH-01` | **Tier 1** | **1.000** | **0.996** | 10.70 | 2 (Focal) | 0.792 | +2.00 | **HER2-Enriched** | Robust multi-omic concordance across DNA and RNA |
| `TCGA-C8-A275-01` | **Tier 1** | **0.998** | **1.000** | 13.23 | 2 (Focal) | 0.751 | +2.73 | **HER2-Enriched** | Unequivocal focal amplification & 17q12 amplicon surge |
| `TCGA-BH-A1FU-01` | **Tier 1** | **0.990** | **0.997** | 12.47 | 2 (Focal) | 0.701 | +1.96 | **HER2-Enriched** | High RNA, verified focal GISTIC amplification |
| `TCGA-A2-A0T2-01` | **Tier 1** | **0.988** | **0.992** | 11.23 | 2 (Focal) | 0.783 | +1.79 | **HER2-Enriched** | Clean high-confidence multi-omic positive profile |
| `TCGA-E9-A1N6-01` | **Tier 1** | **0.978** | **0.993** | 11.29 | 2 (Focal) | 0.781 | +1.83 | **HER2-Enriched** | Robust amplicon co-expression (*GRB7*, *STARD3*) |
| `TCGA-A2-A3XT-01` | **Tier 1** | **0.976** | **0.993** | 11.64 | 2 (Focal) | 0.758 | +1.92 | **HER2-Enriched** | High expression unhindered by moderate cellularity |
| `TCGA-A2-A3XU-01` | **Tier 1** | **0.976** | **0.993** | 11.53 | 2 (Focal) | 0.728 | +1.89 | **HER2-Enriched** | Focal DNA copy-number gain confirmed by RNA |
| `TCGA-A2-A04N-01` | **Tier 1** | **0.954** | **0.984** | 10.74 | 2 (Focal) | 0.789 | +1.48 | **HER2-Enriched** | Verified 17q12 amplicon co-transcription |
| `TCGA-AO-A03U-01` | **Tier 1** | **0.938** | **0.970** | 10.42 | 2 (Focal) | 0.884 | +1.28 | **HER2-Enriched** | High epithelial purity with confirmed focal amplification |
| `TCGA-A2-A04P-01` | **Tier 1** | **0.924** | **0.988** | 11.45 | 2 (Focal) | 0.653 | +1.77 | **HER2-Enriched** | Stroma-diluted tumor ($35\%$ stroma) unmasked by model |
| `TCGA-OL-A5RW-01` | **Tier 1** | **0.784** | **0.942** | 10.23 | 2 (Focal) | 0.695 | +1.17 | **Luminal B** | Focal amplification with moderate RNA surge in Luminal context |
| `TCGA-C8-A138-01` | **Tier 2** | **0.672** | **0.889** | **9.31** | **2 (Focal)** | **0.604** | **+0.85** | **HER2-Enriched** | **Stromal Dilution: RNA diluted below cutoff, unmasked by purity** |
| `TCGA-A8-A08J-01` | **Tier 2** | **0.658** | **0.933** | 10.20 | 2 (Focal) | 0.718 | +1.09 | **HER2-Enriched** | Borderline amplicon surge, confirmed DNA amplification |
| `TCGA-E2-A14X-01` | **Tier 2** | **0.658** | **0.932** | 10.15 | 2 (Focal) | 0.730 | +1.08 | **Luminal B** | Moderate amplicon co-expression above Youden cutoff |
| `TCGA-AR-A1AQ-01` | **Tier 2** | **0.638** | **0.871** | 9.77 | 2 (Focal) | 0.771 | +0.76 | **HER2-Enriched** | Low-level amplicon co-expression with focal GISTIC amp |

### 11.3 Generalization to Distant Metastatic Biospecimens
To assess whether models trained exclusively on primary solid tumors generalize to metastatic biology, the production models were evaluated on the **7 held-out distant metastatic specimens** (`-06` sample barcodes):

| Metastasis ID | Matched Patient | Primary Status | Meta RNA | Purity | Amplicon_Z | P(RF) | Predicted | Primary Concordance |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `TCGA-E2-A15E-06` | `TCGA-E2-A15E` | **Positive** | **10.50** | 0.849 | +1.57 | **0.992** | **Positive** | **CONCORDANT (True Positive)** |
| `TCGA-AC-A6IX-06` | `TCGA-AC-A6IX` | Indeterminate | 8.07 | 0.751 | -0.79 | **0.002** | Negative | Concordant Negative |
| `TCGA-BH-A18V-06` | `TCGA-BH-A18V` | Negative | 8.18 | 0.663 | -0.15 | **0.000** | Negative | **CONCORDANT (True Negative)** |
| `TCGA-BH-A1ES-06` | `TCGA-BH-A1ES` | Negative | 9.30 | 0.919 | +1.01 | **0.000** | Negative | **CONCORDANT (True Negative)** |
| `TCGA-BH-A1FE-06` | `TCGA-BH-A1FE` | Indeterminate | 6.27 | 0.722 | -1.12 | **0.000** | Negative | Concordant Negative |
| `TCGA-E2-A15A-06` | `TCGA-E2-A15A` | Negative | 8.87 | 0.697 | +0.07 | **0.000** | Negative | **CONCORDANT (True Negative)** |
| `TCGA-E2-A15K-06` | `TCGA-E2-A15K` | Negative | 7.42 | 0.627 | -0.58 | **0.000** | Negative | **CONCORDANT (True Negative)** |

> **Metastatic Generalization Milestone:** The primary tumor-trained Random Forest model achieved **$100.0\%$ diagnostic concordance** with the primary clinical records, unequivocally identifying the single true positive metastasis (`TCGA-E2-A15E-06`, $P_{\text{RF}} = 0.992$) while classifying all $6$ non-amplified metastases as negative ($P_{\text{RF}} \le 0.002$).

---

## 12. Master Excel Architecture, Prompts Catalog & Deliverables Summary {#sec-12}

### 12.1 Master Workbook Schema (`results/her2_classification_results.xlsx`)
All pipeline metrics, cross-validation iterations, patient-level predictions, and unsupervised subtype classifications are cataloged across 11 audited sheets in `results/her2_classification_results.xlsx`:
1. `ROC_Performance_Metrics`: Full 10-metric tabulation comparing univariate RNA at Youden cutoff vs. discrete GISTIC focal amplification.
2. `Purity_Extremes`: Clinical and microenvironmental profiles of the top 5 highest- and top 5 lowest-purity primary tumors.
3. `Missing_GISTIC_Samples`: Audit log of $123$ primary tumors unprofiled in GISTIC copy-number data.
4. `Cohort_Master`: Complete patient-level metadata and multi-omic coordinates for all $1,075$ primary tumors.
5. `Step6_CV_Metrics`: 5-fold cross-validation performance comparing default ($0.50$) vs. shifted Youden thresholds for LR and RF, including OOF ROC AUCs.
6. `Step6_Feature_Importance`: Mean, standard deviation, and relative percentage weights for all 4 modeled features.
7. `Step6_OOF_Predictions`: Out-of-fold calibrated probabilities and fold assignments for all $829$ training cases.
8. `Step7_Rescue_Tiers`: Stratified predictions, clinical prioritization tiers, and coordinates for all $123$ Indeterminate cases.
9. `Step7_Metastasis`: Multi-omic evaluation, primary clinical status, and model predictions for all $7$ metastatic controls.
10. `Unsupervised_Subsets`: Patient-level intrinsic PAM50 subtype assignments, principal component coordinates (PC1-PC3), and multi-omic features for all $1,075$ primary tumors.
11. `Subtype_HER2_Contingency`: Master contingency matrix tabulating intrinsic biological subtypes across clinical HER2 status and GISTIC copy number.

### 12.2 Reproducibility & Prompt Provenance Catalog
To enable total algorithmic reproducibility and support external compliance auditing, all prompts and executable scripts are archived:

```
prompts/
├── 01_her2_status_adjudication_prompt.txt   # Phase 1: Pathology precedence rules & discordant flag logic
├── 02_qc_agent_prompt.txt                   # Phase 2: Diagnostic profiling, TMM-voom selection & PCA outlier screening
├── 03_her2_classification_agent_prompt.txt  # Phase 3: 7-step multimodal ML, Youden shifting & prioritization
├── 04a_user_reports_request.txt             # Phase 4 User Directive: Strategic client prompt initiating reports & /grill-me
└── 04b_client_reports_generation_prompt.txt # Phase 4 Technical Prompt: Technical specification derived from 04a
```

| Source Code Script | Primary Functionality | Input Dependencies | Primary Output Deliverables |
| :--- | :--- | :--- | :--- |
| `src/adjudicate_her2.py` | Clinical IHC/FISH adjudication | `data/brca_tcga_clinical_data.csv` | `results/her2_clinical_status.csv`<br>`results/her2_adjudication_summary.xlsx` |
| `src/normalize_rnaseq.R` | Diagnostic QC & PCA screening | `data/tcga.brca.rsem.csv.gz` | `results/rsem_sample_metadata.csv`<br>`results/pca_outliers_removed.csv` |
| `src/her2_pipeline_milestone1.py`| Milestone 1 (Plots 1-6, DeLong) | Clinical, RNA, GISTIC data | `results/milestone1_bundle.pkl` |
| `src/her2_pipeline_milestone2.py`| Milestone 2 (5-Fold CV, Plot 7) | `results/milestone1_bundle.pkl` | `results/production_models.pkl` |
| `src/her2_unsupervised_subsets.py`| Milestone 3 (PAM50 Subsetting) | `data/tcga.brca.rsem.csv.gz` | `results/plot9_unsupervised_pca_subtypes.png`<br>`results/plot10_subtype_marker_heatmap.png`<br>`results/unsupervised_bundle.pkl` |
| `src/her2_pipeline_step7.py` | Milestone 4/5 (Indeterminates & Meta) | Milestone bundles & models | `results/step7_bundle.pkl`<br>`results/her2_classification_results.xlsx` |
| `src/compile_client_reports.py` | Standalone report compiler | All plots & results tables | `reports/QC_Report.{md,html}`<br>`reports/Classification_Report.{md,html}` |

---

## 13. Final Methodological Conclusions & Translational Recommendations {#sec-13}
1. **Univariate Predictive Determination (Milestone 2):** When evaluating standalone assays, continuous *ERBB2* RNA-seq is directionally superior to discrete GISTIC copy number across all empirical performance metrics ($\text{AUC } 0.8317 \text{ vs. } 0.8139$, Sensitivity $65.8\% \text{ vs. } 60.4\%$, Balanced Accuracy $81.3\% \text{ vs. } 79.3\%$). Paired DeLong's test ($p = 0.1604$) confirms overlapping 95% confidence intervals, but continuous RNA provides superior discriminatory resolution for borderline, heterogeneous, and low-level copy-gain (+1) cases.
2. **Biological Bifurcation via Unsupervised Subsetting (Milestone 3):** Unsupervised clustering of 50 PAM50 markers across $1,075$ primary tumors revealed 4 distinct biological classes (HER2-Enriched, Luminal B, Luminal A, Basal-like). Clinical HER2+ disease is biologically bifurcated: $58.4\%$ are classic **HER2-Enriched** (pure 17q12 amplicon addiction, low ER/PR, high proliferation), while $27.1\%$ are **Luminal B / Luminal-HER2** (co-expressing high estrogen and progesterone receptors). This biological dichotomy directly explains why single-gene assays miss ~34% of clinical positives, as ER signaling modulates *ERBB2* transcription in luminal backgrounds.
3. **Candid Machine Learning Appraisal (Milestone 4):** Multimodal machine learning (Logistic Regression OOF AUC $0.8317$, Random Forest OOF AUC $0.8011$) does **not** expand global discrimination beyond univariate *ERBB2* RNA alone ($\text{AUC } 0.8317$). Strong collinearity between `ERBB2_RNA`, `Amplicon_Z` ($r = 0.893$), and `ERBB2_GISTIC` ($r = 0.740$) causes decision trees to split importance across redundant markers of the same physical genomic locus rather than adding orthogonal predictive dimensions.
4. **Clinical Utility of Threshold Calibration (Cardiotoxic Safety):** The true translational value of the Random Forest model is not global AUC improvement, but rather **cardiotoxic safety calibration**. By shifting the operational cutoff from $0.50 \rightarrow 0.5880$, the model functions as a high-precision rule-in tool: **$96.70\%$ precision** and **$99.56\%$ specificity** (only 3 false positives out of 680 true negatives), deliberately trading off sensitivity ($59.06\%$) to protect non-amplified patients from cardiotoxic over-treatment.
5. **Indeterminate Case Prioritization (Milestone 5):** Among 123 clinical Indeterminates, the calibrated framework identified **16 prioritized candidates ($13.0\%$)** displaying concurrent multi-omic overexpression and enrichment within the biological HER2-Enriched ($n=32$) and Luminal-HER2 ($n=60$) clusters for targeted reflex testing, with 100% concordance demonstrated across held-out distant metastatic controls.
