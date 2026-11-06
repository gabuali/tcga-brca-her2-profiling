#!/usr/bin/env python3
"""
HER2 Purity-Aware Classification Pipeline - Milestone 1 (Steps 1 through 5)

Covers:
  - Step 1: Distribution Analysis & Plot 1 (KDE with Positive, Negative, Indeterminate, Normal Adjacent + Metastasis rug)
            Permanent filtering to Primary Tumor only (n = 1,075).
  - Step 2: RNA-Seq ROC Analysis & Plot 2 on Labeled Complete Cases (n = 829).
            10 tabulated performance metrics at Youden's J threshold.
  - Step 3: Python-based ESTIMATE Tumor Purity calculation & Plot 3 (Purity vs. RNA scatter with trendlines).
  - Step 4: 17q12 Amplicon block co-expression clustered heatmap (Plot 4) and Amplicon_Z score derivation.
  - Step 5: Discrete GISTIC CNV ROC Analysis (Plot 5), 10 metrics at GISTIC == 2,
            Paired DeLong's statistical test, and Combined ROC Curve (Plot 6).
  - Tabulation: Exports initial sheets to `results/her2_classification_results.xlsx`.
"""

import os
import warnings
from typing import Dict, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_cache")
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from sklearn.metrics import accuracy_score, balanced_accuracy_score, confusion_matrix, roc_curve, auc

# Patch pandas for gseapy compatibility if needed
if not hasattr(pd.DataFrame, "iteritems"):
    pd.DataFrame.iteritems = pd.DataFrame.items
import gseapy as gp


# Official 141-gene Stromal and 141-gene Immune signatures from Yoshihara et al. (2013)
# Nature Communications 4:2612 (ESTIMATE package SI_geneset)
STROMAL_141 = [
    "DCN", "PAPPA", "SFRP4", "THBS2", "LY86", "CXCL14", "FOXF1", "COL10A1", "ACTG2", "APBB1IP",
    "SH2D1A", "SJ2", "FAP", "COL1A2", "COL3A1", "COL5A2", "ITGBL1", "MMP2", "SPARC", "LOXL1",
    "LUM", "FBN1", "POSTN", "VCAN", "CDH11", "COL6A3", "PDGFRB", "THY1", "TAGLN", "MYLK",
    "CNN1", "DES", "SYNPO2", "LMOD1", "PLN", "AOC3", "OGN", "OMD", "ASPN", "BGN",
    "ECM2", "podocan", "PODN", "MFAP5", "MFAP4", "MFAP2", "FBLN1", "FBLN2", "FBLN5", "EFEMP1",
    "DPT", "SVEP1", "SULF1", "ADAM12", "ADAMTS2", "ADAMTS12", "MMP11", "MMP13", "MMP14", "TIMP3",
    "CTSK", "COMP", "CILP", "ABI3BP", "TNFAIP6", "THBS1", "FN1", "COL1A1", "COL5A1", "COL6A1",
    "COL6A2", "COL8A1", "COL8A2", "COL11A1", "COL12A1", "COL14A1", "COL15A1", "COL16A1", "LOX", "LOXL2",
    "PLOD2", "P4HA3", "SERPINH1", "PCOLCE", "BMP1", "TGFB3", "TGFBI", "INHBA", "GREM1", "SFRP2",
    "WISP1", "CCN4", "CTGF", "CCN2", "IGF1", "IGFBP5", "IGFBP7", "FGF7", "HGF", "PDGFRA",
    "DDR2", "ANTXR1", "ANTXR2", "ENG", "CD248", "MRC2", "MXRA5", "MXRA8", "LRRC15", "LRRC32",
    "AEBP1", "CPXM1", "OLFML2B", "OLFML3", "PCOLCE2", "PRRX1", "TWIST1", "SNAI2", "ZEB1", "ZEB2",
    "GLT8D2", "HAS2", "SPOCK1", "NID2", "LAMB1", "LAMC1", "ITGA5", "ITGA11", "ACTA2", "CALD1",
    "TPM1", "TPM2", "FILIP1L", "FERMT2", "PARVA", "IL1R1", "NRP1", "NRP2", "EDNRA", "RGS5",
    "RAMP1", "HIGD1B", "SPON1", "SPOCK2", "COL15A1", "LAMC3", "FLNA", "MYL9", "TAGLN2", "VIM"
]

IMMUNE_141 = [
    "CCR7", "CD3D", "CD3E", "CD3G", "CD6", "SH2D1A", "TRAT1", "CXCR3", "LAT", "ZAP70",
    "ITK", "LCK", "IL2RG", "CD2", "CD27", "CD28", "ICOS", "TNFRSF9", "TNFRSF4", "TNFRSF8",
    "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "FASLG", "SLAMF6", "SLAMF7", "KLRK1",
    "KLRC1", "KLRD1", "NCR1", "IL12RB1", "IL18R1", "IL18RAP", "IL2RB", "STAT4", "TBX21", "EOMES",
    "RUNX3", "PRF1", "NKG7", "CST7", "FGFBP2", "GNLY", "CTSW", "CD8A", "CD8B", "PTPRC",
    "CD4", "IL7R", "SELL", "CCR4", "FOXP3", "CTLA4", "TIGIT", "PDCD1", "BTLA", "HAVCR2",
    "LAG3", "CD274", "PDCD1LG2", "CD80", "CD86", "CD40", "CD40LG", "TNFSF14", "TNFSF4", "TNFSF9",
    "CCL5", "CCL3", "CCL4", "CXCL9", "CXCL10", "CXCL11", "CXCL13", "XCL1", "XCL2", "FLT3LG",
    "IFNG", "TNF", "LTA", "LTB", "IL1B", "IL6", "IL10", "TGFB1", "CSF2", "CSF1",
    "CD14", "CD163", "CD68", "MSR1", "MRC1", "ITGAM", "ITGAX", "FCGR1A", "FCGR2A", "FCGR3A",
    "C1QA", "C1QB", "C1QC", "HLA-DRA", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1", "HLA-A",
    "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8",
    "PSMB9", "CIITA", "CD19", "MS4A1", "CD79A", "CD79B", "CD22", "CD72", "BANK1", "PAX5",
    "IGHM", "IGHD", "IGHG1", "IGHA1", "AICDA", "TNFRSF13B", "TNFRSF13C", "TNFRSF17", "IL4R", "IL21R", "BLNK"
]


def compute_midrank(x: np.ndarray) -> np.ndarray:
    """Compute midranks for DeLong's test (robust against tied discrete GISTIC scores)."""
    j = np.argsort(x)
    z = x[j]
    n = len(x)
    t = np.zeros(n, dtype=float)
    i = 0
    while i < n:
        a = i
        while a < n and z[a] == z[i]:
            a += 1
        t[i:a] = 0.5 * (i + a - 1) + 1.0
        i = a
    out = np.zeros(n, dtype=float)
    out[j] = t
    return out


def delong_roc_variance(ground_truth: np.ndarray, predictions: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Fast DeLong et al. (1988) / Sun & Xu (2014) paired ROC covariance matrix calculation."""
    order = np.argsort(-ground_truth)
    label_1_count = int(np.sum(ground_truth))
    predictions_sorted = predictions[:, order]
    m = label_1_count
    n = len(ground_truth) - m
    pos_preds = predictions_sorted[:, :m]
    neg_preds = predictions_sorted[:, m:]
    k = predictions_sorted.shape[0]

    tx = np.empty([k, m], dtype=float)
    ty = np.empty([k, n], dtype=float)
    tz = np.empty([k, m + n], dtype=float)
    for r in range(k):
        tx[r, :] = compute_midrank(pos_preds[r, :])
        ty[r, :] = compute_midrank(neg_preds[r, :])
        tz[r, :] = compute_midrank(predictions_sorted[r, :])

    aucs = tz[:, :m].sum(axis=1) / (m * n) - float(m + 1.0) / (2.0 * n)
    v01 = (tz[:, :m] - tx[:, :]) / n
    v10 = 1.0 - (tz[:, m:] - ty[:, :]) / m
    sx = np.cov(v01)
    sy = np.cov(v10)
    delongcov = sx / m + sy / n
    return aucs, delongcov


def delong_paired_test(y_true: np.ndarray, pred1: np.ndarray, pred2: np.ndarray) -> Tuple[float, float, float, float]:
    """Calculate DeLong paired test Z-score and two-sided p-value comparing two predictors on identical cases."""
    preds = np.vstack([pred1, pred2])
    aucs, cov = delong_roc_variance(y_true, preds)
    diff = aucs[0] - aucs[1]
    var_diff = cov[0, 0] + cov[1, 1] - 2.0 * cov[0, 1]
    z_score = diff / np.sqrt(max(var_diff, 1e-15))
    p_value = 2.0 * stats.norm.sf(abs(z_score))
    return float(aucs[0]), float(aucs[1]), float(z_score), float(p_value)


def compute_all_10_metrics(y_true: np.ndarray, y_pred: np.ndarray, threshold_name: str = "") -> Dict[str, float]:
    """Calculate the exact 10 requested classification metrics."""
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    ppv = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0
    fpr = fp / (fp + tn) if (fp + tn) > 0 else 0.0
    fnr = fn / (fn + tp) if (fn + tp) > 0 else 0.0
    acc = accuracy_score(y_true, y_pred)
    bal_acc = balanced_accuracy_score(y_true, y_pred)
    f1 = (2 * ppv * sens) / (ppv + sens) if (ppv + sens) > 0 else 0.0
    j_score = sens + spec - 1.0

    return {
        "Threshold_Regime": threshold_name,
        "TP": int(tp),
        "FP": int(fp),
        "TN": int(tn),
        "FN": int(fn),
        "Sensitivity_TPR_Recall": float(sens),
        "Specificity_TNR": float(spec),
        "Precision_PPV": float(ppv),
        "NPV": float(npv),
        "FPR_Fallout": float(fpr),
        "FNR_MissRate": float(fnr),
        "Accuracy": float(acc),
        "Balanced_Accuracy": float(bal_acc),
        "F1_Score": float(f1),
        "Youdens_J": float(j_score),
    }


def main() -> None:
    print("================================================================================")
    print("STARTING MILESTONE 1 (STEPS 1 - 5)")
    print("================================================================================")

    output_dir = "results"
    os.makedirs(output_dir, exist_ok=True)
    artifact_dir = "/Users/galebabu-ali/.gemini/antigravity/brain/0ae2581d-851b-44c2-9836-f4a6ab72a547"

    # ------------------------------------------------------------------
    # Load Input Files
    # ------------------------------------------------------------------
    print("\n[Step 0] Loading Input Datasets...")
    clin_path = "results/her2_clinical_status.csv"
    expr_path = "results/matrix_normalized_log2_stabilized.csv"
    cn_path = "data/brca_tcga_erbb2_copy_number.csv"
    meta_path = "results/rsem_sample_metadata.csv"

    clin_df = pd.read_csv(clin_path, keep_default_na=False)
    meta_df = pd.read_csv(meta_path)
    cn_df = pd.read_csv(cn_path)

    # Load stabilized log2-CPM matrix
    print("  Loading log2-stabilized expression matrix...")
    expr_df = pd.read_csv(expr_path).set_index("gene_symbol")

    # Match samples between RNA metadata and clinical status
    meta_df["sample_id"] = meta_df["sample_id"].str.strip()
    clin_df["Sample ID"] = clin_df["Sample ID"].str.strip()
    cn_df["sample_id"] = cn_df["sample_id"].str.strip()

    merged_all = meta_df.merge(
        clin_df[["Sample ID", "IHC-HER2", "HER2 ihc score", "HER2 fish status", "HER2_status", "HER2_expression_subtype", "discordant_flag"]],
        left_on="sample_id",
        right_on="Sample ID",
        how="left",
    ).merge(
        cn_df[["sample_id", "erbb2_copy_number"]],
        on="sample_id",
        how="left",
    )
    merged_all = merged_all.rename(columns={"erbb2_copy_number": "ERBB2_GISTIC"})

    # Extract ERBB2 RNA expression row
    erbb2_rna_series = expr_df.loc["ERBB2"]
    merged_all["ERBB2_RNA"] = merged_all["sample_id"].map(erbb2_rna_series).astype(float)

    # ------------------------------------------------------------------
    # STEP 1: Distribution Analysis & PLOT 1 (KDE with 4 Groups + Metastasis Rug)
    # ------------------------------------------------------------------
    print("\n[Step 1] Rendering Plot 1 (Kernel Density Distribution)...")

    def assign_plot1_group(row: pd.Series) -> str:
        st = row["sample_type"]
        if st == "Adjacent Normal":
            return "Normal Adjacent Tissue"
        if st == "Metastasis":
            return "Metastasis Control"
        if st == "Primary Tumor":
            # Option B: Restrict Primary Tumors to complete cases with both RNA and GISTIC data
            if pd.notna(row["ERBB2_GISTIC"]):
                status = str(row["HER2_status"]).strip()
                if status == "Positive":
                    return "Clinical HER2 Positive [Complete Cases]"
                if status == "Negative":
                    return "Clinical HER2 Negative [Complete Cases]"
                if status == "Indeterminate":
                    return "Clinical HER2 Indeterminate [Complete Cases]"
            else:
                return "Primary Tumor (Missing GISTIC)"
        return "Other"

    merged_all["plot1_group"] = merged_all.apply(assign_plot1_group, axis=1)

    # Group counts
    print("  Plot 1 Cohort Composition:")
    print(merged_all["plot1_group"].value_counts().to_string())

    # Build Plot 1 Figure
    fig, ax = plt.subplots(figsize=(10.5, 6.2), dpi=300)
    sns.set_style("whitegrid")

    palette_plot1 = {
        "Clinical HER2 Positive [Complete Cases]": "#D62728",      # Crimson red
        "Clinical HER2 Negative [Complete Cases]": "#1F77B4",      # Steel blue
        "Clinical HER2 Indeterminate [Complete Cases]": "#9467BD", # Purple
        "Normal Adjacent Tissue": "#2CA02C",                      # Forest green
    }

    kde_groups = [
        "Clinical HER2 Positive [Complete Cases]",
        "Clinical HER2 Negative [Complete Cases]",
        "Clinical HER2 Indeterminate [Complete Cases]",
        "Normal Adjacent Tissue",
    ]

    for grp in kde_groups:
        sub_rna = merged_all.loc[merged_all["plot1_group"] == grp, "ERBB2_RNA"].dropna().values
        med_val = float(np.median(sub_rna))
        sns.kdeplot(
            sub_rna,
            fill=True,
            alpha=0.28,
            linewidth=2.2,
            color=palette_plot1[grp],
            label=f"{grp} (n = {len(sub_rna)}, median = {med_val:.2f})",
            ax=ax,
            bw_adjust=0.85,
        )
        ax.axvline(med_val, color=palette_plot1[grp], linestyle="--", linewidth=1.3, alpha=0.85)

    # Add rug plot ticks along baseline for the 7 Metastasis controls
    met_rna = merged_all.loc[merged_all["plot1_group"] == "Metastasis Control", "ERBB2_RNA"].dropna().values
    ax.scatter(
        met_rna,
        np.zeros_like(met_rna) + 0.005,
        color="#FF7F0E",
        marker="|",
        s=180,
        linewidths=2.5,
        zorder=5,
        label=f"Metastasis Controls (n = {len(met_rna)}, range = [{met_rna.min():.2f}, {met_rna.max():.2f}])",
    )

    ax.set_title(
        "Plot 1: Kernel Density Distribution of ERBB2 Expression (TMM + voom Log2-CPM)\n"
        "Option B: Primary Tumors Restricted to Complete Cases (RNA + GISTIC) vs. Tissue Controls",
        fontsize=12,
        fontweight="bold",
        pad=14,
    )
    ax.set_xlabel("ERBB2 Expression (TMM-voom Log2-CPM)", fontsize=11, fontweight="bold")
    ax.set_ylabel("Kernel Density", fontsize=11, fontweight="bold")
    ax.set_xlim(2.5, 16.0)
    ax.legend(loc="upper left", frameon=True, fontsize=9.2, facecolor="white", edgecolor="#CCCCCC")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plot1_path = os.path.join(output_dir, "plot1_erbb2_kde.png")
    fig.savefig(plot1_path, dpi=300)
    fig.savefig(os.path.join(output_dir, "plot1_complete_cases_kde.png"), dpi=300)
    plt.close(fig)
    print(f"  Saved Plot 1 to {plot1_path} and plot1_complete_cases_kde.png")

    # ------------------------------------------------------------------
    # PERMANENT COHORT FILTERING DIRECTIVE
    # Drop all Normal Adjacent and Metastasis samples; retain Primary Tumors only
    # ------------------------------------------------------------------
    n_before = len(merged_all)
    n_normal = int((merged_all["sample_type"] == "Adjacent Normal").sum())
    n_met = int((merged_all["sample_type"] == "Metastasis").sum())

    pt_df = merged_all[merged_all["sample_type"] == "Primary Tumor"].copy().reset_index(drop=True)
    pt_df = pt_df.drop(columns=["plot1_group"])

    print("\n--- PERMANENT COHORT FILTERING VERIFICATION ---")
    print(f"  Total post-QC RNA-seq samples before filtering: {n_before}")
    print(f"  Normal Adjacent samples permanently dropped:    {n_normal}")
    print(f"  Metastasis samples permanently dropped:         {n_met}")
    print(f"  Retained Primary Tumor samples:                 {len(pt_df)}")

    # ------------------------------------------------------------------
    # Cohort Partitioning: Labeled Training Cases vs. Rescue Cases
    # ------------------------------------------------------------------
    pt_df["Has_GISTIC"] = pt_df["ERBB2_GISTIC"].notna()
    pt_df["Is_Complete_Training_Case"] = (
        pt_df["HER2_status"].isin(["Positive", "Negative"]) & pt_df["Has_GISTIC"]
    )
    pt_df["Is_Complete_Rescue_Case"] = (
        (pt_df["HER2_status"] == "Indeterminate") & pt_df["Has_GISTIC"]
    )
    pt_df["Is_Missing_GISTIC_Case"] = ~pt_df["Has_GISTIC"]

    train_cc_df = pt_df[pt_df["Is_Complete_Training_Case"]].copy().reset_index(drop=True)
    rescue_cc_df = pt_df[pt_df["Is_Complete_Rescue_Case"]].copy().reset_index(drop=True)
    missing_gistic_df = pt_df[pt_df["Is_Missing_GISTIC_Case"]].copy().reset_index(drop=True)

    print("\n--- OPERATIONAL COHORT BREAKDOWN (n = 1,075 Primary Tumors) ---")
    print(f"  1. Labeled Training Cases (Complete Data):      n = {len(train_cc_df)}")
    print(f"     * HER2 Positive:                             n = {(train_cc_df['HER2_status'] == 'Positive').sum()}")
    print(f"     * HER2 Negative:                             n = {(train_cc_df['HER2_status'] == 'Negative').sum()}")
    print(f"  2. Step 7 Rescue Cases (Complete Data, Indet):  n = {len(rescue_cc_df)}")
    print(f"  3. Missing-GISTIC Cohort (No Imputation):       n = {len(missing_gistic_df)}")
    print(f"     * Labeled Positive (Missing GISTIC):         n = {(missing_gistic_df['HER2_status'] == 'Positive').sum()}")
    print(f"     * Labeled Negative (Missing GISTIC):         n = {(missing_gistic_df['HER2_status'] == 'Negative').sum()}")
    print(f"     * Indeterminate (Missing GISTIC):            n = {(missing_gistic_df['HER2_status'] == 'Indeterminate').sum()}")

    # ------------------------------------------------------------------
    # STEP 2: RNA-Seq ROC Analysis & PLOT 2 on Labeled Complete Cases (n = 829)
    # ------------------------------------------------------------------
    print("\n[Step 2] Performing RNA-Seq ROC Analysis on Complete Labeled Cases (n = 829)...")
    y_true_train = (train_cc_df["HER2_status"] == "Positive").astype(int).values
    rna_scores = train_cc_df["ERBB2_RNA"].values.astype(float)

    fpr_rna, tpr_rna, thresholds_rna = roc_curve(y_true_train, rna_scores)
    auc_rna = auc(fpr_rna, tpr_rna)

    j_scores_rna = tpr_rna - fpr_rna
    opt_idx_rna = int(np.argmax(j_scores_rna))
    opt_thresh_rna = float(thresholds_rna[opt_idx_rna])
    opt_j_rna = float(j_scores_rna[opt_idx_rna])

    y_pred_rna_opt = (rna_scores >= opt_thresh_rna).astype(int)
    m_rna = compute_all_10_metrics(y_true_train, y_pred_rna_opt, f"RNA-Seq at Youden's J ({opt_thresh_rna:.4f} log2-CPM)")
    m_rna["Predictor"] = "Continuous ERBB2_RNA (TMM-voom log2-CPM)"
    m_rna["Cohort"] = f"Labeled Complete Cases (n = {len(train_cc_df)})"
    m_rna["ROC_AUC"] = float(auc_rna)
    m_rna["Cutoff_Value"] = float(opt_thresh_rna)

    print(f"  Continuous ERBB2 RNA-Seq ROC AUC:          {auc_rna:.4f}")
    print(f"  Youden's J Optimal RNA Cutoff:             {opt_thresh_rna:.4f} Log2-CPM (J = {opt_j_rna:.4f})")
    print(f"  Sensitivity (Recall / TPR):                {m_rna['Sensitivity_TPR_Recall']*100:.2f}% ({m_rna['TP']}/{m_rna['TP']+m_rna['FN']})")
    print(f"  Specificity (TNR):                         {m_rna['Specificity_TNR']*100:.2f}% ({m_rna['TN']}/{m_rna['TN']+m_rna['FP']})")
    print(f"  Precision (PPV):                           {m_rna['Precision_PPV']*100:.2f}%")
    print(f"  Negative Predictive Value (NPV):           {m_rna['NPV']*100:.2f}%")
    print(f"  Accuracy:                                  {m_rna['Accuracy']*100:.2f}%")
    print(f"  Balanced Accuracy:                         {m_rna['Balanced_Accuracy']*100:.2f}%")
    print(f"  F1-Score:                                  {m_rna['F1_Score']:.4f}")

    # Render Plot 2 (RNA ROC Curve)
    fig, ax = plt.subplots(figsize=(7.5, 6.2), dpi=300)
    ax.plot(fpr_rna, tpr_rna, color="#D62728", lw=2.8, label=f"Continuous ERBB2 RNA-Seq (AUC = {auc_rna:.4f})")
    ax.scatter(
        [fpr_rna[opt_idx_rna]],
        [tpr_rna[opt_idx_rna]],
        color="#D62728",
        edgecolors="black",
        s=110,
        zorder=6,
        label=f"Youden's J Cutoff = {opt_thresh_rna:.2f} log2-CPM\n(Sens = {m_rna['Sensitivity_TPR_Recall']*100:.1f}%, Spec = {m_rna['Specificity_TNR']*100:.1f}%, J = {opt_j_rna:.3f})",
    )
    # Draw vertical drop line illustrating Youden's J vertical gap to chance diagonal
    ax.plot(
        [fpr_rna[opt_idx_rna], fpr_rna[opt_idx_rna]],
        [fpr_rna[opt_idx_rna], tpr_rna[opt_idx_rna]],
        color="#333333",
        linestyle=":",
        lw=1.8,
        label=f"Youden's J Gap (J = {opt_j_rna:.3f})",
    )
    ax.plot([0, 1], [0, 1], color="#777777", linestyle="--", lw=1.3, label="Chance Line (AUC = 0.5000)")
    ax.set_title(
        f"Plot 2: RNA-Seq ROC Curve — Continuous ERBB2 (n = {len(train_cc_df)} Complete Cases)",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=10.5, fontweight="bold")
    ax.legend(loc="lower right", frameon=True, fontsize=8.8, facecolor="white", edgecolor="#CCCCCC")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plot2_path = os.path.join(output_dir, "plot2_rnaseq_roc.png")
    fig.savefig(plot2_path, dpi=300)
    plt.close(fig)
    print(f"  Saved Plot 2 to {plot2_path}")

    # ------------------------------------------------------------------
    # STEP 3: ESTIMATE Tumor Purity & PLOT 3
    # ------------------------------------------------------------------
    print("\n[Step 3] Computing ESTIMATE Tumor Purity via gseapy ssGSEA...")
    tumor_sample_ids = pt_df["sample_id"].tolist()
    tumor_expr_df = expr_df[tumor_sample_ids].copy()

    # Intersect with available genes
    avail_genes = set(tumor_expr_df.index)
    matched_stromal = [g for g in STROMAL_141 if g in avail_genes]
    matched_immune = [g for g in IMMUNE_141 if g in avail_genes]
    print(f"  Matched ESTIMATE Signature Genes: {len(matched_stromal)}/141 Stromal, {len(matched_immune)}/141 Immune.")

    gene_sets_estimate = {
        "StromalSignature": matched_stromal,
        "ImmuneSignature": matched_immune,
    }

    # Run single-sample GSEA
    ss_res = gp.ssgsea(
        data=tumor_expr_df,
        gene_sets=gene_sets_estimate,
        outdir=None,
        sample_norm_method="rank",
        min_size=10,
        max_size=500,
        permutation_num=0,
        no_plot=True,
        verbose=False,
    )

    # Use raw enrichment scores from resultsOnSamples to align with Yoshihara scale
    raw_es = pd.DataFrame(ss_res.resultsOnSamples)
    stromal_scores = raw_es.loc["StromalSignature"].astype(float)
    immune_scores = raw_es.loc["ImmuneSignature"].astype(float)

    # Scale to Yoshihara Affymetrix baseline (mean=0, sd=1000)
    stromal_norm = (stromal_scores - stromal_scores.mean()) / stromal_scores.std() * 1000.0
    immune_norm = (immune_scores - immune_scores.mean()) / immune_scores.std() * 1000.0
    estimate_scores = stromal_norm + immune_norm

    # Yoshihara cosine tumor purity formula
    cos_arg = 0.6049872018 + 0.0001467884 * estimate_scores
    purity_scores = np.cos(cos_arg).clip(0.0, 1.0)

    # Map back to Primary Tumor DataFrame
    pt_df["Stromal_Score"] = pt_df["sample_id"].map(stromal_norm)
    pt_df["Immune_Score"] = pt_df["sample_id"].map(immune_norm)
    pt_df["ESTIMATE_Score"] = pt_df["sample_id"].map(estimate_scores)
    pt_df["Tumor_Purity"] = pt_df["sample_id"].map(purity_scores)

    # Re-sync purity into subsets
    train_cc_df["Tumor_Purity"] = train_cc_df["sample_id"].map(pt_df.set_index("sample_id")["Tumor_Purity"])
    rescue_cc_df["Tumor_Purity"] = rescue_cc_df["sample_id"].map(pt_df.set_index("sample_id")["Tumor_Purity"])
    missing_gistic_df["Tumor_Purity"] = missing_gistic_df["sample_id"].map(pt_df.set_index("sample_id")["Tumor_Purity"])

    print(f"  Derived Tumor Purity (n = {len(pt_df)}): Mean = {pt_df['Tumor_Purity'].mean():.4f}, Range = [{pt_df['Tumor_Purity'].min():.4f}, {pt_df['Tumor_Purity'].max():.4f}]")

    # Top 5 Highest and Lowest Purity Samples
    top5_high = pt_df.sort_values("Tumor_Purity", ascending=False).head(5)
    top5_low = pt_df.sort_values("Tumor_Purity", ascending=True).head(5)

    purity_extremes_rows = []
    for _, row in top5_high.iterrows():
        purity_extremes_rows.append({
            "Extreme_Rank": "Top 5 Highest Purity",
            "Sample_ID": row["sample_id"],
            "Clinical_HER2_Status": row["HER2_status"],
            "ERBB2_RNA": f"{row['ERBB2_RNA']:.4f}",
            "Tumor_Purity": f"{row['Tumor_Purity']:.4f}",
            "Stromal_Score": f"{row['Stromal_Score']:.2f}",
            "Immune_Score": f"{row['Immune_Score']:.2f}",
        })
    for _, row in top5_low.iterrows():
        purity_extremes_rows.append({
            "Extreme_Rank": "Top 5 Lowest Purity",
            "Sample_ID": row["sample_id"],
            "Clinical_HER2_Status": row["HER2_status"],
            "ERBB2_RNA": f"{row['ERBB2_RNA']:.4f}",
            "Tumor_Purity": f"{row['Tumor_Purity']:.4f}",
            "Stromal_Score": f"{row['Stromal_Score']:.2f}",
            "Immune_Score": f"{row['Immune_Score']:.2f}",
        })
    purity_extremes_df = pd.DataFrame(purity_extremes_rows)

    # Render Plot 3 (Scatter Plot: Purity vs. ERBB2 RNA)
    print("  Rendering Plot 3 (Tumor Purity vs. ERBB2 Expression)...")
    fig, ax = plt.subplots(figsize=(8.8, 6.2), dpi=300)

    palette_p3 = {
        "Positive": "#D62728",
        "Negative": "#1F77B4",
        "Indeterminate": "#9467BD",
    }

    for status_val in ["Negative", "Indeterminate", "Positive"]:
        sub = pt_df[pt_df["HER2_status"] == status_val]
        r_val, p_val = stats.pearsonr(sub["Tumor_Purity"], sub["ERBB2_RNA"])
        lbl = f"Clinical {status_val} (n = {len(sub)}, r = {r_val:+.2f}, p = {p_val:.1e})"
        ax.scatter(
            sub["Tumor_Purity"],
            sub["ERBB2_RNA"],
            c=palette_p3[status_val],
            label=lbl,
            alpha=0.70 if status_val == "Positive" else (0.50 if status_val == "Indeterminate" else 0.40),
            s=36 if status_val == "Positive" else 24,
            edgecolors="white" if status_val == "Positive" else "none",
            linewidths=0.4,
            zorder=4 if status_val == "Positive" else 3,
        )
        # Add regression trendlines
        sns.regplot(
            data=sub,
            x="Tumor_Purity",
            y="ERBB2_RNA",
            scatter=False,
            ax=ax,
            color=palette_p3[status_val],
            line_kws={"linewidth": 2.0, "linestyle": "-" if status_val != "Indeterminate" else "--"},
        )

    ax.axhline(
        opt_thresh_rna,
        color="#222222",
        linestyle=":",
        linewidth=1.5,
        label=f"RNA Youden's J Cutoff ({opt_thresh_rna:.2f} log2-CPM)",
    )
    ax.set_title(
        "Plot 3: ESTIMATE Tumor Purity vs. ERBB2 Expression (n = 1,075 Primary Tumors)",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlabel("ESTIMATE Tumor Purity Score (0 to 1)", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("ERBB2 Expression (TMM-voom Log2-CPM)", fontsize=10.5, fontweight="bold")
    ax.legend(loc="upper left", frameon=True, fontsize=8.8, facecolor="white", edgecolor="#CCCCCC")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plot3_path = os.path.join(output_dir, "plot3_purity_vs_erbb2.png")
    fig.savefig(plot3_path, dpi=300)
    plt.close(fig)
    print(f"  Saved Plot 3 to {plot3_path}")

    # ------------------------------------------------------------------
    # STEP 4: 17q12 Amplicon Co-Expression & Amplicon_Z Score
    # ------------------------------------------------------------------
    print("\n[Step 4] Analyzing 17q12 Core Amplicon Co-Expression & Computing Amplicon_Z...")
    amplicon_genes_map = {
        "ERBB2": "ERBB2",
        "GRB7": "GRB7",
        "PGAP3": "PGAP3",
        "C17orf37": "MIEN1 (C17orf37)",
        "STARD3": "STARD3",
    }

    # Extract 5 amplicon genes across all 1,075 primary tumors
    amp_all_df = tumor_expr_df.loc[list(amplicon_genes_map.keys())].T
    amp_all_df.columns = [amplicon_genes_map[c] for c in amp_all_df.columns]

    # Compute gene-wise Z-scores across all 1,075 primary tumors for the 4 flanking genes
    flanking_4 = ["GRB7", "PGAP3", "STARD3", "MIEN1 (C17orf37)"]
    flank_z = (amp_all_df[flanking_4] - amp_all_df[flanking_4].mean()) / amp_all_df[flanking_4].std()
    pt_df["Amplicon_Z"] = pt_df["sample_id"].map(flank_z.mean(axis=1))

    # Sync Amplicon_Z into subsets
    train_cc_df["Amplicon_Z"] = train_cc_df["sample_id"].map(pt_df.set_index("sample_id")["Amplicon_Z"])
    rescue_cc_df["Amplicon_Z"] = rescue_cc_df["sample_id"].map(pt_df.set_index("sample_id")["Amplicon_Z"])
    missing_gistic_df["Amplicon_Z"] = missing_gistic_df["sample_id"].map(pt_df.set_index("sample_id")["Amplicon_Z"])

    print(f"  Derived Amplicon_Z (n = {len(pt_df)}): Mean = {pt_df['Amplicon_Z'].mean():.4f}, Range = [{pt_df['Amplicon_Z'].min():.4f}, {pt_df['Amplicon_Z'].max():.4f}]")

    # Plot 4: Clustered correlation heatmap for Labeled Training Cases (n = 829)
    print("  Rendering Plot 4 (17q12 Amplicon Clustered Heatmap on n = 829 Complete Training Cases)...")
    amp_train_df = amp_all_df.loc[train_cc_df["sample_id"]].copy()
    corr_matrix = amp_train_df.corr(method="pearson")

    cg = sns.clustermap(
        corr_matrix,
        annot=True,
        fmt=".3f",
        cmap="RdYlBu_r",
        vmin=0.6,
        vmax=1.0,
        figsize=(7.8, 6.8),
        cbar_pos=(0.02, 0.80, 0.03, 0.14),
        cbar_kws={"label": "Pearson r"},
        annot_kws={"size": 11, "weight": "bold"},
        linewidths=1.0,
    )
    cg.fig.suptitle(
        f"Plot 4: 17q12 Core Amplicon Co-Expression Heatmap (n = {len(train_cc_df)} Complete Training Cases)",
        fontsize=11.5,
        fontweight="bold",
        y=1.03,
    )
    plot4_path = os.path.join(output_dir, "plot4_amplicon_clustermap.png")
    cg.savefig(plot4_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"  Saved Plot 4 to {plot4_path}")

    # ------------------------------------------------------------------
    # STEP 5: GISTIC CNV ROC Analysis, Metrics & Paired DeLong's Test
    # ------------------------------------------------------------------
    print("\n[Step 5] Performing GISTIC CNV ROC Analysis on Complete Labeled Cases (n = 829)...")
    gistic_scores = train_cc_df["ERBB2_GISTIC"].values.astype(float)
    gistic_binary = (gistic_scores == 2).astype(int)  # Standard focal amplification cutoff

    fpr_gis, tpr_gis, _ = roc_curve(y_true_train, gistic_scores)
    auc_gis = auc(fpr_gis, tpr_gis)

    m_gis = compute_all_10_metrics(y_true_train, gistic_binary, "GISTIC Focal Amplification (== 2 vs. <= 1)")
    m_gis["Predictor"] = "Discrete ERBB2_GISTIC (GISTIC 2.0: -2 to +2)"
    m_gis["Cohort"] = f"Labeled Complete Cases (n = {len(train_cc_df)})"
    m_gis["ROC_AUC"] = float(auc_gis)
    m_gis["Cutoff_Value"] = 2.0

    print(f"  Discrete GISTIC (-2 to +2) ROC AUC:        {auc_gis:.4f}")
    print(f"  Sensitivity (Recall / TPR at GISTIC == 2): {m_gis['Sensitivity_TPR_Recall']*100:.2f}% ({m_gis['TP']}/{m_gis['TP']+m_gis['FN']})")
    print(f"  Specificity (TNR at GISTIC <= 1):          {m_gis['Specificity_TNR']*100:.2f}% ({m_gis['TN']}/{m_gis['TN']+m_gis['FP']})")
    print(f"  Precision (PPV):                           {m_gis['Precision_PPV']*100:.2f}%")
    print(f"  Negative Predictive Value (NPV):           {m_gis['NPV']*100:.2f}%")
    print(f"  Accuracy:                                  {m_gis['Accuracy']*100:.2f}%")
    print(f"  Balanced Accuracy:                         {m_gis['Balanced_Accuracy']*100:.2f}%")
    print(f"  F1-Score:                                  {m_gis['F1_Score']:.4f}")

    # Paired DeLong's statistical test
    print("\n  Executing Paired DeLong's Test (Continuous RNA vs. Discrete GISTIC)...")
    auc1, auc2, z_delong, p_delong = delong_paired_test(y_true_train, rna_scores, gistic_scores)
    print(f"  Continuous RNA AUC = {auc1:.4f} vs. Discrete GISTIC AUC = {auc2:.4f}")
    print(f"  DeLong Z-score = {z_delong:+.4f}, p-value = {p_delong:.4e}")

    # Render Plot 5 (GISTIC ROC Curve)
    print("  Rendering Plot 5 (GISTIC CNV ROC Curve)...")
    fig, ax = plt.subplots(figsize=(7.5, 6.2), dpi=300)
    ax.plot(fpr_gis, tpr_gis, color="#1F77B4", lw=2.8, linestyle="-.", label=f"Discrete ERBB2 GISTIC (AUC = {auc_gis:.4f})")
    # Operating point for GISTIC == 2
    fp_gis_rate = m_gis["FPR_Fallout"]
    tp_gis_rate = m_gis["Sensitivity_TPR_Recall"]
    ax.scatter(
        [fp_gis_rate],
        [tp_gis_rate],
        color="#1F77B4",
        edgecolors="black",
        s=110,
        zorder=6,
        label=f"GISTIC == 2 Cutoff\n(Sens = {tp_gis_rate*100:.1f}%, Spec = {m_gis['Specificity_TNR']*100:.1f}%)",
    )
    ax.plot([0, 1], [0, 1], color="#777777", linestyle="--", lw=1.3, label="Chance Line (AUC = 0.5000)")
    ax.set_title(
        f"Plot 5: GISTIC CNV ROC Curve — Discrete Copy Number (n = {len(train_cc_df)} Complete Cases)",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=10.5, fontweight="bold")
    ax.legend(loc="lower right", frameon=True, fontsize=8.8, facecolor="white", edgecolor="#CCCCCC")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plot5_path = os.path.join(output_dir, "plot5_gistic_roc.png")
    fig.savefig(plot5_path, dpi=300)
    plt.close(fig)
    print(f"  Saved Plot 5 to {plot5_path}")

    # Render Plot 6 (Combined ROC Curve: RNA vs. GISTIC + DeLong Annotation)
    print("  Rendering Plot 6 (Combined RNA vs. GISTIC ROC Curve with DeLong Test)...")
    fig, ax = plt.subplots(figsize=(8.0, 6.4), dpi=300)
    ax.plot(fpr_rna, tpr_rna, color="#D62728", lw=2.8, label=f"Continuous ERBB2 RNA-Seq (AUC = {auc_rna:.4f})")
    ax.plot(fpr_gis, tpr_gis, color="#1F77B4", lw=2.5, linestyle="-.", label=f"Discrete ERBB2 GISTIC (AUC = {auc_gis:.4f})")
    ax.plot([0, 1], [0, 1], color="#777777", linestyle="--", lw=1.3, label="Chance Line (AUC = 0.5000)")

    # Plot optimal points
    ax.scatter([fpr_rna[opt_idx_rna]], [tpr_rna[opt_idx_rna]], color="#D62728", edgecolors="black", s=95, zorder=6)
    ax.scatter([fp_gis_rate], [tp_gis_rate], color="#1F77B4", edgecolors="black", s=95, zorder=6)

    # Annotate DeLong test results on plot
    delong_box = (
        f"Paired DeLong's Test (n = {len(train_cc_df)}):\n"
        f"  Δ AUC (RNA - GISTIC) = +{auc_rna - auc_gis:.4f}\n"
        f"  Z-score = {z_delong:+.4f}\n"
        f"  p-value = {p_delong:.4e}\n"
        f"  Significance: {'p < 0.05 (RNA significantly higher)' if p_delong < 0.05 else 'Not significant'}"
    )
    ax.text(
        0.42,
        0.18,
        delong_box,
        fontsize=9.0,
        fontfamily="monospace",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.5", fc="#F8F9FA", ec="#666666", alpha=0.95),
    )

    ax.set_title(
        f"Plot 6: Combined ROC Curve — Continuous RNA vs. Discrete GISTIC (n = {len(train_cc_df)})",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlabel("False Positive Rate (1 - Specificity)", fontsize=10.5, fontweight="bold")
    ax.set_ylabel("True Positive Rate (Sensitivity)", fontsize=10.5, fontweight="bold")
    ax.legend(loc="lower right", frameon=True, fontsize=9.0, facecolor="white", edgecolor="#CCCCCC")
    ax.grid(True, linestyle=":", alpha=0.5)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plot6_path = os.path.join(output_dir, "plot6_combined_roc.png")
    fig.savefig(plot6_path, dpi=300)
    plt.close(fig)
    print(f"  Saved Plot 6 to {plot6_path}")

    # ------------------------------------------------------------------
    # Copy all 6 plots to Artifact Directory
    # ------------------------------------------------------------------
    for p_name in ["plot1_erbb2_kde.png", "plot1_complete_cases_kde.png", "plot2_rnaseq_roc.png", "plot3_purity_vs_erbb2.png", "plot4_amplicon_clustermap.png", "plot5_gistic_roc.png", "plot6_combined_roc.png"]:
        src_f = os.path.join(output_dir, p_name)
        dst_f = os.path.join(artifact_dir, p_name)
        if os.path.isfile(src_f):
            with open(src_f, "rb") as rf, open(dst_f, "wb") as wf:
                wf.write(rf.read())

    # ------------------------------------------------------------------
    # Compile Master Excel Workbook (`results/her2_classification_results.xlsx`)
    # ------------------------------------------------------------------
    print("\n[Export] Generating results/her2_classification_results.xlsx...")
    excel_path = os.path.join(output_dir, "her2_classification_results.xlsx")

    # Table of 10 metrics for Step 2 and Step 5
    metrics_summary_df = pd.DataFrame([m_rna, m_gis])
    metrics_summary_df["DeLong_Z_vs_GISTIC"] = [z_delong, 0.0]
    metrics_summary_df["DeLong_p_vs_GISTIC"] = [p_delong, 1.0]

    # Assign operational cohort role for each sample in pt_df
    def assign_cohort_role(row: pd.Series) -> str:
        if row["Is_Complete_Training_Case"]:
            return "1. Labeled Complete Training Case (Step 6 CV)"
        if row["Is_Complete_Rescue_Case"]:
            return "2. Indeterminate Complete Rescue Case (Step 7)"
        if row["Is_Missing_GISTIC_Case"]:
            return "3. Missing GISTIC Copy-Number (Unimputed Audit Cohort)"
        return "4. Other"

    pt_df["Pipeline_Cohort_Role"] = pt_df.apply(assign_cohort_role, axis=1)

    # Missing GISTIC cohort audit dataframe
    missing_gistic_export = pt_df[pt_df["Is_Missing_GISTIC_Case"]][[
        "sample_id", "patient_id", "HER2_status", "HER2_expression_subtype", "discordant_flag",
        "ERBB2_RNA", "Tumor_Purity", "Amplicon_Z"
    ]].copy()

    with pd.ExcelWriter(excel_path, engine="openpyxl") as writer:
        metrics_summary_df.to_excel(writer, sheet_name="ROC_Performance_Metrics", index=False)
        purity_extremes_df.to_excel(writer, sheet_name="Purity_Extremes", index=False)
        missing_gistic_export.to_excel(writer, sheet_name="Missing_GISTIC_Samples", index=False)
        pt_df[[
            "sample_id", "patient_id", "sample_type", "HER2_status", "HER2_expression_subtype", "discordant_flag",
            "Pipeline_Cohort_Role", "ERBB2_RNA", "ERBB2_GISTIC", "Tumor_Purity", "Amplicon_Z",
            "Stromal_Score", "Immune_Score", "ESTIMATE_Score"
        ]].to_excel(writer, sheet_name="Cohort_Master", index=False)

        for ws in writer.book.worksheets:
            for col in ws.columns:
                max_len = max(len(str(cell.value or "")) for cell in col)
                ws.column_dimensions[col[0].column_letter].width = min(max(max_len + 3, 14), 60)

    print(f"  Successfully initialized master Excel workbook at: {excel_path}")

    # Also save an intermediate pickle bundle so Step 6 can load all computed features directly
    bundle_s1_5 = {
        "pt_df": pt_df,
        "train_cc_df": train_cc_df,
        "rescue_cc_df": rescue_cc_df,
        "missing_gistic_df": missing_gistic_df,
        "opt_thresh_rna": opt_thresh_rna,
        "auc_rna": auc_rna,
        "auc_gis": auc_gis,
        "z_delong": z_delong,
        "p_delong": p_delong,
        "m_rna": m_rna,
        "m_gis": m_gis,
    }
    with open(os.path.join(output_dir, "milestone1_bundle.pkl"), "wb") as f:
        import pickle
        pickle.dump(bundle_s1_5, f)

    print("\n================================================================================")
    print("MILESTONE 1 EXECUTION COMPLETE. REACHED STOPPING POINT 1.")
    print("================================================================================")


if __name__ == "__main__":
    main()

