#!/usr/bin/env python3
"""
HER2 Purity-Aware Classification Pipeline - Step 7 (Milestone 3)

Covers:
  - Step 7.1: Score the 123 Indeterminate Primary Tumors with complete data.
  - Step 7.2: Recovery Tiering:
      - Tier 1: High-Confidence Rescue (RF Prob >= 0.70)
      - Tier 2: Moderate-Confidence Rescue (Optimal Threshold 0.5880 <= RF Prob < 0.70)
      - Tier 3: Negative / Discard (RF Prob < 0.5880)
      - Reporting Logistic Regression Probabilities & Tiers side-by-side.
  - Step 7.3: Render Plot 8 (Multi-panel 2D Facet Plot colored by Tier Assignment).
  - Step 7.4: Bonus Discovery Step: Classify 7 excluded Metastasis control samples.
  - Step 7.5: Finalize master Excel workbook (results/her2_classification_results.xlsx).
"""

import os
import sys
import pickle
import warnings

sys.path.insert(0, os.path.abspath("src"))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_cache")
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


def main() -> None:
    print("================================================================================")
    print("STARTING STEP 7: RESCUE TIERING, PLOT 8 & METASTASIS INSPECTION")
    print("================================================================================")

    output_dir = "results"
    artifact_dir = "/Users/galebabu-ali/.gemini/antigravity/brain/0ae2581d-851b-44c2-9836-f4a6ab72a547"

    # ------------------------------------------------------------------
    # Load Bundles & Models
    # ------------------------------------------------------------------
    print("\n[Step 7.1] Loading Precomputed Features & Production Models...")
    bundle_path = os.path.join(output_dir, "milestone1_bundle.pkl")
    prod_path = os.path.join(output_dir, "production_models.pkl")

    with open(bundle_path, "rb") as f:
        m1_bundle = pickle.load(f)
    with open(prod_path, "rb") as f:
        prod_bundle = pickle.load(f)

    rescue_df = m1_bundle["rescue_cc_df"].copy().reset_index(drop=True)
    rf_model = prod_bundle["rf_model"]
    lr_model = prod_bundle["lr_model"]
    scaler = prod_bundle["scaler"]
    features = prod_bundle["features"]
    opt_thresh_rf = float(prod_bundle["opt_thresh_rf"])
    opt_thresh_lr = float(prod_bundle["opt_thresh_lr"])
    opt_thresh_rna = float(m1_bundle["opt_thresh_rna"])

    print(f"  Rescue Cohort: n = {len(rescue_df)} Indeterminate Primary Tumors with complete data.")
    print(f"  Features: {features}")
    print(f"  Optimal RF Probability Threshold: {opt_thresh_rf:.4f}")
    print(f"  Optimal LR Probability Threshold: {opt_thresh_lr:.4f}")

    # ------------------------------------------------------------------
    # Step 7.2: Predict Probabilities & Assign Confidence Tiers
    # ------------------------------------------------------------------
    print("\n[Step 7.2] Generating Predictions and Assigning Confidence Tiers...")
    X_rescue = rescue_df[features].values
    X_rescue_scaled = scaler.transform(X_rescue)

    rf_probs = rf_model.predict_proba(X_rescue)[:, 1]
    lr_probs = lr_model.predict_proba(X_rescue_scaled)[:, 1]

    rescue_df["RF_Probability"] = np.round(rf_probs, 4)
    rescue_df["LR_Probability"] = np.round(lr_probs, 4)

    def assign_rf_tier(p: float) -> str:
        if p >= 0.70:
            return "Tier 1: High-Confidence Rescue"
        elif p >= opt_thresh_rf:
            return "Tier 2: Moderate-Confidence Rescue"
        else:
            return "Tier 3: Negative / Discard"

    def assign_lr_tier(p: float) -> str:
        if p >= 0.70:
            return "Tier 1: High-Confidence Rescue"
        elif p >= opt_thresh_lr:
            return "Tier 2: Moderate-Confidence Rescue"
        else:
            return "Tier 3: Negative / Discard"

    rescue_df["Rescue_Tier_RF"] = rescue_df["RF_Probability"].apply(assign_rf_tier)
    rescue_df["Rescue_Tier_LR"] = rescue_df["LR_Probability"].apply(assign_lr_tier)

    print("\n  Random Forest Tier Distribution (Primary Classifier):")
    rf_counts = rescue_df["Rescue_Tier_RF"].value_counts()
    for tier, count in rf_counts.items():
        print(f"    * {tier:<36}: n = {count:3d} ({count / len(rescue_df) * 100:.1f}%)")

    print("\n  Logistic Regression Tier Distribution (Side-by-Side):")
    lr_counts = rescue_df["Rescue_Tier_LR"].value_counts()
    for tier, count in lr_counts.items():
        print(f"    * {tier:<36}: n = {count:3d} ({count / len(rescue_df) * 100:.1f}%)")

    print("\n  Tier Concordance Matrix (RF vs. LR):")
    concordance_df = pd.crosstab(rescue_df["Rescue_Tier_RF"], rescue_df["Rescue_Tier_LR"])
    print(concordance_df.to_string())

    # ------------------------------------------------------------------
    # Step 7.3: PLOT 8 - Multi-Panel 2D Facet Plot
    # ------------------------------------------------------------------
    print("\n[Step 7.3] Rendering Plot 8 (Multi-Panel Facet Plot of Rescued Cases)...")
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.5), dpi=300)
    sns.set_style("whitegrid")

    palette_tier = {
        "Tier 1: High-Confidence Rescue": "#D62728",     # Crimson Red
        "Tier 2: Moderate-Confidence Rescue": "#FF7F0E", # Vivid Orange
        "Tier 3: Negative / Discard": "#1F77B4",         # Steel Blue
    }

    tier_order = [
        "Tier 3: Negative / Discard",
        "Tier 2: Moderate-Confidence Rescue",
        "Tier 1: High-Confidence Rescue",
    ]

    tier_counts = rescue_df["Rescue_Tier_RF"].value_counts().to_dict()

    # Panel A: Tumor Purity vs. ERBB2 RNA
    ax = axes[0]
    for tier in tier_order:
        sub = rescue_df[rescue_df["Rescue_Tier_RF"] == tier]
        if len(sub) == 0:
            continue
        n_t = tier_counts.get(tier, len(sub))
        ax.scatter(
            sub["Tumor_Purity"],
            sub["ERBB2_RNA"],
            c=palette_tier[tier],
            label=f"{tier} (n = {n_t})",
            s=80 if "Tier 1" in tier else (65 if "Tier 2" in tier else 35),
            marker="*" if "Tier 1" in tier else ("D" if "Tier 2" in tier else "o"),
            alpha=0.90 if "Tier 1" in tier else (0.85 if "Tier 2" in tier else 0.50),
            edgecolors="black" if "Tier 3" not in tier else "none",
            linewidths=0.8,
            zorder=5 if "Tier 1" in tier else (4 if "Tier 2" in tier else 3),
        )
    ax.axhline(opt_thresh_rna, color="#333333", linestyle=":", lw=1.5, label=f"RNA Youden ({opt_thresh_rna:.2f})")
    ax.set_title("Panel A: Tumor Purity vs. ERBB2 RNA", fontsize=11, fontweight="bold")
    ax.set_xlabel("ESTIMATE Tumor Purity Score", fontsize=10, fontweight="bold")
    ax.set_ylabel("ERBB2 Expression (TMM-voom Log2-CPM)", fontsize=10, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8.2, frameon=True, facecolor="white")
    ax.grid(True, linestyle=":", alpha=0.5)

    # Panel B: Amplicon_Z vs. ERBB2 RNA
    ax = axes[1]
    for tier in tier_order:
        sub = rescue_df[rescue_df["Rescue_Tier_RF"] == tier]
        if len(sub) == 0:
            continue
        ax.scatter(
            sub["Amplicon_Z"],
            sub["ERBB2_RNA"],
            c=palette_tier[tier],
            label=f"{tier}",
            s=80 if "Tier 1" in tier else (65 if "Tier 2" in tier else 35),
            marker="*" if "Tier 1" in tier else ("D" if "Tier 2" in tier else "o"),
            alpha=0.90 if "Tier 1" in tier else (0.85 if "Tier 2" in tier else 0.50),
            edgecolors="black" if "Tier 3" not in tier else "none",
            linewidths=0.8,
            zorder=5 if "Tier 1" in tier else (4 if "Tier 2" in tier else 3),
        )
    ax.axhline(opt_thresh_rna, color="#333333", linestyle=":", lw=1.5, label=f"RNA Youden ({opt_thresh_rna:.2f})")
    ax.axvline(0.0, color="#777777", linestyle="--", lw=1.0)
    ax.set_title("Panel B: 17q12 Amplicon_Z vs. ERBB2 RNA", fontsize=11, fontweight="bold")
    ax.set_xlabel("17q12 Amplicon Mean Z-Score (4 Genes)", fontsize=10, fontweight="bold")
    ax.set_ylabel("ERBB2 Expression (TMM-voom Log2-CPM)", fontsize=10, fontweight="bold")
    ax.legend(loc="upper left", fontsize=8.2, frameon=True, facecolor="white")
    ax.grid(True, linestyle=":", alpha=0.5)

    # Panel C: ERBB2 GISTIC vs. ERBB2 RNA
    ax = axes[2]
    for tier in tier_order:
        sub = rescue_df[rescue_df["Rescue_Tier_RF"] == tier]
        if len(sub) == 0:
            continue
        # Add slight jitter to GISTIC on x-axis for visibility
        x_jitter = sub["ERBB2_GISTIC"] + np.random.normal(0, 0.05, size=len(sub))
        ax.scatter(
            x_jitter,
            sub["ERBB2_RNA"],
            c=palette_tier[tier],
            label=f"{tier}",
            s=80 if "Tier 1" in tier else (65 if "Tier 2" in tier else 35),
            marker="*" if "Tier 1" in tier else ("D" if "Tier 2" in tier else "o"),
            alpha=0.90 if "Tier 1" in tier else (0.85 if "Tier 2" in tier else 0.50),
            edgecolors="black" if "Tier 3" not in tier else "none",
            linewidths=0.8,
            zorder=5 if "Tier 1" in tier else (4 if "Tier 2" in tier else 3),
        )
    ax.axhline(opt_thresh_rna, color="#333333", linestyle=":", lw=1.5, label=f"RNA Youden ({opt_thresh_rna:.2f})")
    ax.axvline(1.5, color="#D62728", linestyle="--", lw=1.2, label="Focal Amp Cutoff (GISTIC == 2)")
    ax.set_title("Panel C: ERBB2 GISTIC Copy Number vs. RNA", fontsize=11, fontweight="bold")
    ax.set_xlabel("ERBB2 GISTIC Score (-2 to +2, jittered)", fontsize=10, fontweight="bold")
    ax.set_ylabel("ERBB2 Expression (TMM-voom Log2-CPM)", fontsize=10, fontweight="bold")
    ax.set_xticks([-2, -1, 0, 1, 2])
    ax.legend(loc="upper left", fontsize=8.2, frameon=True, facecolor="white")
    ax.grid(True, linestyle=":", alpha=0.5)

    fig.suptitle(
        f"Plot 8: Rescued Indeterminate Cases (n = {len(rescue_df)}) Stratified by Purity-Aware Confidence Tiers\n"
        f"Tier 1: High-Confidence (RF Prob >= 0.70, n = {tier_counts.get('Tier 1: High-Confidence Rescue', 0)}) | "
        f"Tier 2: Moderate-Confidence ({opt_thresh_rf:.4f} <= Prob < 0.70, n = {tier_counts.get('Tier 2: Moderate-Confidence Rescue', 0)}) | "
        f"Tier 3: Negative / Discard (Prob < {opt_thresh_rf:.4f}, n = {tier_counts.get('Tier 3: Negative / Discard', 0)})",
        fontsize=12,
        fontweight="bold",
        y=1.03,
    )

    plt.tight_layout()
    plot8_path = os.path.join(output_dir, "plot8_rescue_facet_plot.png")
    fig.savefig(plot8_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved Plot 8 to {plot8_path}")

    # Copy Plot 8 to artifact directory
    dst_plot8 = os.path.join(artifact_dir, "plot8_rescue_facet_plot.png")
    with open(plot8_path, "rb") as rf_in, open(dst_plot8, "wb") as wf_out:
        wf_out.write(rf_in.read())
    print(f"  Synced Plot 8 to artifact directory: {dst_plot8}")

    # ------------------------------------------------------------------
    # Step 7.4: Bonus Discovery Step - Metastasis Samples (n = 7)
    # ------------------------------------------------------------------
    print("\n[Step 7.4] Evaluating Excluded Metastasis Samples (n = 7)...")
    meta_path = "results/rsem_sample_metadata.csv"
    expr_path = "results/matrix_normalized_log2_stabilized.csv"
    cn_path = "data/brca_tcga_erbb2_copy_number.csv"
    clin_path = "results/her2_clinical_status.csv"

    meta_df = pd.read_csv(meta_path)
    expr_df = pd.read_csv(expr_path).set_index("gene_symbol")
    cn_df = pd.read_csv(cn_path)
    clin_df = pd.read_csv(clin_path, keep_default_na=False)

    met_df = meta_df[meta_df["sample_type"] == "Metastasis"].copy().reset_index(drop=True)
    met_df["ERBB2_RNA"] = met_df["sample_id"].map(expr_df.loc["ERBB2"]).astype(float)

    # Amplicon Z for metastasis relative to primary tumors
    flanking_4 = ["GRB7", "PGAP3", "STARD3", "C17orf37"]
    pt_ids = m1_bundle["pt_df"]["sample_id"].tolist()
    pt_amp_means = expr_df.loc[flanking_4, pt_ids].mean(axis=1)
    pt_amp_stds = expr_df.loc[flanking_4, pt_ids].std(axis=1)

    met_amp_df = expr_df.loc[flanking_4, met_df["sample_id"]].T
    met_z = (met_amp_df - pt_amp_means) / pt_amp_stds
    met_df["Amplicon_Z"] = met_z.mean(axis=1).values

    # Tumor Purity for metastasis samples
    # Yoshihara cosine ESTIMATE scoring
    if not hasattr(pd.DataFrame, "iteritems"):
        pd.DataFrame.iteritems = pd.DataFrame.items
    from her2_pipeline_milestone1 import IMMUNE_141, STROMAL_141
    import gseapy as gp

    avail_genes = set(expr_df.index)
    matched_stromal = [g for g in STROMAL_141 if g in avail_genes]
    matched_immune = [g for g in IMMUNE_141 if g in avail_genes]
    gene_sets_estimate = {"StromalSignature": matched_stromal, "ImmuneSignature": matched_immune}

    ss_res = gp.ssgsea(
        data=expr_df[met_df["sample_id"]],
        gene_sets=gene_sets_estimate,
        outdir=None,
        sample_norm_method="rank",
        min_size=10,
        max_size=500,
        permutation_num=0,
        no_plot=True,
        verbose=False,
    )
    raw_es = pd.DataFrame(ss_res.resultsOnSamples)
    stromal_scores = raw_es.loc["StromalSignature"].astype(float)
    immune_scores = raw_es.loc["ImmuneSignature"].astype(float)

    # Use Affymetrix reference baseline
    stromal_norm = (stromal_scores - stromal_scores.mean()) / max(stromal_scores.std(), 1e-6) * 1000.0
    immune_norm = (immune_scores - immune_scores.mean()) / max(immune_scores.std(), 1e-6) * 1000.0
    est_scores = stromal_norm + immune_norm
    cos_arg = 0.6049872018 + 0.0001467884 * est_scores
    met_df["Tumor_Purity"] = np.cos(cos_arg).clip(0.0, 1.0).values

    # Match patient's primary tumor GISTIC and clinical label
    cn_df["patient_id"] = cn_df["sample_id"].str[:12]
    clin_df["patient_id"] = clin_df["Sample ID"].str[:12]

    cn_pt_dict = cn_df.set_index("patient_id")["erbb2_copy_number"].to_dict()
    clin_pt_dict = clin_df.set_index("patient_id")["HER2_status"].to_dict()
    ihc_pt_dict = clin_df.set_index("patient_id")["HER2 ihc score"].to_dict()
    fish_pt_dict = clin_df.set_index("patient_id")["HER2 fish status"].to_dict()

    met_df["Patient_Primary_Tumor_GISTIC"] = met_df["patient_id"].map(cn_pt_dict)
    met_df["Patient_Primary_Tumor_Clinical_HER2"] = met_df["patient_id"].map(clin_pt_dict)
    met_df["Patient_Primary_Tumor_IHC"] = met_df["patient_id"].map(ihc_pt_dict)
    met_df["Patient_Primary_Tumor_FISH"] = met_df["patient_id"].map(fish_pt_dict)

    # For prediction: evaluate with matched GISTIC, filling the 1 missing with diploid baseline (0.0)
    met_df["ERBB2_GISTIC_For_Model"] = met_df["Patient_Primary_Tumor_GISTIC"].fillna(0.0)
    met_df["GISTIC_Data_Origin"] = np.where(
        met_df["Patient_Primary_Tumor_GISTIC"].notna(),
        "Matched Primary Tumor GISTIC",
        "Missing in GISTIC (Assumed Diploid 0.0)",
    )

    X_met = met_df[["ERBB2_RNA", "ERBB2_GISTIC_For_Model", "Tumor_Purity", "Amplicon_Z"]].values
    X_met_scaled = scaler.transform(X_met)

    met_df["RF_Probability"] = np.round(rf_model.predict_proba(X_met)[:, 1], 4)
    met_df["LR_Probability"] = np.round(lr_model.predict_proba(X_met_scaled)[:, 1], 4)

    met_df["RF_Predicted_Status"] = np.where(met_df["RF_Probability"] >= opt_thresh_rf, "Positive", "Negative")
    met_df["LR_Predicted_Status"] = np.where(met_df["LR_Probability"] >= opt_thresh_lr, "Positive", "Negative")
    met_df["Rescue_Tier_RF"] = met_df["RF_Probability"].apply(assign_rf_tier)

    print("\n  Metastasis Predictions (7 Samples):")
    met_summary_cols = [
        "sample_id", "patient_id", "Patient_Primary_Tumor_Clinical_HER2", "ERBB2_RNA",
        "Patient_Primary_Tumor_GISTIC", "Tumor_Purity", "Amplicon_Z",
        "RF_Probability", "RF_Predicted_Status", "Rescue_Tier_RF", "LR_Probability", "LR_Predicted_Status"
    ]
    print(met_df[met_summary_cols].to_string(index=False))

    # ------------------------------------------------------------------
    # Step 7.5: Update Master Excel Workbook
    # ------------------------------------------------------------------
    print("\n[Step 7.5] Updating Master Excel Workbook results/her2_classification_results.xlsx...")
    excel_path = os.path.join(output_dir, "her2_classification_results.xlsx")

    # Format Step7_Rescue_Tiers export
    rescue_export = rescue_df[[
        "sample_id", "patient_id", "HER2_status", "HER2_expression_subtype", "discordant_flag",
        "ERBB2_RNA", "ERBB2_GISTIC", "Tumor_Purity", "Amplicon_Z",
        "RF_Probability", "Rescue_Tier_RF", "LR_Probability", "Rescue_Tier_LR"
    ]].sort_values("RF_Probability", ascending=False).reset_index(drop=True)

    # Format Step7_Metastasis export
    met_export = met_df[[
        "sample_id", "patient_id", "sample_type", "Patient_Primary_Tumor_Clinical_HER2",
        "Patient_Primary_Tumor_IHC", "Patient_Primary_Tumor_FISH",
        "ERBB2_RNA", "Patient_Primary_Tumor_GISTIC", "GISTIC_Data_Origin",
        "Tumor_Purity", "Amplicon_Z",
        "RF_Probability", "RF_Predicted_Status", "Rescue_Tier_RF",
        "LR_Probability", "LR_Predicted_Status"
    ]].reset_index(drop=True)

    with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        rescue_export.to_excel(writer, sheet_name="Step7_Rescue_Tiers", index=False)
        met_export.to_excel(writer, sheet_name="Step7_Metastasis", index=False)

        # Format column widths
        for sheet_name in ["Step7_Rescue_Tiers", "Step7_Metastasis"]:
            ws = writer.sheets[sheet_name]
            for col in ws.columns:
                max_len = max(len(str(cell.value or "")) for cell in col)
                ws.column_dimensions[col[0].column_letter].width = min(max(max_len + 3, 14), 45)

    print(f"  Successfully wrote sheets 'Step7_Rescue_Tiers' and 'Step7_Metastasis' to: {excel_path}")

    # Serialize final Step 7 results bundle
    step7_bundle = {
        "rescue_export": rescue_export,
        "met_export": met_export,
        "rf_counts": rf_counts,
        "lr_counts": lr_counts,
        "concordance_df": concordance_df,
    }
    with open(os.path.join(output_dir, "step7_bundle.pkl"), "wb") as f:
        pickle.dump(step7_bundle, f)

    print("\n================================================================================")
    print("STEP 7 EXECUTION COMPLETE. FINAL RESCUE DELIVERABLES GENERATED.")
    print("================================================================================")


if __name__ == "__main__":
    main()
