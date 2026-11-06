#!/usr/bin/env python3
"""
HER2 Purity-Aware Classification Pipeline - Milestone 2 (Step 6)

Covers:
  - 5-Fold Stratified Cross-Validation on Labeled Complete Cases (n = 829).
  - Feature Space: [ERBB2_RNA, ERBB2_GISTIC, Tumor_Purity, Amplicon_Z].
  - Cost-Sensitive Learning:
      1. Logistic Regression (class_weight='balanced')
      2. Random Forest Classifier (500 trees, class_weight='balanced')
  - Out-of-fold probability collection & Youden's J probability threshold shifting.
  - Plot 7: Random Forest Feature Importance Bar Chart.
  - Production model training on 100% of complete-case data (n = 829).
  - Export to results/her2_classification_results.xlsx:
      - Sheet 'Step6_CV_Metrics'
      - Sheet 'Step6_Feature_Importance'
      - Sheet 'Step6_OOF_Predictions'
  - Model serialization to results/production_models.pkl.
"""

import os
import pickle
import warnings
from typing import Dict, List, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_cache")
warnings.filterwarnings("ignore")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    accuracy_score,
    balanced_accuracy_score,
    confusion_matrix,
    f1_score,
    precision_score,
    recall_score,
    roc_curve,
    auc,
)
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler


def compute_metrics_dict(
    y_true: np.ndarray, y_prob: np.ndarray, threshold: float, model_name: str, regime_name: str
) -> Dict[str, float]:
    """Calculate standard classification metrics at a given probability threshold."""
    y_pred = (y_prob >= threshold).astype(int)
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    sens = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    spec = tn / (tn + fp) if (tn + fp) > 0 else 0.0
    prec = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    npv = tn / (tn + fn) if (tn + fn) > 0 else 0.0
    acc = accuracy_score(y_true, y_pred)
    bal_acc = balanced_accuracy_score(y_true, y_pred)
    f1_pos = f1_score(y_true, y_pred, pos_label=1, zero_division=0)
    f1_neg = f1_score(y_true, y_pred, pos_label=0, zero_division=0)
    f1_macro = f1_score(y_true, y_pred, average="macro", zero_division=0)

    return {
        "Model": model_name,
        "Threshold_Regime": regime_name,
        "Threshold_Value": float(threshold),
        "TP": int(tp),
        "FP": int(fp),
        "TN": int(tn),
        "FN": int(fn),
        "Sensitivity_Recall_Positive": float(sens),
        "Specificity_TNR_Negative": float(spec),
        "Precision_PPV_Positive": float(prec),
        "NPV_Negative": float(npv),
        "Accuracy": float(acc),
        "Balanced_Accuracy": float(bal_acc),
        "F1_Positive_HER2Pos": float(f1_pos),
        "F1_Negative_HER2Neg": float(f1_neg),
        "F1_Macro": float(f1_macro),
        "Youdens_J": float(sens + spec - 1.0),
    }


def main() -> None:
    print("================================================================================")
    print("STARTING MILESTONE 2 (STEP 6: PURITY-AWARE ML MODELING)")
    print("================================================================================")

    output_dir = "results"
    artifact_dir = "/Users/galebabu-ali/.gemini/antigravity/brain/0ae2581d-851b-44c2-9836-f4a6ab72a547"
    bundle_path = os.path.join(output_dir, "milestone1_bundle.pkl")

    if not os.path.isfile(bundle_path):
        raise FileNotFoundError(f"Missing milestone 1 bundle at: {bundle_path}")

    with open(bundle_path, "rb") as f:
        bundle = pickle.load(f)

    train_cc_df = bundle["train_cc_df"].copy().reset_index(drop=True)
    features = ["ERBB2_RNA", "ERBB2_GISTIC", "Tumor_Purity", "Amplicon_Z"]

    print(f"\n[Step 6.1] Preparing Complete-Case Feature Matrix (n = {len(train_cc_df)})...")
    print(f"  Features ({len(features)}): {features}")
    print(f"  Target Distribution: Positive = {(train_cc_df['HER2_status'] == 'Positive').sum()}, Negative = {(train_cc_df['HER2_status'] == 'Negative').sum()} (Ratio = 1 : {(train_cc_df['HER2_status'] == 'Negative').sum() / (train_cc_df['HER2_status'] == 'Positive').sum():.2f})")

    X = train_cc_df[features].values
    y = (train_cc_df["HER2_status"] == "Positive").astype(int).values

    # ------------------------------------------------------------------
    # 5-Fold Stratified Cross-Validation Loop
    # ------------------------------------------------------------------
    print("\n[Step 6.2] Running 5-Fold Stratified Cross-Validation...")
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)

    oof_lr_probs = np.zeros(len(y), dtype=float)
    oof_rf_probs = np.zeros(len(y), dtype=float)
    fold_indices = np.zeros(len(y), dtype=int)

    rf_fold_importances = []
    fold_perf_rows = []

    for fold, (train_idx, val_idx) in enumerate(skf.split(X, y), 1):
        X_train, y_train = X[train_idx], y[train_idx]
        X_val, y_val = X[val_idx], y[val_idx]
        fold_indices[val_idx] = fold

        # 1. Standardize features inside CV fold for Logistic Regression
        scaler = StandardScaler()
        X_train_scaled = scaler.fit_transform(X_train)
        X_val_scaled = scaler.transform(X_val)

        # 2. Train Logistic Regression with cost-sensitive class weights
        lr_fold = LogisticRegression(
            class_weight="balanced",
            max_iter=1000,
            solver="lbfgs",
            random_state=42 + fold,
        )
        lr_fold.fit(X_train_scaled, y_train)
        val_lr_prob = lr_fold.predict_proba(X_val_scaled)[:, 1]
        oof_lr_probs[val_idx] = val_lr_prob

        # 3. Train Random Forest (500 trees, class_weight='balanced')
        rf_fold = RandomForestClassifier(
            n_estimators=500,
            class_weight="balanced",
            random_state=42 + fold,
            n_jobs=-1,
        )
        rf_fold.fit(X_train, y_train)
        val_rf_prob = rf_fold.predict_proba(X_val)[:, 1]
        oof_rf_probs[val_idx] = val_rf_prob

        rf_fold_importances.append(rf_fold.feature_importances_)

        # Per-fold ROC AUC
        fpr_lr_f, tpr_lr_f, _ = roc_curve(y_val, val_lr_prob)
        auc_lr_f = auc(fpr_lr_f, tpr_lr_f)

        fpr_rf_f, tpr_rf_f, _ = roc_curve(y_val, val_rf_prob)
        auc_rf_f = auc(fpr_rf_f, tpr_rf_f)

        fold_perf_rows.append({
            "Fold": fold,
            "N_Train": len(train_idx),
            "N_Val": len(val_idx),
            "Val_Positives": int(y_val.sum()),
            "Val_Negatives": int(len(val_idx) - y_val.sum()),
            "LR_AUC": float(auc_lr_f),
            "RF_AUC": float(auc_rf_f),
        })

    fold_perf_df = pd.DataFrame(fold_perf_rows)
    print("  Cross-Validation Fold Summary:")
    print(fold_perf_df.to_string(index=False))
    print(f"  Mean LR AUC: {fold_perf_df['LR_AUC'].mean():.4f} +/- {fold_perf_df['LR_AUC'].std():.4f}")
    print(f"  Mean RF AUC: {fold_perf_df['RF_AUC'].mean():.4f} +/- {fold_perf_df['RF_AUC'].std():.4f}")

    # ------------------------------------------------------------------
    # Step 6.3: Probability Threshold Shifting via Youden's J
    # ------------------------------------------------------------------
    print("\n[Step 6.3] Calculating Out-of-Fold ROC Curves & Youden's J Probability Thresholds...")
    # Logistic Regression OOF ROC & Youden
    fpr_lr, tpr_lr, thresholds_lr = roc_curve(y, oof_lr_probs)
    auc_lr_oof = auc(fpr_lr, tpr_lr)
    j_lr = tpr_lr - fpr_lr
    opt_idx_lr = int(np.argmax(j_lr))
    opt_thresh_lr = float(thresholds_lr[opt_idx_lr])
    opt_j_lr = float(j_lr[opt_idx_lr])

    # Random Forest OOF ROC & Youden
    fpr_rf, tpr_rf, thresholds_rf = roc_curve(y, oof_rf_probs)
    auc_rf_oof = auc(fpr_rf, tpr_rf)
    j_rf = tpr_rf - fpr_rf
    opt_idx_rf = int(np.argmax(j_rf))
    opt_thresh_rf = float(thresholds_rf[opt_idx_rf])
    opt_j_rf = float(j_rf[opt_idx_rf])

    print(f"  LR Out-of-Fold AUC: {auc_lr_oof:.4f} | Optimal Shifted Threshold: {opt_thresh_lr:.4f} (Max Youden J = {opt_j_lr:.4f})")
    print(f"  RF Out-of-Fold AUC: {auc_rf_oof:.4f} | Optimal Shifted Threshold: {opt_thresh_rf:.4f} (Max Youden J = {opt_j_rf:.4f})")

    # Evaluate classification metrics across default (0.50) and shifted thresholds
    m_lr_default = compute_metrics_dict(y, oof_lr_probs, 0.50, "Logistic Regression", "Default Threshold (0.50)")
    m_lr_shifted = compute_metrics_dict(y, oof_lr_probs, opt_thresh_lr, "Logistic Regression", f"Shifted Youden Threshold ({opt_thresh_lr:.4f})")
    m_rf_default = compute_metrics_dict(y, oof_rf_probs, 0.50, "Random Forest (500 Trees)", "Default Threshold (0.50)")
    m_rf_shifted = compute_metrics_dict(y, oof_rf_probs, opt_thresh_rf, "Random Forest (500 Trees)", f"Shifted Youden Threshold ({opt_thresh_rf:.4f})")

    metrics_comp_df = pd.DataFrame([m_lr_default, m_lr_shifted, m_rf_default, m_rf_shifted])
    metrics_comp_df["ROC_AUC_OOF"] = [auc_lr_oof, auc_lr_oof, auc_rf_oof, auc_rf_oof]

    print("\n--- OUT-OF-FOLD PERFORMANCE METRICS COMPARISON ---")
    display_cols = ["Model", "Threshold_Regime", "Threshold_Value", "TP", "FP", "TN", "FN", "Sensitivity_Recall_Positive", "Specificity_TNR_Negative", "Precision_PPV_Positive", "F1_Positive_HER2Pos", "Balanced_Accuracy", "ROC_AUC_OOF"]
    print(metrics_comp_df[display_cols].to_string(index=False))

    # ------------------------------------------------------------------
    # Step 6.4: PLOT 7 - Random Forest Feature Importance Bar Chart
    # ------------------------------------------------------------------
    print("\n[Step 6.4] Rendering Plot 7 (Random Forest Feature Importance)...")
    rf_importances_arr = np.array(rf_fold_importances)  # shape: (5, 4)
    imp_mean = rf_importances_arr.mean(axis=0)
    imp_std = rf_importances_arr.std(axis=0)

    feat_imp_df = pd.DataFrame({
        "Feature": features,
        "Mean_Importance": imp_mean,
        "Std_Importance": imp_std,
        "Mean_Percent": imp_mean * 100.0,
        "Std_Percent": imp_std * 100.0,
    }).sort_values("Mean_Importance", ascending=True).reset_index(drop=True)

    # Feature descriptions
    feat_desc = {
        "ERBB2_RNA": "Direct Target Expression (TMM-voom Log2-CPM)",
        "Amplicon_Z": "17q12 Flanking Co-Expression (GRB7, PGAP3, STARD3, MIEN1)",
        "Tumor_Purity": "ESTIMATE Purity Score (Stromal / Immune Infiltration)",
        "ERBB2_GISTIC": "Discrete GISTIC Copy Number (-2 to +2)",
    }
    feat_imp_df["Feature_Description"] = feat_imp_df["Feature"].map(feat_desc)

    print("  Feature Importances across 5 Folds:")
    for _, r in feat_imp_df.sort_values("Mean_Importance", ascending=False).iterrows():
        print(f"    * {r['Feature']:<14}: {r['Mean_Percent']:.2f}% +/- {r['Std_Percent']:.2f}% ({r['Feature_Description']})")

    # Build Plot 7 Figure
    fig, ax = plt.subplots(figsize=(10.5, 5.6), dpi=300)
    sns.set_style("whitegrid")

    palette_bars = ["#6BAED6", "#74C476", "#FD8D3C", "#E6550D"]
    bars = ax.barh(
        y=np.arange(len(feat_imp_df)),
        width=feat_imp_df["Mean_Percent"],
        xerr=feat_imp_df["Std_Percent"],
        color=palette_bars,
        edgecolor="#333333",
        linewidth=1.2,
        capsize=5,
        height=0.55,
        alpha=0.90,
    )

    ax.set_yticks(np.arange(len(feat_imp_df)))
    ax.set_yticklabels(
        [feat for feat in feat_imp_df["Feature"]],
        fontsize=10.5,
        fontweight="bold",
    )
    ax.set_xlabel("Mean Feature Importance / MDI (%) across 5 Stratified Folds", fontsize=10.5, fontweight="bold")
    ax.set_title(
        "Plot 7: Random Forest Feature Importance (500 Trees, class_weight='balanced')\n"
        f"Relative Contribution across 4-Feature Space (n = {len(train_cc_df)} Complete Training Cases)",
        fontsize=11.5,
        fontweight="bold",
        pad=12,
    )
    ax.set_xlim(0, 50)

    # Annotate bar values
    for i, bar in enumerate(bars):
        w = bar.get_width()
        err = feat_imp_df.loc[i, "Std_Percent"]
        ax.text(
            w + err + 0.8,
            bar.get_y() + bar.get_height() / 2.0,
            f"{w:.2f}% ± {err:.2f}%",
            va="center",
            ha="left",
            fontsize=9.5,
            fontweight="bold",
            color="#222222",
        )

    ax.grid(True, linestyle=":", alpha=0.5, axis="x")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    plot7_path = os.path.join(output_dir, "plot7_rf_feature_importance.png")
    fig.savefig(plot7_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved Plot 7 to {plot7_path}")

    # Copy Plot 7 to artifact directory
    dst_plot7 = os.path.join(artifact_dir, "plot7_rf_feature_importance.png")
    with open(plot7_path, "rb") as rf_in, open(dst_plot7, "wb") as wf_out:
        wf_out.write(rf_in.read())
    print(f"  Synced Plot 7 to artifact directory: {dst_plot7}")

    # ------------------------------------------------------------------
    # Step 6.5: Train Final Production Models on 100% of Complete Dataset (n = 829)
    # ------------------------------------------------------------------
    print("\n[Step 6.5] Training Final Production Models on 100% of Complete Cases (n = 829)...")
    prod_scaler = StandardScaler()
    X_scaled_full = prod_scaler.fit_transform(X)

    prod_lr = LogisticRegression(
        class_weight="balanced",
        max_iter=1000,
        solver="lbfgs",
        random_state=42,
    )
    prod_lr.fit(X_scaled_full, y)

    prod_rf = RandomForestClassifier(
        n_estimators=500,
        class_weight="balanced",
        random_state=42,
        n_jobs=-1,
    )
    prod_rf.fit(X, y)

    print("  Production models successfully trained.")
    print(f"  Production RF Feature Importances: {dict(zip(features, np.round(prod_rf.feature_importances_*100, 2)))}")

    # Serialize production models for Step 7
    prod_bundle = {
        "features": features,
        "scaler": prod_scaler,
        "lr_model": prod_lr,
        "rf_model": prod_rf,
        "opt_thresh_lr": opt_thresh_lr,
        "opt_thresh_rf": opt_thresh_rf,
        "auc_lr_oof": auc_lr_oof,
        "auc_rf_oof": auc_rf_oof,
        "metrics_comp_df": metrics_comp_df,
        "feat_imp_df": feat_imp_df,
    }
    prod_bundle_path = os.path.join(output_dir, "production_models.pkl")
    with open(prod_bundle_path, "wb") as f:
        pickle.dump(prod_bundle, f)
    print(f"  Saved serialized production models to: {prod_bundle_path}")

    # ------------------------------------------------------------------
    # Step 6.6: Export Out-of-Fold Predictions & Update Master Excel
    # ------------------------------------------------------------------
    print("\n[Step 6.6] Exporting OOF Predictions & Updating Master Excel Workbook...")
    oof_df = train_cc_df[[
        "sample_id", "patient_id", "HER2_status", "HER2_expression_subtype", "discordant_flag",
        "ERBB2_RNA", "ERBB2_GISTIC", "Tumor_Purity", "Amplicon_Z"
    ]].copy()

    oof_df["CV_Fold"] = fold_indices
    oof_df["LR_OOF_Prob"] = np.round(oof_lr_probs, 4)
    oof_df["RF_OOF_Prob"] = np.round(oof_rf_probs, 4)
    oof_df["LR_Pred_Default_0.50"] = np.where(oof_lr_probs >= 0.50, "Positive", "Negative")
    oof_df["LR_Pred_Shifted_Youden"] = np.where(oof_lr_probs >= opt_thresh_lr, "Positive", "Negative")
    oof_df["RF_Pred_Default_0.50"] = np.where(oof_rf_probs >= 0.50, "Positive", "Negative")
    oof_df["RF_Pred_Shifted_Youden"] = np.where(oof_rf_probs >= opt_thresh_rf, "Positive", "Negative")

    excel_path = os.path.join(output_dir, "her2_classification_results.xlsx")

    # Load existing workbook sheets and write new ones
    with pd.ExcelWriter(excel_path, engine="openpyxl", mode="a", if_sheet_exists="replace") as writer:
        metrics_comp_df.to_excel(writer, sheet_name="Step6_CV_Metrics", index=False)
        feat_imp_df.sort_values("Mean_Importance", ascending=False).to_excel(writer, sheet_name="Step6_Feature_Importance", index=False)
        oof_df.to_excel(writer, sheet_name="Step6_OOF_Predictions", index=False)

        # Re-format column widths
        for sheet_name in ["Step6_CV_Metrics", "Step6_Feature_Importance", "Step6_OOF_Predictions"]:
            ws = writer.sheets[sheet_name]
            for col in ws.columns:
                max_len = max(len(str(cell.value or "")) for cell in col)
                ws.column_dimensions[col[0].column_letter].width = min(max(max_len + 3, 14), 50)

    print(f"  Successfully updated master Excel workbook at: {excel_path}")
    print("\n================================================================================")
    print("MILESTONE 2 EXECUTION COMPLETE. REACHED STOPPING POINT 3.")
    print("================================================================================")


if __name__ == "__main__":
    main()
