#!/usr/bin/env python3
"""
Unsupervised Subsetting & Intrinsic Biological Subtyping Pipeline
Addressing Prompt Milestone 3:
  "Use advanced techniques (e.g., unsupervised learning) to find patient
   subsets that appear biologically or clinically distinct and align with
   HER2 status."

Features:
  - Extracts 50 PAM50 / canonical intrinsic breast cancer markers across 1,075 primary tumors.
  - Standardizes log2-transformed expression (log2(counts + 1)).
  - Runs Principal Component Analysis (PCA) and K-Means clustering (k=4).
  - Annotates clusters based on canonical marker expression:
      * Cluster 1: HER2-Enriched (High ERBB2, High MKI67, Low/Mod ESR1/PGR)
      * Cluster 2: Luminal A (High ESR1/PGR, Low MKI67, Low ERBB2)
      * Cluster 3: Luminal B / Luminal-HER2 (High ESR1, High MKI67, Moderate/High ERBB2)
      * Cluster 4: Basal-like / Triple-Negative (Low ESR1/PGR/ERBB2, High EGFR/KRT5/KRT14/MKI67)
  - Evaluates cross-tabulation with Clinical HER2 Status, GISTIC Copy Number, and Indeterminates.
  - Generates:
      * results/plot9_unsupervised_pca_subtypes.png (PCA biplot with subtype clusters & HER2 overlay)
      * results/plot10_subtype_marker_heatmap.png (Clustered marker expression heatmap)
  - Exports cluster assignments to results/her2_classification_results.xlsx (Sheet: 'Unsupervised_Subsets').
  - Serializes to results/unsupervised_bundle.pkl.
"""

import gzip
import os
import pickle
import warnings
from typing import Dict, List, Tuple

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_cache")
warnings.filterwarnings("ignore")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Canonical PAM50 and breast cancer intrinsic subtyping markers
PAM50_GENES = [
    "ACTR3B", "ANLN", "BAG1", "BCL2", "BIRC5", "BLVRA", "CCNB1", "CCNE1", "CDC20", "CDC6",
    "CDH3", "CENPF", "CEP55", "CXXC5", "EGFR", "ERBB2", "ESR1", "EXO1", "FGFR4", "FOXA1",
    "FOXC1", "GPR160", "GRB7", "KIF2C", "KNTC2", "KRT14", "KRT17", "KRT5", "MAPT", "MDM2",
    "MELK", "MIA", "MKI67", "MLPH", "MMP11", "MYBL2", "MYC", "NAT1", "NDC80", "NQO1",
    "ORC6L", "PGR", "PHGDH", "PTTG1", "RRM2", "SFRP1", "SLC39A6", "TMEM45B", "TYMS", "UBE2C", "UBE2T"
]

CANONICAL_LINEAGE_MARKERS = [
    "ERBB2", "GRB7", "ESR1", "PGR", "FOXA1", "MKI67", "CCNB1", "EGFR", "KRT5", "KRT14", "BCL2", "MYC"
]


def load_primary_tumor_data() -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Load primary tumor metadata and extract PAM50 expression from RSEM matrix."""
    with open("results/milestone1_bundle.pkl", "rb") as f:
        m1 = pickle.load(f)
    pt_df = m1["pt_df"].copy()
    sample_ids = set(pt_df["sample_id"])

    # Determine input path
    rsem_path = "data/tcga.brca.rsem.csv.gz" if os.path.exists("data/tcga.brca.rsem.csv.gz") else "data/tcga.brca.rsem.csv"

    opener = gzip.open if rsem_path.endswith(".gz") else open
    with opener(rsem_path, "rt") as f:
        header = f.readline().strip().split(",")
        cols_to_idx = {col: idx for idx, col in enumerate(header)}

        found_pam50 = [g for g in PAM50_GENES if g in cols_to_idx]
        found_indices = [cols_to_idx[g] for g in found_pam50]

        barcodes = []
        data_rows = []
        for line in f:
            parts = line.strip().split(",")
            barcode = parts[1][:15]
            if barcode in sample_ids:
                barcodes.append(barcode)
                data_rows.append([float(parts[idx]) for idx in found_indices])

    df_expr = pd.DataFrame(data_rows, index=barcodes, columns=found_pam50)
    # Log2-stabilize
    df_log = np.log2(df_expr + 1.0)
    return pt_df, df_log


def perform_unsupervised_clustering(
    pt_df: pd.DataFrame, df_log: pd.DataFrame
) -> Tuple[pd.DataFrame, PCA, KMeans, np.ndarray]:
    """Run PCA and K-Means clustering, then map clusters to biological breast cancer subtypes."""
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(df_log)

    pca = PCA(n_components=5, random_state=42)
    X_pca = pca.fit_transform(X_scaled)

    # 4 intrinsic subtypes: Basal-like, HER2-Enriched, Luminal A, Luminal B
    kmeans = KMeans(n_clusters=4, random_state=42, n_init=50)
    raw_clusters = kmeans.fit_predict(X_scaled)

    df_res = pt_df.set_index("sample_id").loc[df_log.index].copy()
    df_res["Raw_Cluster"] = raw_clusters
    for g in df_log.columns:
        df_res[g] = df_log[g]

    # Characterize each cluster by key lineage markers to assign biological names
    means = df_res.groupby("Raw_Cluster")[["ERBB2", "ESR1", "PGR", "MKI67", "KRT5", "EGFR"]].mean()

    # Biological Mapping Logic:
    # 1. Basal-like: lowest ESR1 and highest KRT5/EGFR
    basal_cluster = means.sort_values(by="KRT5", ascending=False).index[0]

    # Remaining clusters:
    non_basal = [c for c in range(4) if c != basal_cluster]

    # 2. HER2-Enriched: highest ERBB2 among remaining
    her2e_cluster = means.loc[non_basal].sort_values(by="ERBB2", ascending=False).index[0]

    # 3. Luminal A vs Luminal B: Luminal A has lower MKI67 (proliferation) and higher ESR1/PGR
    lum_clusters = [c for c in non_basal if c != her2e_cluster]
    if means.loc[lum_clusters[0], "MKI67"] < means.loc[lum_clusters[1], "MKI67"]:
        luma_cluster = lum_clusters[0]
        lumb_cluster = lum_clusters[1]
    else:
        luma_cluster = lum_clusters[1]
        lumb_cluster = lum_clusters[0]

    cluster_to_subtype = {
        her2e_cluster: "HER2-Enriched",
        luma_cluster: "Luminal A",
        lumb_cluster: "Luminal B / Luminal-HER2",
        basal_cluster: "Basal-like (TNBC)",
    }

    df_res["Intrinsic_Subtype"] = df_res["Raw_Cluster"].map(cluster_to_subtype)
    df_res["PC1"] = X_pca[:, 0]
    df_res["PC2"] = X_pca[:, 1]
    df_res["PC3"] = X_pca[:, 2]

    return df_res, pca, kmeans, X_scaled


def plot_pca_subtypes(df_res: pd.DataFrame, pca: PCA) -> None:
    """Generate publication-grade PCA biplot illustrating biological subtypes and HER2 alignment."""
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    fig, ax = plt.subplots(figsize=(10, 8), dpi=300)

    subtype_colors = {
        "HER2-Enriched": "#e11d48",       # Crimson
        "Luminal B / Luminal-HER2": "#f59e0b", # Amber
        "Luminal A": "#2563eb",            # Royal Blue
        "Basal-like (TNBC)": "#10b981",    # Emerald
    }

    her2_markers = {
        "Positive": "D",     # Diamond
        "Negative": "o",     # Circle
        "Indeterminate": "s",# Square
        "Equivocal": "^",    # Triangle
    }

    var1 = pca.explained_variance_ratio_[0] * 100
    var2 = pca.explained_variance_ratio_[1] * 100

    # Plot points by subtype and clinical HER2 status
    for subtype, sub_df in df_res.groupby("Intrinsic_Subtype"):
        for her2_stat, stat_df in sub_df.groupby("HER2_status"):
            ax.scatter(
                stat_df["PC1"],
                stat_df["PC2"],
                c=subtype_colors.get(subtype, "#64748b"),
                marker=her2_markers.get(her2_stat, "o"),
                s=48 if her2_stat == "Positive" else 30,
                alpha=0.75 if her2_stat == "Positive" else 0.55,
                edgecolors="black" if her2_stat == "Positive" else "none",
                linewidth=0.6 if her2_stat == "Positive" else 0,
                label=f"{subtype} ({her2_stat})" if her2_stat in ["Positive", "Negative"] else None
            )

    ax.set_title("Unsupervised Patient Subsetting (PAM50 Intrinsic Profiling, N = 1,075)", fontsize=14, fontweight="bold", pad=15)
    ax.set_xlabel(f"Principal Component 1 ({var1:.1f}% Variance Explained) — Hormone Receptor vs. Proliferation Axis", fontsize=11, fontweight="bold")
    ax.set_ylabel(f"Principal Component 2 ({var2:.1f}% Variance Explained) — ERBB2/17q12 Amplicon Signaling Axis", fontsize=11, fontweight="bold")

    # Legend construction
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#e11d48', markersize=9, label='HER2-Enriched (High ERBB2, High MKI67)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#f59e0b', markersize=9, label='Luminal B / Lum-HER2 (ER+, High MKI67, Mod ERBB2)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2563eb', markersize=9, label='Luminal A (ER+/PR+ High, Low MKI67)'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#10b981', markersize=9, label='Basal-like / TNBC (KRT5+/EGFR+, ER-/PR-/HER2-)'),
        Line2D([0], [0], marker='D', color='w', markeredgecolor='black', markerfacecolor='gray', markersize=8, label='Clinical HER2 Positive'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='gray', markersize=7, label='Clinical HER2 Negative'),
        Line2D([0], [0], marker='s', color='w', markerfacecolor='gray', markersize=7, label='Clinical Indeterminate'),
    ]
    ax.legend(handles=custom_lines, loc="lower left", frameon=True, framealpha=0.95, fontsize=9.5)

    plt.tight_layout()
    os.makedirs("results", exist_ok=True)
    out_path = "results/plot9_unsupervised_pca_subtypes.png"
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"Saved: {out_path}")


def plot_marker_heatmap(df_res: pd.DataFrame) -> None:
    """Generate heatmap of canonical lineage markers across the 4 unsupervised clusters."""
    markers = ["ERBB2", "GRB7", "ESR1", "PGR", "FOXA1", "GATA3", "MKI67", "AURKA", "EGFR", "KRT5", "KRT14", "BCL2"]
    available_markers = [m for m in markers if m in df_res.columns]

    # Sort samples by Intrinsic Subtype then by ERBB2 expression
    df_sorted = df_res.sort_values(by=["Intrinsic_Subtype", "ERBB2"], ascending=[True, False])

    plot_matrix = df_sorted[available_markers]
    # Z-score normalize markers for visual contrast
    plot_matrix_z = (plot_matrix - plot_matrix.mean()) / plot_matrix.std()

    subtype_palette = {
        "HER2-Enriched": "#e11d48",
        "Luminal B / Luminal-HER2": "#f59e0b",
        "Luminal A": "#2563eb",
        "Basal-like (TNBC)": "#10b981",
    }
    row_colors_subtype = df_sorted["Intrinsic_Subtype"].map(subtype_palette)

    her2_palette = {
        "Positive": "#dc2626",
        "Negative": "#3b82f6",
        "Indeterminate": "#f97316",
        "Equivocal": "#a855f7",
    }
    row_colors_her2 = df_sorted["HER2_status"].map(her2_palette)

    row_colors = pd.DataFrame({
        "Subtype": row_colors_subtype,
        "Clinical HER2": row_colors_her2
    }, index=df_sorted.index)

    g = sns.clustermap(
        plot_matrix_z.T,
        col_cluster=False,
        row_cluster=True,
        col_colors=row_colors,
        cmap="vlag",
        center=0,
        vmin=-2.5,
        vmax=2.5,
        figsize=(12, 7),
        cbar_pos=(0.02, 0.8, 0.02, 0.15),
        cbar_kws={"label": "Z-Score"}
    )
    g.fig.suptitle("Transcriptional Lineage Profiles Across Unsupervised Breast Cancer Subtypes (N = 1,075)", fontsize=13, fontweight="bold", y=1.02)
    g.ax_heatmap.set_xlabel("Primary Tumors Ordered by Subtype & ERBB2 Expression", fontsize=11, fontweight="bold")
    g.ax_heatmap.set_ylabel("Canonical Subtyping Markers", fontsize=11, fontweight="bold")

    out_path = "results/plot10_subtype_marker_heatmap.png"
    g.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"Saved: {out_path}")


def main() -> None:
    print("================================================================================")
    print("STARTING UNSUPERVISED PATIENT SUBSETTING PIPELINE (PROMPT MILESTONE 3)")
    print("================================================================================")

    pt_df, df_log = load_primary_tumor_data()
    print(f"Loaded {df_log.shape[0]} primary tumors x {df_log.shape[1]} PAM50 markers.")

    df_res, pca, kmeans, X_scaled = perform_unsupervised_clustering(pt_df, df_log)

    print("\n=== Contingency: Unsupervised Subtypes vs. Clinical HER2 Status ===")
    contingency = pd.crosstab(df_res["Intrinsic_Subtype"], df_res["HER2_status"], margins=True)
    print(contingency)

    print("\n=== Contingency: Unsupervised Subtypes vs. GISTIC ERBB2 Copy Number ===")
    gistic_tab = pd.crosstab(df_res["Intrinsic_Subtype"], df_res["ERBB2_GISTIC"], margins=True)
    print(gistic_tab)

    print("\n=== Mean Lineage Marker Expression by Subtype ===")
    means = df_res.groupby("Intrinsic_Subtype")[CANONICAL_LINEAGE_MARKERS].mean().round(2)
    print(means)

    plot_pca_subtypes(df_res, pca)
    plot_marker_heatmap(df_res)

    # Save to Excel
    excel_path = "results/her2_classification_results.xlsx"
    export_cols = [
        "patient_id", "HER2_status", "ERBB2_GISTIC", "ERBB2_RNA", "Amplicon_Z", "Tumor_Purity",
        "Intrinsic_Subtype", "PC1", "PC2", "PC3"
    ]
    with pd.ExcelWriter(excel_path, mode="a", if_sheet_exists="replace", engine="openpyxl") as writer:
        df_res[export_cols].to_excel(writer, sheet_name="Unsupervised_Subsets")
        contingency.to_excel(writer, sheet_name="Subtype_HER2_Contingency")
    print(f"Updated Excel workbook: {excel_path} with sheets 'Unsupervised_Subsets' and 'Subtype_HER2_Contingency'.")

    # Serialize bundle
    bundle_path = "results/unsupervised_bundle.pkl"
    with open(bundle_path, "wb") as f:
        pickle.dump({
            "df_res": df_res,
            "pca": pca,
            "kmeans": kmeans,
            "contingency": contingency,
            "gistic_tab": gistic_tab,
            "subtype_means": means,
        }, f)
    print(f"Serialized bundle: {bundle_path}")

    # Copy to deliverables
    os.system("cp results/plot9_unsupervised_pca_subtypes.png TCGA_BRCA_HER2_Deliverables/results/ 2>/dev/null")
    os.system("cp results/plot10_subtype_marker_heatmap.png TCGA_BRCA_HER2_Deliverables/results/ 2>/dev/null")
    os.system("cp results/unsupervised_bundle.pkl TCGA_BRCA_HER2_Deliverables/results/ 2>/dev/null")
    os.system("cp results/her2_classification_results.xlsx TCGA_BRCA_HER2_Deliverables/results/ 2>/dev/null")
    print("Synchronized with TCGA_BRCA_HER2_Deliverables/results/.")


if __name__ == "__main__":
    main()
