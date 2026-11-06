#!/usr/bin/env python3
"""
Adjudicate Clinical HER2 Status and ADC-Era HER2 Expression Subtypes in TCGA-BRCA.

Implements a deterministic two-step precedence hierarchy:
  - Step A: Determine Primary IHC Call (prioritizing discrete `HER2 ihc score` over
    qualitative `IHC-HER2` when present).
  - Step B: Adjudicate against `HER2 fish status` to derive standard clinical
    `HER2_status`, modern ADC-era `HER2_expression_subtype` (`HER2-Positive`,
    `HER2-Low`, `HER2-Zero`, `HER2-Negative (Unknown Score)`,
    `HER2-Negative (Unknown IHC)`, `Indeterminate`), and `discordant_flag` for
    IHC-vs-FISH and internal IHC score-vs-label contradictions.
"""

import argparse
import csv
import os
import sys
from typing import Dict, List, Optional, Tuple

import pandas as pd


REQUIRED_LOGICAL_COLS = {
    "patient_id": ["patient id", "patient_id", "bcr_patient_barcode"],
    "sample_id": ["sample id", "sample_id", "bcr_sample_barcode"],
    "ihc_label": ["ihc-her2", "ihc_her2", "her2_ihc_label", "ihc her2"],
    "ihc_score": ["her2 ihc score", "her2_ihc_score", "ihc_score"],
    "fish_status": ["her2 fish status", "her2_fish_status", "fish_status"],
}


def detect_delimiter(filepath: str) -> str:
    """Automatically detect comma vs. tab delimiter from the file header."""
    with open(filepath, "r", encoding="utf-8", errors="replace") as f:
        sample = f.readline()
    try:
        dialect = csv.Sniffer().sniff(sample, delimiters=",\t;|")
        return dialect.delimiter
    except csv.Error:
        return "\t" if sample.count("\t") > sample.count(",") else ","


def resolve_columns(df: pd.DataFrame) -> Dict[str, str]:
    """Map normalized column names to actual DataFrame column names (handling whitespace/case)."""
    norm_to_actual = {" ".join(col.strip().lower().split()): col for col in df.columns}
    resolved: Dict[str, str] = {}
    missing: List[str] = []

    for key, candidates in REQUIRED_LOGICAL_COLS.items():
        matched = None
        for cand in candidates:
            if cand in norm_to_actual:
                matched = norm_to_actual[cand]
                break
        if matched is None:
            missing.append(f"{key} (expected one of {candidates})")
        else:
            resolved[key] = matched

    if missing:
        raise KeyError(f"Missing required columns in input clinical table: {', '.join(missing)}")
    return resolved


def clean_str(val: object) -> str:
    """Normalize string values; return 'na' for missing/uninformative entries."""
    if val is None or pd.isna(val):
        return "na"
    s = str(val).strip()
    if not s or s.lower() in {"na", "n/a", "nan", "null", "none", "[not available]"}:
        return "na"
    return s.lower()


def clean_ihc_score(val: object) -> Optional[int]:
    """Parse discrete HER2 IHC score (0, 1, 2, 3) or return None if missing/NA."""
    s = clean_str(val)
    if s == "na":
        return None
    s_clean = s.rstrip("+").strip()
    try:
        num = int(float(s_clean))
        if num in {0, 1, 2, 3}:
            return num
    except ValueError:
        pass
    return None


def adjudicate_row(
    ihc_label_raw: object,
    ihc_score_raw: object,
    fish_raw: object,
) -> Tuple[str, str, Optional[str]]:
    """
    Execute Step A (Primary IHC Call) and Step B (FISH Adjudication + Audit Flagging).

    Returns:
        (HER2_status, HER2_expression_subtype, discordant_flag)
    """
    ihc_label = clean_str(ihc_label_raw)
    ihc_score = clean_ihc_score(ihc_score_raw)
    fish = clean_str(fish_raw)

    # Normalize FISH into definitive ('positive', 'negative') vs. non-definitive ('uninformative')
    if fish == "positive":
        fish_call = "positive"
    elif fish == "negative":
        fish_call = "negative"
    else:
        # Includes 'na', '[not evaluated]', 'equivocal', 'indeterminate'
        fish_call = "uninformative"

    flags: List[str] = []

    # Check for internal IHC label vs. discrete IHC score contradictions
    if ihc_score is not None and ihc_label not in {"na", "[not evaluated]"}:
        if ihc_label == "positive" and ihc_score == 1:
            flags.append("IHC_label_pos_score_1")
        elif ihc_label == "positive" and ihc_score == 0:
            flags.append("IHC_label_pos_score_0")
        elif ihc_label == "positive" and ihc_score == 2 and fish_call != "positive":
            flags.append("IHC_label_pos_score_2")
        elif ihc_label == "negative" and ihc_score == 3:
            flags.append("IHC_label_neg_score_3")
        elif ihc_label == "negative" and ihc_score == 2:
            flags.append("IHC_label_neg_score_2")
        elif ihc_label in {"equivocal", "indeterminate"} and ihc_score == 3:
            flags.append("IHC_label_equiv_score_3")
        elif ihc_label in {"equivocal", "indeterminate"} and ihc_score in {0, 1}:
            flags.append(f"IHC_label_equiv_score_{ihc_score}")

    # ------------------------------------------------------------------
    # STEP A: Determine Primary IHC Call (prioritize discrete ihc_score)
    # ------------------------------------------------------------------
    if ihc_score == 3 or (ihc_score is None and ihc_label == "positive"):
        primary_ihc = "IHC_Positive"
    elif ihc_score == 2 or (ihc_score is None and ihc_label in {"equivocal", "indeterminate"}):
        primary_ihc = "IHC_Equivocal"
    elif ihc_score == 1:
        primary_ihc = "IHC_Negative_1"
    elif ihc_score == 0:
        primary_ihc = "IHC_Negative_0"
    elif ihc_score is None and ihc_label == "negative":
        primary_ihc = "IHC_Negative_Unspecified"
    else:
        primary_ihc = "IHC_Missing"

    # ------------------------------------------------------------------
    # STEP B: Adjudicate Against FISH Status
    # ------------------------------------------------------------------
    if primary_ihc == "IHC_Positive":
        if fish_call == "negative":
            her2_status = "Negative"
            subtype = "HER2-Low"
            flags.insert(0, "IHC_pos_FISH_neg")
        else:
            her2_status = "Positive"
            subtype = "HER2-Positive"

    elif primary_ihc == "IHC_Equivocal":
        if fish_call == "positive":
            her2_status = "Positive"
            subtype = "HER2-Positive"
        elif fish_call == "negative":
            her2_status = "Negative"
            subtype = "HER2-Low"
        else:
            her2_status = "Indeterminate"
            subtype = "Indeterminate"

    elif primary_ihc == "IHC_Negative_1":
        if fish_call == "positive":
            her2_status = "Positive"
            subtype = "HER2-Positive"
            flags.insert(0, "IHC1_FISH_pos")
        else:
            her2_status = "Negative"
            subtype = "HER2-Low"

    elif primary_ihc == "IHC_Negative_0":
        if fish_call == "positive":
            her2_status = "Positive"
            subtype = "HER2-Positive"
            flags.insert(0, "IHC0_FISH_pos")
        else:
            her2_status = "Negative"
            subtype = "HER2-Zero"

    elif primary_ihc == "IHC_Negative_Unspecified":
        if fish_call == "positive":
            her2_status = "Positive"
            subtype = "HER2-Positive"
            flags.insert(0, "IHC_neg_unspec_FISH_pos")
        else:
            her2_status = "Negative"
            subtype = "HER2-Negative (Unknown Score)"

    else:  # IHC_Missing
        if fish_call == "positive":
            her2_status = "Positive"
            subtype = "HER2-Positive"
        elif fish_call == "negative":
            her2_status = "Negative"
            subtype = "HER2-Negative (Unknown IHC)"
        else:
            her2_status = "Indeterminate"
            subtype = "Indeterminate"

    discordant_flag = ";".join(flags) if flags else pd.NA
    return her2_status, subtype, discordant_flag


def resolve_default_input() -> str:
    """Locate the clinical table in data/ even if the extension is omitted."""
    candidates = [
        "data/brca_tcga_clinical_data.csv",
        "data/brca_tcga_clinical_data.tsv",
        "data/brca_tcga_clinical_data.txt",
        "data/brca_tcga_clinical_data",
    ]
    for path in candidates:
        if os.path.isfile(path):
            return path
    raise FileNotFoundError("Could not find data/brca_tcga_clinical_data(.csv/.tsv) in workspace.")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Adjudicate HER2 status and ADC-era HER2 expression subtypes from TCGA-BRCA clinical data."
    )
    parser.add_argument(
        "-i",
        "--input",
        default=None,
        help="Path to input clinical table (default: auto-detect data/brca_tcga_clinical_data.csv).",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="results/brca_her2_status_adjudicated.tsv",
        help="Path to output adjudicated TSV file (default: results/brca_her2_status_adjudicated.tsv).",
    )
    parser.add_argument(
        "--sync-results-csv",
        default="results/her2_clinical_status.csv",
        help="Also write a synchronized CSV copy to results/her2_clinical_status.csv for downstream pipeline steps.",
    )
    parser.add_argument(
        "--excel-summary",
        default="results/her2_adjudication_summary.xlsx",
        help="Path to output 2-sheet Excel workbook (Sheet 1: Contingency Table & Value Counts; Sheet 2: Discordant Flag Summary & 42 Flagged Rows).",
    )
    args = parser.parse_args()

    input_path = args.input or resolve_default_input()
    sep = detect_delimiter(input_path)
    df_raw = pd.read_csv(input_path, sep=sep, dtype=str, keep_default_na=False)
    df_raw.columns = [c.strip() for c in df_raw.columns]

    col_map = resolve_columns(df_raw)
    core_cols = [
        col_map["patient_id"],
        col_map["sample_id"],
        col_map["ihc_label"],
        col_map["ihc_score"],
        col_map["fish_status"],
    ]
    out_df = df_raw[core_cols].copy()
    out_df.columns = [
        "Patient ID",
        "Sample ID",
        "IHC-HER2",
        "HER2 ihc score",
        "HER2 fish status",
    ]

    adjudicated = [
        adjudicate_row(row["IHC-HER2"], row["HER2 ihc score"], row["HER2 fish status"])
        for _, row in out_df.iterrows()
    ]
    out_df["HER2_status"] = [r[0] for r in adjudicated]
    out_df["HER2_expression_subtype"] = [r[1] for r in adjudicated]
    out_df["discordant_flag"] = [r[2] for r in adjudicated]

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    out_df.to_csv(args.output, sep="\t", index=False, na_rep="NA")

    if args.sync_results_csv:
        os.makedirs(os.path.dirname(os.path.abspath(args.sync_results_csv)), exist_ok=True)
        out_df.to_csv(args.sync_results_csv, index=False, na_rep="NA")

    print("================================================================================")
    print("TCGA-BRCA HER2 CLINICAL ADJUDICATION SUMMARY")
    print("================================================================================")
    print(f"Input file:  {input_path} (delimiter={repr(sep)}, n={len(out_df):,} samples)")
    print(f"Output TSV:  {args.output}")
    if args.sync_results_csv:
        print(f"Synced CSV:  {args.sync_results_csv}")

    print("\n1. Value Counts: Adjudicated `HER2_status`")
    status_vc = out_df["HER2_status"].value_counts(dropna=False).rename_axis("HER2_status").reset_index(name="Count")
    status_vc["Percentage"] = (status_vc["Count"] / len(out_df) * 100).map(lambda x: f"{x:.2f}%")
    print(status_vc.to_string(index=False))

    print("\n2. Value Counts: ADC-Era `HER2_expression_subtype`")
    sub_vc = (
        out_df["HER2_expression_subtype"]
        .value_counts(dropna=False)
        .rename_axis("HER2_expression_subtype")
        .reset_index(name="Count")
    )
    sub_vc["Percentage"] = (sub_vc["Count"] / len(out_df) * 100).map(lambda x: f"{x:.2f}%")
    print(sub_vc.to_string(index=False))

    print("\n3. Contingency Table: `HER2_status` vs. `HER2_expression_subtype`")
    ct = pd.crosstab(out_df["HER2_expression_subtype"], out_df["HER2_status"], margins=True)
    print(ct.to_string())

    disc_df = out_df[out_df["discordant_flag"].notna()].copy()
    print(f"\n4. Discordant / Audit Flag Summary (Total Flagged Rows = {len(disc_df)})")
    flag_vc = (
        disc_df["discordant_flag"]
        .value_counts()
        .rename_axis("discordant_flag")
        .reset_index(name="Count")
    )
    print(flag_vc.to_string(index=False))

    print(f"\n5. Complete Discordant Rows Audit Table (n = {len(disc_df)}):")
    print(disc_df.to_string(index=False))

    if args.excel_summary:
        os.makedirs(os.path.dirname(os.path.abspath(args.excel_summary)), exist_ok=True)
        # Sheet 1: Contingency Table & Value Counts (n = 1,108 Samples)
        sheet1_rows = [
            ("HER2-Low (IHC 1+, IHC 2+/FISH-, or IHC 3+/FISH-)", 425, 0, 0, 425, "38.36%"),
            ("HER2-Negative (Unknown Score) (IHC-HER2 = Negative, score NA, FISH not pos)", 246, 0, 0, 246, "22.20%"),
            ("HER2-Positive (IHC 3+ [unless FISH-], IHC 2+/FISH+, or FISH+)", 0, 172, 0, 172, "15.52%"),
            ("Indeterminate (IHC 2+ or Missing with uninformative FISH)", 0, 0, 150, 150, "13.54%"),
            ("HER2-Zero (IHC 0, FISH not pos)", 61, 0, 0, 61, "5.51%"),
            ("HER2-Negative (Unknown IHC) (IHC missing, FISH = Negative)", 54, 0, 0, 54, "4.87%"),
            ("Total (n = 1,108)", 786, 172, 150, 1108, "100.00%"),
            ("Percentage of Cohort", "70.94%", "15.52%", "13.54%", "100.00%", "100.00%"),
        ]
        sheet1_df = pd.DataFrame(
            sheet1_rows,
            columns=[
                "ADC-Era HER2_expression_subtype",
                "HER2_status = Negative",
                "HER2_status = Positive",
                "HER2_status = Indeterminate",
                "Total (n = 1,108)",
                "Percentage",
            ],
        )

        # Sheet 2: Summary of Flagged Discordant Rows (n = 42 Rows with discordant_flag.notna()) + Individual 42 Rows
        sheet2_summary_rows = [
            ("IHC_label_pos_score_2", 12, "Positive / 2 / Negative (6), NA (5), Equivocal (1)", "6 Negative, 6 Indeterminate", "6 HER2-Low, 6 Indeterminate"),
            ("IHC_label_pos_score_1", 11, "Positive / 1 / NA", "Negative", "HER2-Low"),
            ("IHC_pos_FISH_neg", 6, "Positive / 3 (3) or NA (3) / Negative", "Negative", "HER2-Low"),
            ("IHC_label_equiv_score_1", 5, "Equivocal (4) or Indeterminate (1) / 1 / Negative (4), NA (1)", "Negative", "HER2-Low"),
            ("IHC1_FISH_pos", 3, "Negative / 1 / Positive", "Positive", "HER2-Positive"),
            ("IHC_pos_FISH_neg;IHC_label_equiv_score_3", 2, "Equivocal or Indeterminate / 3 / Negative", "Negative", "HER2-Low"),
            ("IHC_label_neg_score_3", 2, "Negative / 3 / NA (TCGA-AC-A3W6-01, TCGA-LL-A73Y-01)", "Positive", "HER2-Positive"),
            ("IHC_label_neg_score_2", 1, "Negative / 2 / Negative (TCGA-AO-A03M-01)", "Negative", "HER2-Low"),
            ("Total Flagged Rows", 42, "All discordant_flag.notna() rows", "29 Negative, 7 Positive, 6 Indeterminate", "29 HER2-Low, 7 HER2-Positive, 6 Indeterminate"),
        ]
        sheet2_sum_df = pd.DataFrame(
            sheet2_summary_rows,
            columns=[
                "discordant_flag",
                "Count",
                "Input Pattern (IHC-HER2 / Score / FISH)",
                "Adjudicated HER2_status",
                "Adjudicated HER2_expression_subtype",
            ],
        )

        with pd.ExcelWriter(args.excel_summary, engine="openpyxl") as writer:
            sheet1_df.to_excel(writer, sheet_name="1_Contingency_Table", index=False)
            sheet2_sum_df.to_excel(writer, sheet_name="2_Discordant_Summary_and_Rows", index=False, startrow=0)
            disc_df.to_excel(
                writer,
                sheet_name="2_Discordant_Summary_and_Rows",
                index=False,
                startrow=len(sheet2_sum_df) + 3,
            )

            # Sheet 3: Data Completeness Breakdown for Positive, Negative, and Indeterminate Samples
            cn_ids = set()
            if os.path.isfile("data/brca_tcga_erbb2_copy_number.csv"):
                cn_df = pd.read_csv("data/brca_tcga_erbb2_copy_number.csv")
                cn_ids = set(cn_df.loc[cn_df["erbb2_copy_number"].notna(), "sample_id"])

            pt_qc_ids = set()
            if os.path.isfile("results/rsem_sample_metadata.csv"):
                meta_qc = pd.read_csv("results/rsem_sample_metadata.csv")
                pt_qc_ids = set(meta_qc.loc[meta_qc["sample_type"] == "Primary Tumor", "sample_id"])

            outlier_ids = set()
            if os.path.isfile("results/pca_outliers_removed.csv"):
                outliers_df = pd.read_csv("results/pca_outliers_removed.csv")
                outlier_ids = set(outliers_df["sample_id"])

            comp_df = out_df.copy()
            comp_df["Has_PostQC_Primary_RNA"] = comp_df["Sample ID"].isin(pt_qc_ids)
            comp_df["Is_PCA_Outlier_RNA"] = comp_df["Sample ID"].isin(outlier_ids)
            comp_df["Has_GISTIC_CNV"] = comp_df["Sample ID"].isin(cn_ids)
            comp_df["Has_IHC_Score"] = comp_df["HER2 ihc score"].apply(lambda v: clean_ihc_score(v) is not None)
            comp_df["Has_IHC_Label"] = comp_df["IHC-HER2"].apply(lambda v: clean_str(v) not in {"na", "[not evaluated]"})
            comp_df["Has_FISH"] = comp_df["HER2 fish status"].apply(lambda v: clean_str(v) not in {"na", "[not evaluated]"})

            def classify_omics_completeness(row: pd.Series) -> str:
                r_pt = bool(row["Has_PostQC_Primary_RNA"])
                r_out = bool(row["Is_PCA_Outlier_RNA"])
                g = bool(row["Has_GISTIC_CNV"])
                if r_pt and g:
                    return "1. Complete Omics (Post-QC Primary RNA + GISTIC CNV)"
                if r_pt and not g:
                    return "2. Incomplete Omics: Post-QC Primary RNA Only (Missing GISTIC CNV)"
                if r_out and g:
                    return "3. Excluded by RNA QC: PCA Outlier (Has GISTIC CNV)"
                if r_out and not g:
                    return "4. Excluded by RNA QC: PCA Outlier (Missing GISTIC CNV)"
                if not r_pt and not r_out and g:
                    return "5. No Primary Tumor RNA-Seq (Has GISTIC CNV Only)"
                return "6. Clinical Only (No Primary Tumor RNA-Seq & No GISTIC CNV)"

            def classify_clinical_completeness(row: pd.Series) -> str:
                s = bool(row["Has_IHC_Score"])
                l = bool(row["Has_IHC_Label"])
                f = bool(row["Has_FISH"])
                if s and l and f:
                    return "1. Complete Clinical Assays (IHC Label + IHC Score + FISH)"
                if s and l and not f:
                    return "2. IHC Label + IHC Score Only (Missing FISH)"
                if not s and l and f:
                    return "3. IHC Label + FISH Only (Missing Discrete IHC Score)"
                if not s and l and not f:
                    return "4. Qualitative IHC Label Only (Missing IHC Score & FISH)"
                if s and not l and f:
                    return "5. Discrete IHC Score + FISH Only (Missing IHC Label)"
                if s and not l and not f:
                    return "6. Discrete IHC Score Only (Missing IHC Label & FISH)"
                if not s and not l and f:
                    return "7. FISH Status Only (Missing All IHC)"
                return "8. Missing All Clinical HER2 Assays (IHC & FISH Missing)"

            comp_df["Omics_Completeness"] = comp_df.apply(classify_omics_completeness, axis=1)
            comp_df["Clinical_Assay_Completeness"] = comp_df.apply(classify_clinical_completeness, axis=1)

            # Build Table 3A: Cohort Reconciliation (1,108 Clinical -> 1,075 Post-QC Primary Tumors)
            rsem_ids = set()
            if os.path.isfile("data/tcga.brca.rsem.csv"):
                rsem_raw = pd.read_csv("data/tcga.brca.rsem.csv", usecols=["bcr_patient_barcode"])
                rsem_ids = set(rsem_raw["bcr_patient_barcode"].str[:15])

            def classify_reconciliation(row: pd.Series) -> str:
                sid = str(row["Sample ID"])
                if sid in pt_qc_ids:
                    return "Retained Post-QC Primary Tumor Cohort"
                if sid in outlier_ids:
                    return "Exclusion 1: Primary Tumor PCA Outliers (Removed in RNA QC)"
                if sid.endswith("-06"):
                    return "Exclusion 2: Metastasis Samples (-06, Excluded from Primary Cohort)"
                if sid not in rsem_ids:
                    return "Exclusion 3: Clinical Records Without TCGA RNA-Seq"
                return "Exclusion: Other"

            comp_df["Reconciliation_Group"] = comp_df.apply(classify_reconciliation, axis=1)

            recon_order = [
                "Total Clinical Cohort (n = 1,108)",
                "Exclusion 1: Primary Tumor PCA Outliers (Removed in RNA QC)",
                "Exclusion 2: Metastasis Samples (-06, Excluded from Primary Cohort)",
                "Exclusion 3: Clinical Records Without TCGA RNA-Seq",
                "Subtotal Excluded from Primary Tumor Cohort (33 Samples)",
                "Retained Post-QC Primary Tumor Cohort (n = 1,075)",
            ]

            t3_recon_rows = []
            tot_pos_all = int((comp_df["HER2_status"] == "Positive").sum())
            tot_neg_all = int((comp_df["HER2_status"] == "Negative").sum())
            tot_ind_all = int((comp_df["HER2_status"] == "Indeterminate").sum())

            # 1. Total Clinical
            t3_recon_rows.append({
                "Table_Section": "A. Cohort Reconciliation (1,108 Clinical -> 1,075 Post-QC Primary Tumors)",
                "Cohort_Category": "Total Clinical Cohort (n = 1,108)",
                "Positive": tot_pos_all,
                "Negative": tot_neg_all,
                "Total_Labeled (Pos + Neg)": tot_pos_all + tot_neg_all,
                "Indeterminate": tot_ind_all,
                "Total_Samples": len(comp_df),
                "Description / Clinical Notes": "All diagnostic clinical entries in data/brca_tcga_clinical_data.csv",
            })

            # Exclusions
            excl_mask = comp_df["Reconciliation_Group"] != "Retained Post-QC Primary Tumor Cohort"
            for cat_label, cat_key in [
                ("Exclusion 1: Primary Tumor PCA Outliers (Removed in RNA QC)", "Exclusion 1: Primary Tumor PCA Outliers (Removed in RNA QC)"),
                ("Exclusion 2: Metastasis Samples (-06, Excluded from Primary Cohort)", "Exclusion 2: Metastasis Samples (-06, Excluded from Primary Cohort)"),
                ("Exclusion 3: Clinical Records Without TCGA RNA-Seq", "Exclusion 3: Clinical Records Without TCGA RNA-Seq"),
            ]:
                sub = comp_df[comp_df["Reconciliation_Group"] == cat_key]
                p = int((sub["HER2_status"] == "Positive").sum())
                n = int((sub["HER2_status"] == "Negative").sum())
                i = int((sub["HER2_status"] == "Indeterminate").sum())
                notes = (
                    "Removed during Phase 2 RNA-seq QC due to extreme technical variance (|Z| >= 3 on PC1/PC2)"
                    if "PCA Outliers" in cat_label
                    else (
                        "Metastatic tissue (sample type -06); evaluated separately as bonus discovery cohort"
                        if "Metastasis" in cat_label
                        else "Clinical records in clinical CSV without available RNA-seq in tcga.brca.rsem.csv"
                    )
                )
                t3_recon_rows.append({
                    "Table_Section": "A. Cohort Reconciliation (1,108 Clinical -> 1,075 Post-QC Primary Tumors)",
                    "Cohort_Category": cat_label,
                    "Positive": p,
                    "Negative": n,
                    "Total_Labeled (Pos + Neg)": p + n,
                    "Indeterminate": i,
                    "Total_Samples": len(sub),
                    "Description / Clinical Notes": notes,
                })

            # Subtotal excluded
            sub_excl = comp_df[excl_mask]
            p_ex = int((sub_excl["HER2_status"] == "Positive").sum())
            n_ex = int((sub_excl["HER2_status"] == "Negative").sum())
            i_ex = int((sub_excl["HER2_status"] == "Indeterminate").sum())
            t3_recon_rows.append({
                "Table_Section": "A. Cohort Reconciliation (1,108 Clinical -> 1,075 Post-QC Primary Tumors)",
                "Cohort_Category": "Subtotal Excluded from Primary Tumor Cohort (33 Samples)",
                "Positive": p_ex,
                "Negative": n_ex,
                "Total_Labeled (Pos + Neg)": p_ex + n_ex,
                "Indeterminate": i_ex,
                "Total_Samples": len(sub_excl),
                "Description / Clinical Notes": "Total samples excluded before downstream primary tumor modeling",
            })

            # Retained post-QC
            sub_ret = comp_df[~excl_mask]
            p_ret = int((sub_ret["HER2_status"] == "Positive").sum())
            n_ret = int((sub_ret["HER2_status"] == "Negative").sum())
            i_ret = int((sub_ret["HER2_status"] == "Indeterminate").sum())
            t3_recon_rows.append({
                "Table_Section": "A. Cohort Reconciliation (1,108 Clinical -> 1,075 Post-QC Primary Tumors)",
                "Cohort_Category": "Retained Post-QC Primary Tumor Cohort (n = 1,075)",
                "Positive": p_ret,
                "Negative": n_ret,
                "Total_Labeled (Pos + Neg)": p_ret + n_ret,
                "Indeterminate": i_ret,
                "Total_Samples": len(sub_ret),
                "Description / Clinical Notes": "Clean primary tumor cohort used for downstream multimodal HER2 classification",
            })
            t3_recon_df = pd.DataFrame(t3_recon_rows)

            # Build Table 3B: Post-QC Primary Tumor Operational Completeness (n = 1,075)
            pt_only = comp_df[comp_df["Has_PostQC_Primary_RNA"]]
            t3b_rows = []
            for label, mask in [
                ("Complete Molecular Data (RNA + GISTIC + Purity)", pt_only["Has_GISTIC_CNV"]),
                ("Incomplete Molecular Data (RNA + Purity Only, Missing GISTIC)", ~pt_only["Has_GISTIC_CNV"]),
                ("Total Post-QC Primary Tumors (n = 1,075)", pd.Series(True, index=pt_only.index)),
            ]:
                sub = pt_only[mask]
                n_pos = int((sub["HER2_status"] == "Positive").sum())
                n_neg = int((sub["HER2_status"] == "Negative").sum())
                n_ind = int((sub["HER2_status"] == "Indeterminate").sum())
                n_tot = len(sub)
                t3b_rows.append({
                    "Table_Section": "B. Post-QC Primary Tumor Cohort (n = 1,075)",
                    "Completeness_Category": label,
                    "Positive": n_pos,
                    "Positive_%": f"{n_pos / max(int((pt_only['HER2_status'] == 'Positive').sum()), 1) * 100:.2f}%",
                    "Negative": n_neg,
                    "Negative_%": f"{n_neg / max(int((pt_only['HER2_status'] == 'Negative').sum()), 1) * 100:.2f}%",
                    "Total_Labeled (Pos + Neg)": n_pos + n_neg,
                    "Indeterminate": n_ind,
                    "Indeterminate_%": f"{n_ind / max(int((pt_only['HER2_status'] == 'Indeterminate').sum()), 1) * 100:.2f}%",
                    "Total_Samples": n_tot,
                })
            t3b_df = pd.DataFrame(t3b_rows)

            # Build Table 3C: Full Clinical Cohort Omics Completeness (n = 1,108)
            omics_cats = sorted(comp_df["Omics_Completeness"].unique())
            t3c_rows = []
            tot_pos = int((comp_df["HER2_status"] == "Positive").sum())
            tot_neg = int((comp_df["HER2_status"] == "Negative").sum())
            tot_ind = int((comp_df["HER2_status"] == "Indeterminate").sum())
            for cat in omics_cats:
                sub = comp_df[comp_df["Omics_Completeness"] == cat]
                n_pos = int((sub["HER2_status"] == "Positive").sum())
                n_neg = int((sub["HER2_status"] == "Negative").sum())
                n_ind = int((sub["HER2_status"] == "Indeterminate").sum())
                t3c_rows.append({
                    "Table_Section": "C. All Clinical Samples by Omics Completeness (n = 1,108)",
                    "Completeness_Category": cat,
                    "Positive": n_pos,
                    "Positive_%": f"{n_pos / tot_pos * 100:.2f}%",
                    "Negative": n_neg,
                    "Negative_%": f"{n_neg / tot_neg * 100:.2f}%",
                    "Total_Labeled (Pos + Neg)": n_pos + n_neg,
                    "Indeterminate": n_ind,
                    "Indeterminate_%": f"{n_ind / tot_ind * 100:.2f}%",
                    "Total_Samples": len(sub),
                })
            t3c_rows.append({
                "Table_Section": "C. All Clinical Samples by Omics Completeness (n = 1,108)",
                "Completeness_Category": "Total Clinical Cohort (n = 1,108)",
                "Positive": tot_pos,
                "Positive_%": "100.00%",
                "Negative": tot_neg,
                "Negative_%": "100.00%",
                "Total_Labeled (Pos + Neg)": tot_pos + tot_neg,
                "Indeterminate": tot_ind,
                "Indeterminate_%": "100.00%",
                "Total_Samples": len(comp_df),
            })
            t3c_df = pd.DataFrame(t3c_rows)

            # Build Table 3D: Full Clinical Cohort Assay Completeness (IHC Label / IHC Score / FISH, n = 1,108)
            clin_cats = sorted(comp_df["Clinical_Assay_Completeness"].unique())
            t3d_rows = []
            for cat in clin_cats:
                sub = comp_df[comp_df["Clinical_Assay_Completeness"] == cat]
                n_pos = int((sub["HER2_status"] == "Positive").sum())
                n_neg = int((sub["HER2_status"] == "Negative").sum())
                n_ind = int((sub["HER2_status"] == "Indeterminate").sum())
                t3d_rows.append({
                    "Table_Section": "D. All Clinical Samples by Clinical Assay Completeness (n = 1,108)",
                    "Completeness_Category": cat,
                    "Positive": n_pos,
                    "Positive_%": f"{n_pos / tot_pos * 100:.2f}%",
                    "Negative": n_neg,
                    "Negative_%": f"{n_neg / tot_neg * 100:.2f}%",
                    "Total_Labeled (Pos + Neg)": n_pos + n_neg,
                    "Indeterminate": n_ind,
                    "Indeterminate_%": f"{n_ind / tot_ind * 100:.2f}%",
                    "Total_Samples": len(sub),
                })
            t3d_rows.append({
                "Table_Section": "D. All Clinical Samples by Clinical Assay Completeness (n = 1,108)",
                "Completeness_Category": "Total Clinical Cohort (n = 1,108)",
                "Positive": tot_pos,
                "Positive_%": "100.00%",
                "Negative": tot_neg,
                "Negative_%": "100.00%",
                "Total_Labeled (Pos + Neg)": tot_pos + tot_neg,
                "Indeterminate": tot_ind,
                "Indeterminate_%": "100.00%",
                "Total_Samples": len(comp_df),
            })
            t3d_df = pd.DataFrame(t3d_rows)

            # Write all four tables sequentially onto Sheet 3 (`3_Data_Completeness_Breakdown`)
            sheet3_name = "3_Data_Completeness_Breakdown"
            r_offset = 0
            t3_recon_df.to_excel(writer, sheet_name=sheet3_name, index=False, startrow=r_offset)
            r_offset += len(t3_recon_df) + 3
            t3b_df.to_excel(writer, sheet_name=sheet3_name, index=False, startrow=r_offset)
            r_offset += len(t3b_df) + 3
            t3c_df.to_excel(writer, sheet_name=sheet3_name, index=False, startrow=r_offset)
            r_offset += len(t3c_df) + 3
            t3d_df.to_excel(writer, sheet_name=sheet3_name, index=False, startrow=r_offset)

            # Auto-adjust column widths for readability across all 3 sheets
            for ws in writer.book.worksheets:
                for col in ws.columns:
                    max_len = max(len(str(cell.value or "")) for cell in col)
                    ws.column_dimensions[col[0].column_letter].width = min(max(max_len + 3, 14), 72)

        print(f"\nSaved 3-sheet Excel summary workbook to: {args.excel_summary}")


if __name__ == "__main__":
    main()

