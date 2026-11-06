#!/usr/bin/env Rscript
# ==============================================================================
# Adaptive RNA-Seq Diagnostic, QC, Prevalence Filtering, Normalization &
# PCA-Based Sample Outlier Detection Pipeline
# Cohort: TCGA-BRCA (data/tcga.brca.rsem.csv)
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

# ------------------------------------------------------------------------------
# Step 1: Data Parsing & Matrix Separation
# ------------------------------------------------------------------------------
cat("=== Step 1: Data Parsing & Matrix Separation ===\n")

input_path <- "data/tcga.brca.rsem.csv"
if (!file.exists(input_path) && file.exists("data/tcga.brca.rsem.csv.gz")) {
  input_path <- "data/tcga.brca.rsem.csv.gz"
}
output_dir <- "results"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Fast load of the wide RSEM matrix (rows = samples, cols = 3 metadata + genes)
raw_dt <- fread(input_path, data.table = FALSE)

# Extract the 3 metadata columns
meta_cols <- c("patient_id", "bcr_patient_barcode", "sample_type")
metadata <- raw_dt[, meta_cols]

# Derive 15-character TCGA Sample ID (e.g., "TCGA-3C-AAAU-01") from bcr_patient_barcode.
# Note: patient_id (12-char) is NOT unique because patients can have Primary Tumor (-01),
# Adjacent Normal (-11), and Metastasis (-06) samples. Both the 28-char bcr_patient_barcode
# and its 15-char sample_id prefix are 100% unique across all 1,212 rows and match the
# clinical and copy-number tables directly.
metadata$sample_id <- substr(metadata$bcr_patient_barcode, 1, 15)

if (anyDuplicated(metadata$sample_id) > 0) {
  stop("Duplicate sample_id values detected! Check barcode formatting.")
}

metadata <- metadata[, c("sample_id", "patient_id", "bcr_patient_barcode", "sample_type")]
rownames(metadata) <- metadata$sample_id
metadata$sample_type <- factor(metadata$sample_type)

# Extract gene columns into a distinct numeric matrix and transpose:
# Rows = Genes (20,502), Columns = Samples (1,212)
gene_matrix_wide <- as.matrix(raw_dt[, setdiff(colnames(raw_dt), meta_cols)])
storage.mode(gene_matrix_wide) <- "double"

counts_matrix <- t(gene_matrix_wide)
colnames(counts_matrix) <- rownames(metadata)
rm(raw_dt, gene_matrix_wide)

cat(sprintf(
  "Parsed counts_matrix: %d genes (rows) x %d samples (columns)\n\n",
  nrow(counts_matrix), ncol(counts_matrix)
))

# ------------------------------------------------------------------------------
# Step 2: Blind Data Profiling & Diagnostic Logic Gates
# ------------------------------------------------------------------------------
cat("=== Step 2: Blind Data Profiling & Diagnostic Logic Gates ===\n")

# 1. Column Sums: total sum of expression values per sample column
col_sums <- colSums(counts_matrix, na.rm = TRUE)

# 2. Max Value: absolute maximum value across the entire matrix
max_val <- max(counts_matrix, na.rm = TRUE)

# 3. Quantiles: 75th percentile (upper quartile) of non-zero genes for each sample
q75_nonzero <- apply(counts_matrix, 2, function(x) {
  quantile(x[x > 0], probs = 0.75, na.rm = TRUE)
})
q75_cv <- sd(q75_nonzero) / mean(q75_nonzero)

cat("1. Sample Column Sums Summary:\n")
print(summary(col_sums))
cat(sprintf("\n2. Matrix Absolute Maximum Value: %.4f\n\n", max_val))
cat("3. 75th Percentile of Non-Zero Genes per Sample Summary:\n")
print(summary(q75_nonzero))
cat(sprintf("   75th Percentile Coefficient of Variation (CV): %.6f\n\n", q75_cv))

# Evaluate strict diagnostic logic gates:
if (max_val < 25) {
  # GATE A (Log Scaled)
  stop(
    "GATE A TRIGGERED: Maximum value in matrix is < 25 (", round(max_val, 3), "). ",
    "Data is already log-transformed. Halting pipeline: raw count or linear scaling is required."
  )
} else if (all(abs(col_sums - 1e6) < 1e-2)) {
  # GATE B (TPM / CPM Scaled)
  warning(
    "GATE B TRIGGERED: Every sample column sum equals 1,000,000 (1e6). ",
    "Data is pre-normalized to TPM/CPM; DESeq2 cannot be used. Proceeding with 'limma-voom'."
  )
  pipeline_route <- "limma-voom"
} else if (all(q75_nonzero >= 995 & q75_nonzero <= 1005) && q75_cv < 0.01) {
  # GATE C (Upper-Quartile Scaled)
  cat(
    "GATE C TRIGGERED: Column sums vary widely (",
    sprintf("%.2fM to %.2fM", min(col_sums) / 1e6, max(col_sums) / 1e6),
    "), but the 75th percentile of non-zero genes scales uniformly to ~1,000 across all samples ",
    sprintf("(range = %.2f - %.2f, median = %.2f, CV = %.4f%%).\n",
            min(q75_nonzero), max(q75_nonzero), median(q75_nonzero), q75_cv * 100),
    "Data has been Upper-Quartile normalized (TCGA Level 3 RSEM). Proceeding with 'limma-voom' pipeline.\n\n",
    sep = ""
  )
  pipeline_route <- "limma-voom"
} else {
  # GATE D (Raw / Expected Counts)
  cat(
    "GATE D TRIGGERED: Column sums vary widely, max value is large, and 75th percentile ",
    "varies dynamically across samples. Proceeding with standard 'DESeq2' pipeline.\n\n"
  )
  pipeline_route <- "DESeq2"
}

# ------------------------------------------------------------------------------
# Step 3: Implement Stringent Prevalence Filtering
# ------------------------------------------------------------------------------
cat("=== Step 3: Stringent Prevalence Filtering ===\n")

# Retain only genes with expression value >= 15 in at least 200 samples
min_expr <- 15
min_samples <- 200
keep_genes <- rowSums(counts_matrix >= min_expr) >= min_samples
filtered_matrix <- counts_matrix[keep_genes, , drop = FALSE]

cat(sprintf(
  "Retained %d / %d genes (%.1f%%) with expression >= %d in >= %d samples.\n",
  nrow(filtered_matrix), nrow(counts_matrix),
  100 * nrow(filtered_matrix) / nrow(counts_matrix),
  min_expr, min_samples
))
cat(sprintf("ERBB2 retained in filtered_matrix: %s\n\n", "ERBB2" %in% rownames(filtered_matrix)))

# ------------------------------------------------------------------------------
# Step 4: Framework Processing & Normalization (Helper & Initial Pass)
# ------------------------------------------------------------------------------
run_normalization <- function(mat, meta, route, pass_label = "Initial") {
  cat(sprintf("--- Running %s Normalization (%s, n = %d samples) ---\n",
              pass_label, route, ncol(mat)))

  if (route == "limma-voom") {
    suppressPackageStartupMessages({
      library(edgeR)
      library(limma)
    })

    # Initialize edgeR DGEList object
    dge <- DGEList(counts = mat, samples = meta)

    # Calculate TMM normalization factors to correct residual compositional differences
    dge <- calcNormFactors(dge, method = "TMM")
    cat("TMM Normalization Factors Summary:\n")
    print(summary(dge$samples$norm.factors))

    # Matrix 1 (Linear): divide counts columns by effective library sizes
    # and scale back to the average effective library size
    eff_lib_sizes <- dge$samples$lib.size * dge$samples$norm.factors
    avg_eff_lib_size <- mean(eff_lib_sizes)
    depth_scale_factors <- eff_lib_sizes / avg_eff_lib_size
    norm_linear <- sweep(mat, 2, depth_scale_factors, FUN = "/")

    # Matrix 2 (Log2 Scale): stabilize variance using limma::voom
    design_matrix <- model.matrix(~ sample_type, data = meta)
    v <- voom(dge, design = design_matrix, plot = FALSE)
    norm_log2 <- v$E

  } else if (route == "DESeq2") {
    suppressPackageStartupMessages({
      library(DESeq2)
    })

    # Initialize DESeq2 object with rounded integer counts
    dds <- DESeqDataSetFromMatrix(
      countData = round(mat),
      colData = meta,
      design = ~ sample_type
    )

    # Calculate sequencing depth normalization factors
    dds <- estimateSizeFactors(dds)
    cat("DESeq2 Size Factors Summary:\n")
    print(summary(sizeFactors(dds)))

    # Matrix 1 (Linear): sequencing depth-corrected counts on a linear scale
    norm_linear <- counts(dds, normalized = TRUE)

    # Matrix 2 (Log2 Scale): Variance Stabilizing Transformation
    vsd <- vst(dds, blind = FALSE)
    norm_log2 <- assay(vsd)
  }

  list(linear = norm_linear, log2 = norm_log2)
}

cat("=== Step 4: Initial Normalization Pass (Pre-Outlier QC) ===\n")
initial_norm <- run_normalization(filtered_matrix, metadata, pipeline_route, pass_label = "Initial")
cat("\n")

# ------------------------------------------------------------------------------
# Step 5: PCA-Based Sample Outlier Detection (MAD |Z| >= 3) & Visualization
# ------------------------------------------------------------------------------
cat("=== Step 5: PCA-Based Sample Outlier Detection (MAD |Z| >= 3) ===\n")

# Run PCA on the initial log2-stabilized matrix (rows = samples, cols = genes).
# Note: scale. = TRUE weights all 15,632 genes equally so PC1/PC2 capture global
# sample-level technical variation rather than single hyper-variable lineage genes
# (e.g., ESR1, which otherwise causes PC2 to flag 145 Basal-like/ER-negative tumors).
pca_before <- prcomp(t(initial_norm$log2), center = TRUE, scale. = TRUE)
var_exp_before <- (pca_before$sdev^2) / sum(pca_before$sdev^2) * 100

pc1_before <- pca_before$x[, "PC1"]
pc2_before <- pca_before$x[, "PC2"]

# Calculate robust Z-scores using R's mad() function: Z = (PC_i - median(PC_i)) / mad(PC_i)
z_pc1 <- (pc1_before - median(pc1_before)) / mad(pc1_before)
z_pc2 <- (pc2_before - median(pc2_before)) / mad(pc2_before)

# Flag samples with |Z| >= 3 on either PC1 or PC2
is_outlier <- (abs(z_pc1) >= 3) | (abs(z_pc2) >= 3)
outlier_ids <- rownames(metadata)[is_outlier]

cat(sprintf(
  "Detected %d outlier samples (|Z_PC1| >= 3: %d, |Z_PC2| >= 3: %d) out of %d total samples.\n",
  sum(is_outlier), sum(abs(z_pc1) >= 3), sum(abs(z_pc2) >= 3), nrow(metadata)
))
cat("Flagged Outlier Sample IDs:\n")
print(outlier_ids)
cat("\nOutlier Breakdown by Sample Type:\n")
print(table(metadata$sample_type[is_outlier]))
cat("\n")

# Save outlier tracking table for auditability
pca_before_df <- data.frame(
  sample_id = rownames(metadata),
  patient_id = metadata$patient_id,
  sample_type = metadata$sample_type,
  PC1 = pc1_before,
  PC2 = pc2_before,
  Z_PC1 = z_pc1,
  Z_PC2 = z_pc2,
  is_outlier = is_outlier,
  status = ifelse(is_outlier, "Flagged Outlier (|Z| >= 3)", "Retained (|Z| < 3)"),
  stringsAsFactors = FALSE
)

outliers_df <- pca_before_df[pca_before_df$is_outlier, ]
outliers_df <- outliers_df[order(-pmax(abs(outliers_df$Z_PC1), abs(outliers_df$Z_PC2))), ]
fwrite(outliers_df, file.path(output_dir, "pca_outliers_removed.csv"), row.names = FALSE)

# Calculate PC1 and PC2 cutoff coordinates corresponding to |Z| = 3
pc1_cutoffs <- median(pc1_before) + c(-3, 3) * mad(pc1_before)
pc2_cutoffs <- median(pc2_before) + c(-3, 3) * mad(pc2_before)

# Export "Before" PCA plot
p_before <- ggplot(pca_before_df, aes(x = PC1, y = PC2, color = status)) +
  geom_vline(xintercept = pc1_cutoffs, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = pc2_cutoffs, linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_manual(
    values = c(
      "Retained (|Z| < 3)" = "steelblue",
      "Flagged Outlier (|Z| >= 3)" = "#D62728"
    )
  ) +
  labs(
    title = "PCA Before Outlier Removal (Initial voom Log2 Matrix)",
    subtitle = sprintf(
      "n = %d samples | %d flagged outliers with MAD |Z| >= 3 on PC1 or PC2 (dashed lines)",
      nrow(pca_before_df), sum(is_outlier)
    ),
    x = sprintf("PC1 (%.1f%% variance)", var_exp_before[1]),
    y = sprintf("PC2 (%.1f%% variance)", var_exp_before[2]),
    color = "QC Status"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(output_dir, "pca_before_outlier_removal.png"),
  plot = p_before, width = 8, height = 6, dpi = 300
)
cat("Saved results/pca_before_outlier_removal.png\n\n")

# Programmatically remove flagged outlier samples from metadata and filtered count matrix
clean_metadata <- metadata[!is_outlier, , drop = FALSE]
clean_Metadata_levels <- droplevels(clean_metadata$sample_type)
clean_metadata$sample_type <- clean_Metadata_levels
clean_filtered_matrix <- filtered_matrix[, rownames(clean_metadata), drop = FALSE]

# Re-run Step 4 normalization on the clean, outlier-removed cohort
final_norm <- run_normalization(
  clean_filtered_matrix, clean_metadata, pipeline_route, pass_label = "Final Post-QC"
)

# Run PCA on the final, clean post-normalization log2 matrix
pca_after <- prcomp(t(final_norm$log2), center = TRUE, scale. = TRUE)
var_exp_after <- (pca_after$sdev^2) / sum(pca_after$sdev^2) * 100

pca_after_df <- data.frame(
  sample_id = rownames(clean_metadata),
  sample_type = clean_metadata$sample_type,
  PC1 = pca_after$x[, "PC1"],
  PC2 = pca_after$x[, "PC2"],
  stringsAsFactors = FALSE
)

p_after <- ggplot(pca_after_df, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_manual(
    values = c(
      "Primary Tumor" = "#1F77B4",
      "Adjacent Normal" = "#2CA02C",
      "Metastasis" = "#FF7F0E"
    )
  ) +
  labs(
    title = "PCA After Outlier Removal (Final Clean voom Log2 Matrix)",
    subtitle = sprintf(
      "n = %d retained samples (%d Primary Tumor, %d Adjacent Normal, %d Metastasis)",
      nrow(pca_after_df),
      sum(clean_metadata$sample_type == "Primary Tumor"),
      sum(clean_metadata$sample_type == "Adjacent Normal"),
      sum(clean_metadata$sample_type == "Metastasis")
    ),
    x = sprintf("PC1 (%.1f%% variance)", var_exp_after[1]),
    y = sprintf("PC2 (%.1f%% variance)", var_exp_after[2]),
    color = "Sample Type"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(output_dir, "pca_after_outlier_removal.png"),
  plot = p_after, width = 8, height = 6, dpi = 300
)
cat("Saved results/pca_after_outlier_removal.png\n\n")

# ------------------------------------------------------------------------------
# Step 6: Output Export
# ------------------------------------------------------------------------------
cat("=== Step 6: Exporting Final Cleaned Matrices & Metadata ===\n")

meta_out_path   <- file.path(output_dir, "rsem_sample_metadata.csv")
linear_out_path <- file.path(output_dir, "matrix_normalized_linear.csv")
log2_out_path   <- file.path(output_dir, "matrix_normalized_log2_stabilized.csv")

cat(sprintf("Saving cleaned metadata (%d samples) to %s ...\n", nrow(clean_metadata), meta_out_path))
fwrite(clean_metadata, meta_out_path, row.names = FALSE)

cat(sprintf("Saving final linear normalized matrix to %s ...\n", linear_out_path))
linear_dt <- as.data.table(final_norm$linear, keep.rownames = "gene_symbol")
fwrite(linear_dt, linear_out_path)

cat(sprintf("Saving final log2 stabilized matrix to %s ...\n", log2_out_path))
log2_dt <- as.data.table(final_norm$log2, keep.rownames = "gene_symbol")
fwrite(log2_dt, log2_out_path)

cat("\n=== Pipeline Completed Successfully ===\n")
