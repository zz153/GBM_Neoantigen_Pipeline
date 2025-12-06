#!/usr/bin/env Rscript
################################################################################
# 02_filter_mutations.R
# 
# Filters mutations for IDH-wildtype patients, calculates VAF, and identifies
# top mutated genes
# 
# Inputs:
#   - data/raw/TCGA_mutation_maf.rds
#   - data/raw/idhwt_patient_ids.txt
# 
# Outputs:
#   - data/processed/TCGA_mutation_maf_IDHwt.csv
#   - data/processed/Top50_Mutated_Genes.csv
#   - results/tables/Table1_Mutation_Summary.csv
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("STEP 2: FILTERING MUTATIONS FOR IDH-WILDTYPE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────────────────────────

library(dplyr)

# Create output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

cat("Loading data...\n")

# ─────────────────────────────────────────────────────────────────
# 1. Load Data
# ─────────────────────────────────────────────────────────────────

cat("[1] Loading mutation data...\n")

# Check if file exists
if (!file.exists("data/raw/TCGA_mutation_maf.rds")) {
  stop("ERROR: data/raw/TCGA_mutation_maf.rds not found.\n",
       "       Please run scripts/01_download_data.R first.")
}

maf <- readRDS("data/raw/TCGA_mutation_maf.rds")
cat("  ✓ Loaded", nrow(maf), "total mutations\n")

# Load IDH-WT patient IDs
if (!file.exists("data/raw/idhwt_patient_ids.txt")) {
  stop("ERROR: data/raw/idhwt_patient_ids.txt not found.\n",
       "       Please run scripts/01_download_data.R first.")
}

idhwt_patient_ids <- readLines("data/raw/idhwt_patient_ids.txt")
cat("  ✓ Loaded", length(idhwt_patient_ids), "IDH-WT patient IDs\n\n")

# ─────────────────────────────────────────────────────────────────
# 2. Filter for IDH-Wildtype Patients
# ─────────────────────────────────────────────────────────────────

cat("[2] Filtering for IDH-wildtype mutations...\n")

# Extract patient IDs from tumor sample barcodes (first 12 characters)
maf$patient_id <- substr(maf$Tumor_Sample_Barcode, 1, 12)

# Filter for IDH-WT patients
maf_idhwt <- maf[maf$patient_id %in% idhwt_patient_ids, ]

n_idhwt_muts <- nrow(maf_idhwt)
n_idhwt_patients <- length(unique(maf_idhwt$patient_id))

cat("  ✓ IDH-WT mutations:", n_idhwt_muts, "\n")
cat("  ✓ From", n_idhwt_patients, "patients\n")
cat("  ✓ Filtered out", nrow(maf) - n_idhwt_muts, "non-IDH-WT mutations\n\n")

# ─────────────────────────────────────────────────────────────────
# 3. Calculate Variant Allele Frequency (VAF)
# ─────────────────────────────────────────────────────────────────

cat("[3] Calculating Variant Allele Frequency (VAF)...\n")

# Check for required columns
vaf_cols <- c("t_alt_count", "t_depth")
missing_cols <- vaf_cols[!vaf_cols %in% colnames(maf_idhwt)]

if (length(missing_cols) > 0) {
  stop("ERROR: Missing columns for VAF calculation: ", 
       paste(missing_cols, collapse = ", "))
}

# Calculate VAF
maf_idhwt$VAF <- maf_idhwt$t_alt_count / maf_idhwt$t_depth

# Remove mutations with NA or invalid VAF
n_before <- nrow(maf_idhwt)
maf_idhwt <- maf_idhwt[!is.na(maf_idhwt$VAF) & 
                         maf_idhwt$VAF >= 0 & 
                         maf_idhwt$VAF <= 1, ]
n_after <- nrow(maf_idhwt)

cat("  ✓ VAF calculated for", n_after, "mutations\n")
cat("  ✓ Removed", n_before - n_after, "mutations with invalid VAF\n")
cat("  ✓ VAF range:", round(min(maf_idhwt$VAF, na.rm = TRUE), 3), 
    "to", round(max(maf_idhwt$VAF, na.rm = TRUE), 3), "\n")
cat("  ✓ Mean VAF:", round(mean(maf_idhwt$VAF, na.rm = TRUE), 3), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 4. Save Filtered Mutations
# ─────────────────────────────────────────────────────────────────

cat("[4] Saving filtered mutations...\n")

# Select key columns for clean output
key_cols <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position",
              "Variant_Classification", "Variant_Type", "Reference_Allele",
              "Tumor_Seq_Allele2", "HGVSp_Short", "Tumor_Sample_Barcode",
              "patient_id", "t_depth", "t_alt_count", "VAF")

# Keep only columns that exist
available_cols <- key_cols[key_cols %in% colnames(maf_idhwt)]
maf_export <- maf_idhwt[, available_cols]

write.csv(maf_export, 
          "data/processed/TCGA_mutation_maf_IDHwt.csv", 
          row.names = FALSE)

cat("  ✓ Saved:", "data/processed/TCGA_mutation_maf_IDHwt.csv\n")
cat("  ✓ Columns:", ncol(maf_export), "\n")
cat("  ✓ Rows:", nrow(maf_export), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 5. Identify Top Mutated Genes
# ─────────────────────────────────────────────────────────────────

cat("[5] Identifying top 50 mutated genes...\n")

# Count mutations per gene
gene_counts <- sort(table(maf_idhwt$Hugo_Symbol), decreasing = TRUE)

# Create summary dataframe
top_genes_df <- data.frame(
  Gene = names(gene_counts),
  N_mutations = as.numeric(gene_counts),
  stringsAsFactors = FALSE
)

# Count unique patients per gene
top_genes_df$N_patients <- sapply(top_genes_df$Gene, function(g) {
  length(unique(maf_idhwt$patient_id[maf_idhwt$Hugo_Symbol == g]))
})

# Calculate percentage of patients
top_genes_df$Percent_patients <- round(
  100 * top_genes_df$N_patients / n_idhwt_patients, 1
)

# Calculate mean VAF per gene
top_genes_df$Mean_VAF <- sapply(top_genes_df$Gene, function(g) {
  round(mean(maf_idhwt$VAF[maf_idhwt$Hugo_Symbol == g], na.rm = TRUE), 3)
})

# Save top 50
top50 <- head(top_genes_df, 50)
write.csv(top50, 
          "data/processed/Top50_Mutated_Genes.csv", 
          row.names = FALSE)

cat("  ✓ Saved: data/processed/Top50_Mutated_Genes.csv\n\n")

# ─────────────────────────────────────────────────────────────────
# 6. Generate Summary Table
# ─────────────────────────────────────────────────────────────────

cat("[6] Generating mutation summary...\n")

# Overall statistics
summary_stats <- data.frame(
  Metric = c(
    "Total patients (IDH-WT)",
    "Total mutations",
    "Mutations per patient (mean)",
    "Mutations per patient (median)",
    "Genes mutated (total)",
    "Mean VAF",
    "Median VAF"
  ),
  Value = c(
    n_idhwt_patients,
    n_idhwt_muts,
    round(n_idhwt_muts / n_idhwt_patients, 1),
    median(table(maf_idhwt$patient_id)),
    length(unique(maf_idhwt$Hugo_Symbol)),
    round(mean(maf_idhwt$VAF, na.rm = TRUE), 3),
    round(median(maf_idhwt$VAF, na.rm = TRUE), 3)
  )
)

write.csv(summary_stats, 
          "results/tables/Table1_Mutation_Summary.csv", 
          row.names = FALSE)

cat("  ✓ Saved: results/tables/Table1_Mutation_Summary.csv\n\n")

# Print summary
cat("Summary Statistics:\n")
print(summary_stats, row.names = FALSE)

cat("\n")

# ─────────────────────────────────────────────────────────────────
# 7. Display Top 20 Genes
# ─────────────────────────────────────────────────────────────────

cat("[7] Top 20 most frequently mutated genes:\n\n")

top20_display <- head(top_genes_df, 20)
print(top20_display, row.names = FALSE)

cat("\n")

# ─────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat("MUTATION FILTERING COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Files created:\n")
cat("  • data/processed/TCGA_mutation_maf_IDHwt.csv\n")
cat("  • data/processed/Top50_Mutated_Genes.csv\n")
cat("  • results/tables/Table1_Mutation_Summary.csv\n\n")

cat("Next step: Run scripts/03_expression_validation.R\n\n")

cat("═══════════════════════════════════════════════════════════════\n\n")
