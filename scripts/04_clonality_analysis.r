#!/usr/bin/env Rscript
################################################################################
# 04_clonality_analysis.R
# 
# Identifies clonal mutations (VAF > 0.4) in expressed genes
# 
# Inputs:
#   - data/processed/TCGA_mutation_maf_IDHwt.csv
#   - data/processed/expression_summary.csv
# 
# Outputs:
#   - data/processed/clonal_mutations.csv
#   - data/processed/Evolution_Resistant_Candidates.csv
#   - results/tables/Table2_Clonality_Summary.csv
#   - results/figures/Figure2_VAF_Distribution.pdf
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("STEP 4: CLONALITY ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────────────────────────

library(dplyr)
library(ggplot2)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────
# 1. Load Data
# ─────────────────────────────────────────────────────────────────

cat("[1] Loading mutation and expression data...\n")

maf_idhwt <- read.csv("data/processed/TCGA_mutation_maf_IDHwt.csv")
expr_summary <- read.csv("data/processed/expression_summary.csv")

cat("  ✓ Mutations:", nrow(maf_idhwt), "\n")
cat("  ✓ Expression data loaded\n\n")

# ─────────────────────────────────────────────────────────────────
# 2. Identify Expressed Genes (TPM > 10 in tumor)
# ─────────────────────────────────────────────────────────────────

cat("[2] Filtering for expressed genes...\n")

# Get tumor expression only
tumor_expr <- expr_summary[expr_summary$Type == "TCGA_Tumor", ]
expressed_genes <- tumor_expr$Gene[tumor_expr$Mean_TPM > 10]

cat("  ✓ Genes with TPM > 10:", length(expressed_genes), "\n")
cat("  ✓ Genes:", paste(expressed_genes, collapse = ", "), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 3. Filter Mutations for Expressed Genes
# ─────────────────────────────────────────────────────────────────

cat("[3] Filtering mutations in expressed genes...\n")

maf_expressed <- maf_idhwt[maf_idhwt$Hugo_Symbol %in% expressed_genes, ]

cat("  ✓ Mutations in expressed genes:", nrow(maf_expressed), "\n")
cat("  ✓ From", nrow(maf_idhwt), "total mutations\n\n")

# ─────────────────────────────────────────────────────────────────
# 4. Calculate Clonality (VAF-based)
# ─────────────────────────────────────────────────────────────────

cat("[4] Calculating clonality metrics...\n")

# Classify by VAF
maf_expressed$clonality <- ifelse(maf_expressed$VAF > 0.4, "Clonal",
                                  ifelse(maf_expressed$VAF > 0.2, "Subclonal", "Rare"))

# Summary by gene
vaf_summary <- maf_expressed %>%
  group_by(Hugo_Symbol) %>%
  summarise(
    n_mutations = n(),
    n_patients = n_distinct(patient_id),
    mean_VAF = round(mean(VAF, na.rm = TRUE), 3),
    median_VAF = round(median(VAF, na.rm = TRUE), 3),
    n_clonal = sum(VAF > 0.4, na.rm = TRUE),
    pct_clonal = round(100 * n_clonal / n(), 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_clonal))

cat("  ✓ Clonality calculated\n\n")

cat("Clonality Summary by Gene:\n")
print(vaf_summary, row.names = FALSE)
cat("\n")

# ─────────────────────────────────────────────────────────────────
# 5. Filter for Clonal Mutations (VAF > 0.4)
# ─────────────────────────────────────────────────────────────────

cat("[5] Filtering for clonal mutations (VAF > 0.4)...\n")

clonal_mutations <- maf_expressed[maf_expressed$VAF > 0.4, ]

cat("  ✓ Clonal mutations:", nrow(clonal_mutations), "\n")
cat("  ✓ Genes:", length(unique(clonal_mutations$Hugo_Symbol)), "\n")
cat("  ✓ Patients:", length(unique(clonal_mutations$patient_id)), "\n\n")

# Save all clonal mutations
write.csv(clonal_mutations,
          "data/processed/clonal_mutations.csv",
          row.names = FALSE)

cat("  ✓ Saved: data/processed/clonal_mutations.csv\n\n")

# ─────────────────────────────────────────────────────────────────
# 6. Select Final Candidate Genes
# ─────────────────────────────────────────────────────────────────

cat("[6] Selecting final candidate genes...\n")

# Focus on genes with:
# - High expression (TPM > 10)
# - High clonality (>40% mutations are clonal)
# - Multiple patients affected

candidate_genes <- vaf_summary %>%
  filter(pct_clonal > 40, n_patients >= 10) %>%
  arrange(desc(n_clonal))

cat("  ✓ Candidate genes (clonality > 40%, ≥10 patients):\n")
print(candidate_genes, row.names = FALSE)
cat("\n")

# Final candidates: mutations in candidate genes
final_candidates <- clonal_mutations[clonal_mutations$Hugo_Symbol %in% candidate_genes$Hugo_Symbol, ]

cat("  ✓ Final candidate mutations:", nrow(final_candidates), "\n")
cat("  ✓ Genes:", paste(unique(final_candidates$Hugo_Symbol), collapse = ", "), "\n")
cat("  ✓ Patients:", length(unique(final_candidates$patient_id)), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 7. Parse Amino Acid Changes
# ─────────────────────────────────────────────────────────────────

cat("[7] Parsing amino acid changes...\n")

# Function to extract AA info
extract_aa_info <- function(hgvsp) {
  if (is.na(hgvsp) || hgvsp == "") return(list(wt_aa = NA, mt_aa = NA, position = NA))
  
  clean <- gsub("^p\\.", "", hgvsp)
  
  # Missense: p.R130Q
  if (grepl("^[A-Z][0-9]+[A-Z*]$", clean)) {
    wt_aa <- substr(clean, 1, 1)
    mt_aa <- substr(clean, nchar(clean), nchar(clean))
    position <- as.numeric(gsub("[A-Z*]", "", clean))
  } 
  # Frameshift
  else if (grepl("fs", clean)) {
    wt_aa <- substr(clean, 1, 1)
    mt_aa <- "fs"
    position <- as.numeric(gsub("[^0-9]", "", clean))
  } 
  else {
    return(list(wt_aa = NA, mt_aa = NA, position = NA))
  }
  
  return(list(wt_aa = wt_aa, mt_aa = mt_aa, position = position))
}

# Apply to all
aa_info <- do.call(rbind, lapply(final_candidates$HGVSp_Short, extract_aa_info))
final_candidates$wt_aa <- unlist(aa_info[, 1])
final_candidates$mt_aa <- unlist(aa_info[, 2])
final_candidates$position <- as.numeric(aa_info[, 3])

cat("  ✓ Amino acid changes parsed\n\n")

# ─────────────────────────────────────────────────────────────────
# 8. Save Final Candidates
# ─────────────────────────────────────────────────────────────────

cat("[8] Saving evolution-resistant candidates...\n")

write.csv(final_candidates,
          "data/processed/Evolution_Resistant_Candidates.csv",
          row.names = FALSE)

write.csv(vaf_summary,
          "results/tables/Table2_Clonality_Summary.csv",
          row.names = FALSE)

cat("  ✓ Saved: data/processed/Evolution_Resistant_Candidates.csv\n")
cat("  ✓ Saved: results/tables/Table2_Clonality_Summary.csv\n\n")

# ─────────────────────────────────────────────────────────────────
# 9. Generate VAF Distribution Plot
# ─────────────────────────────────────────────────────────────────

cat("[9] Generating VAF distribution plot...\n")

# Only plot candidate genes
plot_data <- maf_expressed[maf_expressed$Hugo_Symbol %in% candidate_genes$Hugo_Symbol, ]

p <- ggplot(plot_data, aes(x = VAF, fill = Hugo_Symbol)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red", size = 1) +
  facet_wrap(~Hugo_Symbol, scales = "free_y") +
  theme_bw() +
  labs(
    title = "Variant Allele Frequency Distribution",
    subtitle = "Red line: VAF = 0.4 (clonality threshold)",
    x = "Variant Allele Frequency (VAF)",
    y = "Number of Mutations",
    fill = "Gene"
  ) +
  theme(legend.position = "none")

ggsave("results/figures/Figure2_VAF_Distribution.pdf", p, width = 10, height = 6)

cat("  ✓ Saved: results/figures/Figure2_VAF_Distribution.pdf\n\n")

# ─────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat("CLONALITY ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Files created:\n")
cat("  • data/processed/clonal_mutations.csv\n")
cat("  • data/processed/Evolution_Resistant_Candidates.csv\n")
cat("  • results/tables/Table2_Clonality_Summary.csv\n")
cat("  • results/figures/Figure2_VAF_Distribution.pdf\n\n")

cat("Evolution-resistant candidates:\n")
cat("  Total mutations:", nrow(final_candidates), "\n")
cat("  Genes:", paste(unique(final_candidates$Hugo_Symbol), collapse = ", "), "\n")
cat("  Patients:", length(unique(final_candidates$patient_id)), "\n\n")

cat("Next step: Run scripts/05_neoantigen_prediction.R\n\n")

cat("═══════════════════════════════════════════════════════════════\n\n")