#!/usr/bin/env Rscript
################################################################################
# 06_survival_analysis.R
# 
# Performs survival analysis by gene mutation status
# 
# Inputs:
#   - data/raw/idhwt_samples.rds
#   - data/processed/TCGA_mutation_maf_IDHwt.csv
# 
# Outputs:
#   - results/figures/Figure3_Survival_Analysis.pdf
#   - results/tables/Table3_Survival_Statistics.csv
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("STEP 6: SURVIVAL ANALYSIS\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────────────────────────

library(survival)
library(survminer)
library(dplyr)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────
# 1. Load Data
# ─────────────────────────────────────────────────────────────────

cat("[1] Loading clinical and mutation data...\n")

# Load clinical data
idhwt_samples <- readRDS("data/raw/idhwt_samples.rds")
maf_idhwt <- read.csv("data/processed/TCGA_mutation_maf_IDHwt.csv")

cat("  ✓ Clinical data:", nrow(idhwt_samples), "samples\n")
cat("  ✓ Mutations:", nrow(maf_idhwt), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 2. Prepare Clinical Data
# ─────────────────────────────────────────────────────────────────

cat("[2] Preparing survival data...\n")

# Filter for primary tumors only
clin_df <- idhwt_samples[idhwt_samples$sample_type == "Primary Tumor", ]

cat("  ✓ Primary tumor samples:", nrow(clin_df), "\n")

# Create survival variables
clin_df$deceased <- clin_df$vital_status == "Dead"
clin_df$OS_months <- ifelse(clin_df$deceased,
                            clin_df$days_to_death / 30.44,
                            clin_df$days_to_last_follow_up / 30.44)

# Remove patients with missing survival data
clin_df <- clin_df[!is.na(clin_df$OS_months) & clin_df$OS_months > 0, ]

cat("  ✓ Patients with survival data:", nrow(clin_df), "\n")
cat("  ✓ Deceased:", sum(clin_df$deceased), "\n")
cat("  ✓ Alive:", sum(!clin_df$deceased), "\n")
cat("  ✓ Median follow-up:", round(median(clin_df$OS_months), 1), "months\n\n")

# ─────────────────────────────────────────────────────────────────
# 3. Add Mutation Status
# ─────────────────────────────────────────────────────────────────

cat("[3] Adding mutation status for candidate genes...\n")

candidate_genes <- c("TP53", "PTEN", "EGFR")

for(gene in candidate_genes) {
  # Get patients with mutations in this gene
  mutated_patients <- unique(maf_idhwt$patient_id[maf_idhwt$Hugo_Symbol == gene])
  
  # Add mutation status to clinical data
  clin_df[[paste0(gene, "_mutated")]] <- clin_df$patient %in% mutated_patients
  
  n_mut <- sum(clin_df[[paste0(gene, "_mutated")]])
  cat("  ✓", gene, ":", n_mut, "mutated patients\n")
}

cat("\n")

# ─────────────────────────────────────────────────────────────────
# 4. Perform Survival Analysis
# ─────────────────────────────────────────────────────────────────

cat("[4] Performing Kaplan-Meier analysis...\n")

survival_results <- data.frame()

for(gene in candidate_genes) {
  cat("  • Analyzing", gene, "...\n")
  
  # Create survival object
  surv_obj <- Surv(time = clin_df$OS_months, event = clin_df$deceased)
  
  # Fit survival curve
  fit <- survfit(surv_obj ~ get(paste0(gene, "_mutated")), data = clin_df)
  
  # Log-rank test
  logrank <- survdiff(surv_obj ~ get(paste0(gene, "_mutated")), data = clin_df)
  
  # Calculate p-value
  pval <- 1 - pchisq(logrank$chisq, length(logrank$n) - 1)
  
  # Get median survival
  median_surv <- summary(fit)$table[, "median"]
  
  # Store results
  survival_results <- rbind(survival_results, data.frame(
    Gene = gene,
    N_mutated = sum(clin_df[[paste0(gene, "_mutated")]]),
    N_wildtype = sum(!clin_df[[paste0(gene, "_mutated")]]),
    Median_OS_mutated = round(median_surv[2], 1),
    Median_OS_wildtype = round(median_surv[1], 1),
    P_value = signif(pval, 3),
    Significant = ifelse(pval < 0.05, "Yes", "No"),
    stringsAsFactors = FALSE
  ))
}

cat("\n")

# ─────────────────────────────────────────────────────────────────
# 5. Save Results
# ─────────────────────────────────────────────────────────────────

cat("[5] Saving survival statistics...\n")

write.csv(survival_results, 
          "results/tables/Table3_Survival_Statistics.csv",
          row.names = FALSE)

cat("  ✓ Saved: results/tables/Table3_Survival_Statistics.csv\n\n")

cat("Survival Analysis Results:\n")
print(survival_results, row.names = FALSE)
cat("\n")

# ─────────────────────────────────────────────────────────────────
# 6. Generate Kaplan-Meier Plots
# ─────────────────────────────────────────────────────────────────

cat("[6] Generating Kaplan-Meier plots...\n")

# Create 2x2 plot layout
pdf("results/figures/SupplementaryFigure2_Survival_Analysis.pdf", width = 12, height = 10)

plot_list <- list()

for(i in 1:length(candidate_genes)) {
  gene <- candidate_genes[i]
  
  surv_obj <- Surv(time = clin_df$OS_months, event = clin_df$deceased)
  fit <- survfit(surv_obj ~ get(paste0(gene, "_mutated")), data = clin_df)
  
  pval <- survival_results$P_value[survival_results$Gene == gene]
  
  p <- ggsurvplot(
    fit,
    data = clin_df,
    pval = TRUE,
    pval.coord = c(0, 0.1),
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    legend.labs = c("Wildtype", "Mutated"),
    legend.title = gene,
    palette = c("#2E9FDF", "#E7B800"),
    title = paste(gene, "Mutation Status"),
    xlab = "Time (months)",
    ylab = "Overall Survival Probability",
    ggtheme = theme_bw()
  )
  
  plot_list[[i]] <- p
  
  cat("  ✓", gene, "plot created (p =", pval, ")\n")
}

# Arrange plots
arranged <- arrange_ggsurvplots(plot_list, ncol = 2, nrow = 2)
print(arranged)

dev.off()

cat("  ✓ Saved: results/figures/Figure3_Survival_Analysis.pdf\n\n")

# ─────────────────────────────────────────────────────────────────
# 7. Cox Proportional Hazards Model
# ─────────────────────────────────────────────────────────────────

cat("[7] Multivariate Cox regression...\n")

cox_model <- coxph(Surv(OS_months, deceased) ~ 
                     TP53_mutated + PTEN_mutated + EGFR_mutated +
                     age_at_diagnosis + gender,
                   data = clin_df)

cat("\nCox Proportional Hazards Model:\n")
print(summary(cox_model))
cat("\n")

# ─────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat("SURVIVAL ANALYSIS COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Files created:\n")
cat("  • results/figures/Figure3_Survival_Analysis.pdf\n")
cat("  • results/tables/Table3_Survival_Statistics.csv\n\n")

cat("Key findings:\n")
for(i in 1:nrow(survival_results)) {
  row <- survival_results[i, ]
  cat("  •", row$Gene, ":\n")
  cat("     - Mutated median OS:", row$Median_OS_mutated, "months\n")
  cat("     - Wildtype median OS:", row$Median_OS_wildtype, "months\n")
  cat("     - P-value:", row$P_value, 
      ifelse(row$Significant == "Yes", "(significant)", "(not significant)"), "\n")
}

cat("\n")
cat("Next step: Run scripts/07_sharing_analysis.R\n\n")

cat("═══════════════════════════════════════════════════════════════\n\n")

