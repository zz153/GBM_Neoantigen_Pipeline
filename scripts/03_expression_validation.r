#!/usr/bin/env Rscript
################################################################################
# 03_expression_validation.R
# 
# Validates tumor-specific expression by comparing:
#   - TCGA Tumor (IDH-WT)
#   - TCGA Normal
#   - GTEx Normal Brain
# 
# Inputs:
#   - data/raw/TCGA_GBM_data.rds
#   - data/raw/idhwt_samples.rds
#   - data/processed/Top50_Mutated_Genes.csv
# 
# Outputs:
#   - data/processed/expression_summary.csv
#   - results/figures/Figure1_Expression_Validation.pdf
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("STEP 3: EXPRESSION VALIDATION (TUMOR VS NORMAL)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────────────────────────

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(recount3)
  library(edgeR)
})

# Create output directories
dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────
# 1. Load TCGA Data
# ─────────────────────────────────────────────────────────────────

cat("[1] Loading TCGA expression data...\n")

if (!file.exists("data/raw/TCGA_GBM_data.rds")) {
  stop("ERROR: data/raw/TCGA_GBM_data.rds not found.\n",
       "       Please run scripts/01_download_data.R first.")
}

tcga_data <- readRDS("data/raw/TCGA_GBM_data.rds")

# Get ALL sample info (not just IDH-WT)
sample_info <- as.data.frame(colData(tcga_data))

cat("  ✓ Loaded expression data\n")
cat("  ✓ Genes:", nrow(tcga_data), "\n")
cat("  ✓ Total samples:", ncol(tcga_data), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 2. Separate Tumor vs Normal Samples
# ─────────────────────────────────────────────────────────────────

cat("[2] Separating sample types...\n")

# Get IDH-WT samples for TUMOR only
idhwt_samples <- readRDS("data/raw/idhwt_samples.rds")
tumor_samples <- idhwt_samples[idhwt_samples$sample_type == "Primary Tumor", ]

# Get ALL TCGA normal samples (regardless of IDH status)
# Normal brain tissue shouldn't have IDH mutations anyway
normal_samples <- sample_info[sample_info$sample_type == "Solid Tissue Normal", ]

cat("  ✓ TCGA Normal samples:", nrow(normal_samples), "\n")
cat("  ✓ TCGA Tumor samples (IDH-WT):", nrow(tumor_samples), "\n\n")

if (nrow(normal_samples) == 0) {
  cat("  ⚠ WARNING: No TCGA normal samples found!\n")
  cat("  ⚠ Will compare tumor to GTEx only\n\n")
}

# ─────────────────────────────────────────────────────────────────
# 3. Extract Expression for Target Genes
# ─────────────────────────────────────────────────────────────────

cat("[3] Extracting expression for target genes...\n")

# Load top mutated genes
top_genes <- read.csv("data/processed/Top50_Mutated_Genes.csv")

# Focus on top 20 genes
top20_genes <- head(top_genes$Gene, 20)
cat("  ✓ Analyzing top 20 genes\n")

# Get TPM values
tpm_all <- assay(tcga_data, "tpm_unstrand")
gene_names <- rowData(tcga_data)$gene_name

# Check expression levels IN TUMOR
expr_levels <- data.frame(
  Gene = top20_genes,
  Mean_TPM = sapply(top20_genes, function(g) {
    gene_idx <- which(gene_names == g)
    if (length(gene_idx) == 0) return(NA)
    mean(tpm_all[gene_idx, rownames(tumor_samples)], na.rm = TRUE)
  }),
  stringsAsFactors = FALSE
)

# Filter for expressed genes (TPM > 10 in tumors)
expressed_genes <- expr_levels$Gene[expr_levels$Mean_TPM > 10 & !is.na(expr_levels$Mean_TPM)]

cat("  ✓ Genes with TPM > 10 in tumors:", length(expressed_genes), "\n")

if (length(expressed_genes) > 0) {
  cat("  ✓ Expressed genes:", paste(expressed_genes, collapse = ", "), "\n\n")
} else {
  cat("  ⚠ No genes with TPM > 10, lowering threshold to TPM > 5\n")
  expressed_genes <- expr_levels$Gene[expr_levels$Mean_TPM > 5 & !is.na(expr_levels$Mean_TPM)]
}

# Default to top 4 if still none found
if (length(expressed_genes) == 0) {
  cat("  ⚠ Using default genes: TP53, PTEN, EGFR, RB1\n\n")
  expressed_genes <- c("TP53", "PTEN", "EGFR", "RB1")
}

# Limit to top 4 expressed genes for cleaner visualization
target_genes <- head(expressed_genes, 4)
cat("  ✓ Selected genes for validation:\n")
for(i in 1:length(target_genes)) {
  mean_tpm <- expr_levels$Mean_TPM[expr_levels$Gene == target_genes[i]]
  cat("     ", i, ".", target_genes[i], " (", round(mean_tpm, 1), " TPM)\n", sep = "")
}
cat("\n")

# ─────────────────────────────────────────────────────────────────
# 4. Build Expression Dataframe (TCGA)
# ─────────────────────────────────────────────────────────────────

cat("[4] Building expression comparison dataframe...\n")

expr_df <- data.frame()

for(gene in target_genes) {
  gene_idx <- which(gene_names == gene)
  
  if (length(gene_idx) == 0) {
    cat("  ⚠ Gene", gene, "not found in expression data\n")
    next
  }
  
  # TCGA Normal
  if (nrow(normal_samples) > 0) {
    normal_expr <- tpm_all[gene_idx, rownames(normal_samples)]
    expr_df <- rbind(expr_df,
                     data.frame(Gene = gene, 
                                Expression = as.numeric(normal_expr), 
                                Type = "TCGA_Normal",
                                stringsAsFactors = FALSE))
    cat("  ✓", gene, "- TCGA Normal:", length(normal_expr), "samples\n")
  }
  
  # TCGA Tumor (IDH-WT)
  tumor_expr <- tpm_all[gene_idx, rownames(tumor_samples)]
  expr_df <- rbind(expr_df,
                   data.frame(Gene = gene, 
                              Expression = as.numeric(tumor_expr), 
                              Type = "TCGA_Tumor",
                              stringsAsFactors = FALSE))
  cat("  ✓", gene, "- TCGA Tumor:", length(tumor_expr), "samples\n")
}

cat("\n  ✓ TCGA expression extracted:", nrow(expr_df), "total datapoints\n\n")

# ─────────────────────────────────────────────────────────────────
# 5. Add GTEx Brain Expression
# ─────────────────────────────────────────────────────────────────

cat("[5] Downloading GTEx brain expression data...\n")
cat("  This may take 5-10 minutes...\n\n")

gtex_success <- FALSE

tryCatch({
  # Get GTEx brain data
  gtex_projects <- available_projects()
  gtex_brain <- gtex_projects[gtex_projects$file_source == "gtex" & 
                                gtex_projects$project == "BRAIN", ]
  
  if (nrow(gtex_brain) == 0) {
    stop("GTEx BRAIN project not found")
  }
  
  cat("  ✓ Found GTEx brain project\n")
  cat("  ✓ Downloading data...\n")
  
  gtex_rse <- create_rse(gtex_brain)
  
  cat("  ✓ GTEx data downloaded\n")
  cat("  ✓ Samples:", ncol(gtex_rse), "\n")
  
  # Calculate TPM from raw counts
  cat("  ✓ Calculating TPM...\n")
  
  counts <- assay(gtex_rse, "raw_counts")
  gene_lengths <- rowData(gtex_rse)$bp_length
  
  # Calculate RPK (Reads Per Kilobase)
  rpk <- counts / (gene_lengths / 1000)
  
  # Calculate TPM
  gtex_tpm <- t(t(rpk) / colSums(rpk)) * 1e6
  
  gtex_genes <- rowData(gtex_rse)$gene_name
  
  cat("  ✓ TPM calculated\n")
  
  # Extract for target genes
  cat("  ✓ Extracting GTEx expression...\n")
  
  for(gene in target_genes) {
    gene_idx <- which(gtex_genes == gene)
    
    if (length(gene_idx) == 0) {
      cat("  ⚠ Gene", gene, "not found in GTEx\n")
      next
    }
    
    gtex_expr <- gtex_tpm[gene_idx, ]
    
    expr_df <- rbind(expr_df,
                     data.frame(Gene = gene, 
                                Expression = as.numeric(gtex_expr), 
                                Type = "GTEx_Brain",
                                stringsAsFactors = FALSE))
    cat("  ✓", gene, "- GTEx Brain:", length(gtex_expr), "samples\n")
  }
  
  cat("  ✓ GTEx expression extracted\n\n")
  gtex_success <- TRUE
  
}, error = function(e) {
  cat("  ⚠ WARNING: Could not download GTEx data\n")
  cat("  ⚠ Error:", e$message, "\n")
  cat("  ⚠ Continuing with TCGA data only\n\n")
})

# ─────────────────────────────────────────────────────────────────
# 6. Calculate Summary Statistics
# ─────────────────────────────────────────────────────────────────

cat("[6] Calculating summary statistics...\n")

expr_summary <- expr_df %>%
  group_by(Gene, Type) %>%
  summarise(
    Mean_TPM = round(mean(Expression, na.rm = TRUE), 2),
    Median_TPM = round(median(Expression, na.rm = TRUE), 2),
    SD_TPM = round(sd(Expression, na.rm = TRUE), 2),
    N_samples = n(),
    .groups = "drop"
  ) %>%
  arrange(Gene, factor(Type, levels = c("GTEx_Brain", "TCGA_Normal", "TCGA_Tumor")))

write.csv(expr_summary, 
          "data/processed/expression_summary.csv", 
          row.names = FALSE)

cat("  ✓ Saved: data/processed/expression_summary.csv\n\n")

cat("Expression Summary Table:\n")
print(expr_summary, row.names = FALSE)

cat("\n")

# ─────────────────────────────────────────────────────────────────
# 7. Generate Expression Comparison Plot
# ─────────────────────────────────────────────────────────────────

cat("[7] Generating expression comparison plot...\n")

# Set factor order for plotting
type_levels <- unique(expr_df$Type)
if (gtex_success) {
  type_levels <- c("GTEx_Brain", "TCGA_Normal", "TCGA_Tumor")
} else {
  type_levels <- c("TCGA_Normal", "TCGA_Tumor")
}

expr_df$Type <- factor(expr_df$Type, levels = type_levels)

# Create plot
p <- ggplot(expr_df, aes(x = Type, y = Expression + 1, fill = Type)) +
  geom_boxplot(outlier.alpha = 0.3, notch = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.1, size = 0.5) +
  facet_wrap(~Gene, scales = "free_y", ncol = 2) +
  scale_y_log10(
    breaks = c(1, 10, 100, 1000),
    labels = c("0", "10", "100", "1000")
  ) +
  scale_fill_manual(
    values = c("GTEx_Brain" = "#2E9FDF", 
               "TCGA_Normal" = "#E7B800",
               "TCGA_Tumor" = "#FC4E07"),
    labels = c("GTEx Brain" = "GTEx Brain", 
               "TCGA_Normal" = "TCGA Normal",
               "TCGA_Tumor" = "TCGA Tumor")
  ) +
  theme_bw(base_size = 12) +
  labs(
    title = "Gene Expression: Tumor vs Normal Brain",
    subtitle = paste("Comparing", nrow(tumor_samples), "tumors,", 
                     nrow(normal_samples), "TCGA normals,",
                     if(gtex_success) paste(ncol(gtex_rse), "GTEx normals") else ""),
    y = "Expression (TPM + 1, log scale)",
    x = "",
    fill = "Sample Type"
  ) +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 11),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9, color = "grey40")
  )

ggsave("results/figures/Figure1_Expression_Validation.pdf", 
       p, 
       width = 10, 
       height = 8)

cat("  ✓ Saved: results/figures/Figure1_Expression_Validation.pdf\n\n")

# ─────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat("EXPRESSION VALIDATION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Files created:\n")
cat("  • data/processed/expression_summary.csv\n")
cat("  • results/figures/Figure1_Expression_Validation.pdf\n\n")

cat("Key findings:\n")
for(gene in target_genes) {
  tumor_mean <- expr_summary$Mean_TPM[expr_summary$Gene == gene & 
                                        expr_summary$Type == "TCGA_Tumor"]
  
  # Compare to TCGA Normal
  if ("TCGA_Normal" %in% expr_summary$Type) {
    tcga_normal_mean <- expr_summary$Mean_TPM[expr_summary$Gene == gene & 
                                                expr_summary$Type == "TCGA_Normal"]
    tcga_fc <- round(tumor_mean / tcga_normal_mean, 1)
    cat("  •", gene, ":\n")
    cat("     - Tumor:", tumor_mean, "TPM\n")
    cat("     - TCGA Normal:", tcga_normal_mean, "TPM (", tcga_fc, "x fold-change)\n", sep = "")
    
    # Compare to GTEx if available
    if ("GTEx_Brain" %in% expr_summary$Type) {
      gtex_mean <- expr_summary$Mean_TPM[expr_summary$Gene == gene & 
                                           expr_summary$Type == "GTEx_Brain"]
      gtex_fc <- round(tumor_mean / gtex_mean, 1)
      cat("     - GTEx Brain:", gtex_mean, "TPM (", gtex_fc, "x fold-change)\n", sep = "")
    }
  } else {
    cat("  •", gene, ":", tumor_mean, "TPM (tumor)\n", sep = "")
  }
}

cat("\n")
cat("Next step: Run scripts/04_clonality_analysis.R\n\n")

cat("═══════════════════════════════════════════════════════════════\n\n")




