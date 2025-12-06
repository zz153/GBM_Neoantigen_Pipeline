#!/usr/bin/env Rscript
################################################################################
# 01_download_data.R
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("STEP 1: DOWNLOADING TCGA-GBM DATA\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

# ─────────────────────────────────────────────────────────────────
# Set working directory (CHANGE THIS PATH FOR YOUR SYSTEM)
# ─────────────────────────────────────────────────────────────────

# Option 1: Manual path (edit this)
project_dir <- "/Users/ranzo85p/Library/CloudStorage/OneDrive-UniversityofOtago/Documents/GBM/GBM_newmodels/GBM_Neoantigen_Pipeline"

# Option 2: Use current directory
# project_dir <- getwd()

setwd(project_dir)

cat("Working directory:", getwd(), "\n\n")

# Create output directory
dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────
# 1. Download Mutation Data
# ─────────────────────────────────────────────────────────────────

cat("[1.1] Querying TCGA-GBM mutation data...\n")

query_mut <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

cat("  ✓ Found", length(getResults(query_mut)$cases), "cases\n\n")

cat("[1.2] Downloading mutation files...\n")
cat("  This may take 5-10 minutes...\n\n")

GDCdownload(query_mut)

cat("\n  ✓ Download complete\n\n")

cat("[1.3] Preparing mutation data...\n")

maf <- GDCprepare(query_mut)

cat("  ✓ Loaded", nrow(maf), "mutations from", 
    length(unique(maf$Tumor_Sample_Barcode)), "patients\n\n")

cat("[1.4] Saving processed mutation data...\n")

write.csv(maf, "data/raw/TCGA_mutation_maf.csv", row.names = FALSE)
saveRDS(maf, "data/raw/TCGA_mutation_maf.rds")

cat("  ✓ Saved: data/raw/TCGA_mutation_maf.csv\n")
cat("  ✓ Saved: data/raw/TCGA_mutation_maf.rds\n\n")

# ─────────────────────────────────────────────────────────────────
# 2. Download Expression Data
# ─────────────────────────────────────────────────────────────────

cat("[2.1] Querying TCGA-GBM expression data...\n")

query_tcga <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

cat("  ✓ Found", nrow(getResults(query_tcga)), "samples\n\n")

cat("[2.2] Downloading expression files...\n")
cat("  This may take 15-25 minutes...\n\n")

GDCdownload(query_tcga)

cat("\n  ✓ Download complete\n\n")

saveRDS(query_tcga, "data/raw/TCGA_GBM_query.rds")

cat("[2.3] Preparing expression data...\n")
cat("  This may take 5-10 minutes...\n\n")

tcga_data <- GDCprepare(query_tcga)

saveRDS(tcga_data, "data/raw/TCGA_GBM_data.rds")

cat("  ✓ Saved: data/raw/TCGA_GBM_data.rds\n\n")

# ─────────────────────────────────────────────────────────────────
# 3. Filter IDH-Wildtype
# ─────────────────────────────────────────────────────────────────

cat("[3] Filtering IDH-wildtype patients...\n")

sample_info <- as.data.frame(colData(tcga_data))
idhwt_samples <- sample_info[sample_info$paper_IDH.status == "WT" & 
                               !is.na(sample_info$paper_IDH.status), ]

cat("  ✓ IDH-WT patients:", length(unique(idhwt_samples$patient)), "\n")

idhwt_patient_ids <- unique(idhwt_samples$patient)
writeLines(idhwt_patient_ids, "data/raw/idhwt_patient_ids.txt")
saveRDS(idhwt_samples, "data/raw/idhwt_samples.rds")

cat("  ✓ Saved: data/raw/idhwt_patient_ids.txt\n")
cat("  ✓ Saved: data/raw/idhwt_samples.rds\n\n")

cat("═══════════════════════════════════════════════════════════════\n")
cat("DOWNLOAD COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Files saved in:\n")
cat("  GDCdata/           (raw downloads)\n")
cat("  data/raw/          (processed files)\n\n")

cat("Next step: Run scripts/02_filter_mutations.R\n\n")