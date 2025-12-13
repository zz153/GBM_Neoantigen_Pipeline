#!/usr/bin/env Rscript
################################################################################
# install_packages.R
# 
# Automated installation of all R package dependencies for GBM Neoantigen Pipeline
# 
# Usage: Rscript install_packages.R
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("GBM NEOANTIGEN PIPELINE - PACKAGE INSTALLATION\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("R Version:", R.version.string, "\n\n")

# ─────────────────────────────────────────────────────────────────
# Check R Version
# ─────────────────────────────────────────────────────────────────

if (R.Version()$major < 4 || (R.Version()$major == 4 && as.numeric(R.Version()$minor) < 3)) {
  cat("⚠ WARNING: R version 4.3.0 or higher is recommended\n")
  cat("  Current version:", R.version.string, "\n")
  cat("  Some packages may not install correctly\n\n")
}

# ─────────────────────────────────────────────────────────────────
# 1. Install BiocManager
# ─────────────────────────────────────────────────────────────────

cat("[1/4] Installing BiocManager...\n")

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
  cat("  ✓ BiocManager installed\n\n")
} else {
  cat("  ✓ BiocManager already installed\n\n")
}

# ─────────────────────────────────────────────────────────────────
# 2. Install Bioconductor Packages
# ─────────────────────────────────────────────────────────────────

cat("[2/4] Installing Bioconductor packages...\n")
cat("  This may take 10-15 minutes...\n\n")

bioc_packages <- c(
  "TCGAbiolinks",           # TCGA data download
  "SummarizedExperiment",   # Genomic data structures
  "recount3",               # GTEx data access
  "edgeR"                   # RNA-seq normalization
)

bioc_failed <- character()

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  • Installing", pkg, "...\n")
    tryCatch({
      BiocManager::install(pkg, update = FALSE, ask = FALSE, force = TRUE)
      cat("    ✓ Success\n")
    }, error = function(e) {
      cat("    ✗ Failed:", e$message, "\n")
      bioc_failed <<- c(bioc_failed, pkg)
    })
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

cat("\n")

# ─────────────────────────────────────────────────────────────────
# 3. Install CRAN Packages
# ─────────────────────────────────────────────────────────────────

cat("[3/4] Installing CRAN packages...\n")
cat("  This may take 5-10 minutes...\n\n")

cran_packages <- c(
  # Core data manipulation
  "dplyr",
  "tidyr",
  
  # Visualization
  "ggplot2",
  "cowplot",
  "gridExtra",
  "ggpubr",
  "RColorBrewer",
  "scales",
  "patchwork",
  "ggrepel",
  
  # Survival analysis
  "survival",
  "survminer",
  
  # Data retrieval
  "httr",
  "jsonlite"
)

cran_failed <- character()

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  • Installing", pkg, "...\n")
    tryCatch({
      install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
      cat("    ✓ Success\n")
    }, error = function(e) {
      cat("    ✗ Failed:", e$message, "\n")
      cran_failed <<- c(cran_failed, pkg)
    })
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

cat("\n")

# ─────────────────────────────────────────────────────────────────
# 4. Verify Installation
# ─────────────────────────────────────────────────────────────────

cat("[4/4] Verifying installation...\n\n")

all_packages <- c(bioc_packages, cran_packages)
verification_failed <- character()

for (pkg in all_packages) {
  suppressWarnings({
    loaded <- require(pkg, character.only = TRUE, quietly = TRUE)
  })
  
  if (!loaded) {
    verification_failed <- c(verification_failed, pkg)
    cat("  ✗", pkg, "- NOT FOUND\n")
  } else {
    cat("  ✓", pkg, "- OK\n")
  }
}

cat("\n")

# ─────────────────────────────────────────────────────────────────
# Summary Report
# ─────────────────────────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat("INSTALLATION SUMMARY\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

total_packages <- length(all_packages)
successful <- total_packages - length(verification_failed)

cat("Total packages:", total_packages, "\n")
cat("Successful:", successful, "\n")
cat("Failed:", length(verification_failed), "\n\n")

if (length(verification_failed) == 0) {
  cat("✓ ALL PACKAGES INSTALLED SUCCESSFULLY!\n\n")
  cat("Next steps:\n")
  cat("  1. Install NetMHCpan-4.2 (see docs/INSTALLATION.md)\n")
  cat("  2. Install NetMHCIIpan-4.3 (see docs/INSTALLATION.md)\n")
  cat("  3. Update paths in scripts/05_neoantigen_prediction.R\n")
  cat("  4. Run: Rscript scripts/01_download_data.R\n\n")
  
} else {
  cat("⚠ WARNING: Some packages failed to install\n\n")
  
  cat("Failed packages:\n")
  for (pkg in verification_failed) {
    cat("  •", pkg, "\n")
  }
  cat("\n")
  
  cat("Manual installation:\n")
  cat("  R\n")
  cat("  BiocManager::install(c('", paste(verification_failed, collapse = "', '"), "'))\n\n")
  
  cat("If problems persist:\n")
  cat("  1. Update R to >= 4.3.0\n")
  cat("  2. Update BiocManager: BiocManager::install(version = '3.18')\n")
  cat("  3. Check system dependencies (see docs/INSTALLATION.md)\n")
  cat("  4. Contact: zohaib.rana@otago.ac.nz\n\n")
}

# ─────────────────────────────────────────────────────────────────
# Print Session Info
# ─────────────────────────────────────────────────────────────────

cat("Session Information:\n")
cat("───────────────────────────────────────────────────────────────\n")
print(sessionInfo())
cat("\n")

cat("═══════════════════════════════════════════════════════════════\n\n")
