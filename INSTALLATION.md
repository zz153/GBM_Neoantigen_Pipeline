# Installation Guide

Complete installation instructions for the GBM Neoantigen Pipeline.

---

## Table of Contents

1. [System Requirements](#system-requirements)
2. [R Installation](#r-installation)
3. [R Package Dependencies](#r-package-dependencies)
4. [NetMHCpan Tools](#netmhcpan-tools)
5. [Verify Installation](#verify-installation)
6. [Troubleshooting](#troubleshooting)

---

## System Requirements

### Minimum Requirements

- **OS**: Linux, macOS, or Windows (WSL recommended)
- **RAM**: 16 GB minimum, 32 GB recommended
- **Storage**: 50 GB free space (TCGA data is ~20 GB)
- **Internet**: Stable connection for TCGA downloads

### Recommended Specifications

- **CPU**: 4+ cores
- **RAM**: 32 GB
- **Storage**: 100 GB SSD

---

## R Installation

### Install R (≥ 4.3.0)

#### macOS

```bash
# Using Homebrew
brew install r

# Or download from CRAN
# https://cran.r-project.org/bin/macosx/
```

#### Linux (Ubuntu/Debian)

```bash
# Add CRAN repository
sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
sudo add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/'

# Install R
sudo apt update
sudo apt install r-base r-base-dev

# Verify
R --version
```

#### Windows

1. Download R from https://cran.r-project.org/bin/windows/base/
2. Run installer
3. Add R to PATH

### Install RStudio (Optional but Recommended)

Download from: https://posit.co/download/rstudio-desktop/

---

## R Package Dependencies

### Method 1: Automatic Installation (Recommended)

Save this as `install_packages.R` in the repository root:

```r
#!/usr/bin/env Rscript
################################################################################
# install_packages.R
# Installs all required R packages for GBM Neoantigen Pipeline
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("INSTALLING R PACKAGES FOR GBM NEOANTIGEN PIPELINE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────
# 1. Install BiocManager
# ─────────────────────────────────────────────────────────────────

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("[1] Installing BiocManager...\n")
  install.packages("BiocManager")
} else {
  cat("[1] BiocManager already installed\n")
}

# ─────────────────────────────────────────────────────────────────
# 2. Install Bioconductor Packages
# ─────────────────────────────────────────────────────────────────

cat("\n[2] Installing Bioconductor packages...\n")

bioc_packages <- c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "recount3",
  "edgeR"
)

for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  • Installing", pkg, "...\n")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

# ─────────────────────────────────────────────────────────────────
# 3. Install CRAN Packages
# ─────────────────────────────────────────────────────────────────

cat("\n[3] Installing CRAN packages...\n")

cran_packages <- c(
  "survival",
  "survminer",
  "ggplot2",
  "dplyr",
  "tidyr",
  "cowplot",
  "gridExtra",
  "httr",
  "jsonlite",
  "ggpubr",
  "RColorBrewer",
  "scales",
  "patchwork",
  "ggrepel"
)

for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  • Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  } else {
    cat("  ✓", pkg, "already installed\n")
  }
}

# ─────────────────────────────────────────────────────────────────
# 4. Verify Installation
# ─────────────────────────────────────────────────────────────────

cat("\n[4] Verifying installation...\n")

all_packages <- c(bioc_packages, cran_packages)
failed <- character()

for (pkg in all_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    failed <- c(failed, pkg)
  }
}

if (length(failed) == 0) {
  cat("\n✓ ALL PACKAGES INSTALLED SUCCESSFULLY!\n\n")
} else {
  cat("\n⚠ WARNING: The following packages failed to install:\n")
  cat(paste("  -", failed, collapse = "\n"), "\n")
  cat("\nPlease install manually:\n")
  cat("  BiocManager::install(c('", paste(failed, collapse = "', '"), "'))\n\n")
}

cat("═══════════════════════════════════════════════════════════════\n\n")
```

Then run:

```bash
Rscript install_packages.R
```

### Method 2: Manual Installation

```r
# Start R
R

# Install BiocManager
install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
  "TCGAbiolinks",
  "SummarizedExperiment",
  "recount3",
  "edgeR"
))

# Install CRAN packages
install.packages(c(
  "survival",
  "survminer",
  "ggplot2",
  "dplyr",
  "tidyr",
  "cowplot",
  "gridExtra",
  "httr",
  "jsonlite",
  "ggpubr",
  "RColorBrewer",
  "scales",
  "patchwork",
  "ggrepel"
))

# Verify
library(TCGAbiolinks)
library(ggplot2)
library(survival)
```

---

## NetMHCpan Tools

### NetMHCpan-4.2 (HLA Class I Prediction)

**CRITICAL**: NetMHCpan tools require academic license registration.

#### Step 1: Request Access

1. Visit: https://services.healthtech.dtu.dk/services/NetMHCpan-4.2/
2. Click "Downloads" → "Request Academic License"
3. Fill in institutional email and details
4. Wait for approval email (typically 1-2 business days)

#### Step 2: Download and Install

```bash
# After receiving download link via email:

# Download
wget <download-link-from-email>

# Extract
tar -xzf netMHCpan-4.2.Linux.tar.gz

# Move to desired location
sudo mv netMHCpan-4.2 /opt/netMHCpan-4.2

# Test
/opt/netMHCpan-4.2/netMHCpan -h
```

#### Step 3: Update Pipeline Path

Edit `scripts/05_neoantigen_prediction.R`:

```r
# Line 218 - Update this path
netmhcpan_path <- "/opt/netMHCpan-4.2/netMHCpan"
```

### NetMHCIIpan-4.3 (HLA Class II Prediction)

#### Step 1: Request Access

1. Visit: https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/
2. Same registration process as NetMHCpan-4.2

#### Step 2: Download and Install

```bash
# After receiving download link:

# Download
wget <download-link-from-email>

# Extract
tar -xzf netMHCIIpan-4.3.Linux.tar.gz

# Move to desired location
sudo mv netMHCIIpan-4.3 /opt/netMHCIIpan-4.3

# Test
/opt/netMHCIIpan-4.3/netMHCIIpan -h
```

#### Step 3: Update Pipeline Path

Edit `scripts/05_neoantigen_prediction.R`:

```r
# Line 290 - Update this path
netmhciipan_path <- "/opt/netMHCIIpan-4.3/netMHCIIpan"
```

### Alternative: Docker Container (Coming Soon)

We're working on a Docker container with pre-installed NetMHCpan tools for easier setup.

---

## Verify Installation

### Test R Packages

```bash
Rscript -e "
  packages <- c('TCGAbiolinks', 'survival', 'ggplot2', 'dplyr')
  for(pkg in packages) {
    if(require(pkg, character.only=TRUE)) {
      cat(sprintf('✓ %s OK\n', pkg))
    } else {
      cat(sprintf('✗ %s FAILED\n', pkg))
    }
  }
"
```

### Test NetMHCpan

```bash
# Test NetMHCpan-4.2
/opt/netMHCpan-4.2/netMHCpan -p test/test_peptides.txt -a HLA-A02:01

# Test NetMHCIIpan-4.3
/opt/netMHCIIpan-4.3/netMHCIIpan -inptype 1 -f test/test_peptides.txt -a DRB1_0101
```

### Run Test Pipeline

```bash
# Quick test (uses small subset of data)
Rscript tests/quick_test.R
```

---

## Troubleshooting

### Issue 1: BiocManager Installation Fails

**Solution**:

```r
# Update R to latest version
# Then try:
options(repos = c(CRAN = "https://cloud.r-project.org"))
install.packages("BiocManager")
```

### Issue 2: TCGAbiolinks Download Timeout

**Solution**:

```r
# In scripts/01_download_data.R, add:
options(timeout = 600)  # 10 minutes
Sys.setenv(VROOM_CONNECTION_SIZE = 131072 * 2)
```

### Issue 3: NetMHCpan "Command not found"

**Solution**:

```bash
# Check installation
ls -la /opt/netMHCpan-4.2/netMHCpan

# Make executable
chmod +x /opt/netMHCpan-4.2/netMHCpan

# Test
/opt/netMHCpan-4.2/netMHCpan -h
```

### Issue 4: Memory Errors

**Solution**:

```r
# Increase R memory (Linux/Mac)
ulimit -s unlimited

# Increase R memory (Windows)
memory.limit(size = 16000)  # 16 GB
```

### Issue 5: SSL Certificate Errors (TCGA Download)

**Solution**:

```r
# In scripts/01_download_data.R, add:
httr::set_config(httr::config(ssl_verifypeer = FALSE))
```

### Issue 6: Missing System Dependencies (Linux)

**Solution**:

```bash
# Ubuntu/Debian
sudo apt-get install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libfontconfig1-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev
```

---

## Getting Help

If you encounter issues not covered here:

1. Check [FAQ.md](FAQ.md)
2. Search [GitHub Issues](https://github.com/zz153/GBM_Neoantigen_Pipeline/issues)
3. Open a new issue with:
   - R version (`R --version`)
   - OS (`uname -a` or `systeminfo`)
   - Error message
   - Script that failed

---

**Last Updated**: December 2024
