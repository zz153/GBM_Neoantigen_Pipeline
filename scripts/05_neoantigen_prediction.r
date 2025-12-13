#!/usr/bin/env Rscript
################################################################################
# 05_neoantigen_prediction.R (FINAL FIX - DRB1 only)
# 
# Predicts HLA-binding neoantigens using NetMHCpan-4.2 and NetMHCIIpan-4.3
# 
# Steps:
#   1. Fetch protein sequences from UniProt
#   2. Generate 21-mer peptides (mutation centered)
#   3. Extract 9-mer peptides (Class I)
#   4. Extract 15-mer peptides (Class II)
#   5. Run NetMHCpan predictions (Class I)
#   6. Run NetMHCIIpan predictions (Class II) - DRB1 alleles only
#   7. Parse and merge results
#   8. Filter strong binders
# 
# Inputs:
#   - data/processed/Evolution_Resistant_Candidates.csv
# 
# Outputs:
#   - data/processed/protein_sequences.csv
#   - results/netmhcpan_output/peptides_9mer.txt
#   - results/netmhcpan_output/peptides_15mer.txt
#   - results/netmhcpan_output/predictions_class1.txt
#   - results/netmhcpan_output/predictions_class2.txt
#   - data/processed/NetMHCpan_All_Predictions.csv
#   - data/processed/Strong_Binders_Final.csv
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("STEP 5: NEOANTIGEN PREDICTION (NetMHCpan + NetMHCIIpan)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────────────────────────

library(dplyr)
library(httr)
library(jsonlite)
library(tidyr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/netmhcpan_output", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────
# 1. Load Clonal Mutations
# ─────────────────────────────────────────────────────────────────

cat("[1] Loading clonal mutations...\n")

final_candidates <- read.csv("data/processed/Evolution_Resistant_Candidates.csv")

cat("  ✓ Loaded", nrow(final_candidates), "clonal mutations\n")
cat("  ✓ Genes:", paste(unique(final_candidates$Hugo_Symbol), collapse = ", "), "\n")
cat("  ✓ Patients:", length(unique(final_candidates$Tumor_Sample_Barcode)), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 2. Fetch Protein Sequences from UniProt
# ─────────────────────────────────────────────────────────────────

cat("[2] Fetching protein sequences from UniProt...\n")

# Function to fetch sequences
fetch_uniprot <- function(gene_name) {
  Sys.sleep(0.3)  # Rate limiting
  
  query <- sprintf('gene:%s AND organism_id:9606 AND reviewed:true', gene_name)
  url <- sprintf(
    'https://rest.uniprot.org/uniprotkb/search?query=%s&fields=accession,gene_primary,sequence&format=json&size=1',
    URLencode(query)
  )
  
  response <- GET(url, timeout(10))
  if (status_code(response) != 200) return(NULL)
  
  data <- fromJSON(content(response, as = "text", encoding = "UTF-8"), 
                   simplifyVector = FALSE)
  
  if (is.null(data$results) || length(data$results) == 0) return(NULL)
  
  entry <- data$results[[1]]
  seq_value <- entry$sequence$value
  
  if (is.null(seq_value)) return(NULL)
  
  return(data.frame(
    gene = gene_name,
    sequence = seq_value,
    stringsAsFactors = FALSE
  ))
}

# Get unique genes
genes <- unique(final_candidates$Hugo_Symbol)
cat("  ✓ Fetching", length(genes), "genes from UniProt...\n")

sequences <- list()
for(gene in genes) {
  seq <- fetch_uniprot(gene)
  if(!is.null(seq)) {
    sequences[[gene]] <- seq
    cat("    - ", gene, ": ", nchar(seq$sequence), " aa\n", sep = "")
  }
}

sequences_df <- do.call(rbind, sequences)

write.csv(sequences_df, "data/processed/protein_sequences.csv", row.names = FALSE)

cat("  ✓ Retrieved", nrow(sequences_df), "sequences\n")
cat("  ✓ Saved: data/processed/protein_sequences.csv\n\n")

# ─────────────────────────────────────────────────────────────────
# 3. Generate 21-mer Peptides (Mutation Centered)
# ─────────────────────────────────────────────────────────────────

cat("[3] Generating 21-mer peptides...\n")

# Function to extract peptide
extract_peptide <- function(sequence, position, mt_aa, peptide_length = 21) {
  if (is.na(sequence) || is.na(position) || is.na(mt_aa)) {
    return(NA)
  }
  
  flank <- floor((peptide_length - 1) / 2)
  start <- max(1, position - flank)
  end <- min(nchar(sequence), position + flank)
  
  peptide <- substr(sequence, start, end)
  
  # Replace WT with mutant AA
  center_pos <- position - start + 1
  if (center_pos > 0 && center_pos <= nchar(peptide)) {
    substr(peptide, center_pos, center_pos) <- mt_aa
  }
  
  # Pad with X if at edges
  if (nchar(peptide) < peptide_length) {
    if (start == 1) {
      padding <- paste(rep("X", peptide_length - nchar(peptide)), collapse = "")
      peptide <- paste0(padding, peptide)
    } else {
      padding <- paste(rep("X", peptide_length - nchar(peptide)), collapse = "")
      peptide <- paste0(peptide, padding)
    }
  }
  
  return(substr(peptide, 1, peptide_length))
}

# Add sequences to mutations
final_candidates <- final_candidates %>%
  left_join(sequences_df, by = c("Hugo_Symbol" = "gene"))

# Generate 21-mers
final_candidates$peptide_21mer <- mapply(
  extract_peptide,
  final_candidates$sequence,
  final_candidates$position,
  final_candidates$mt_aa
)

cat("  ✓ Generated", sum(!is.na(final_candidates$peptide_21mer)), "21-mer peptides\n\n")

# ─────────────────────────────────────────────────────────────────
# 4. Generate 9-mer Peptides (Class I)
# ─────────────────────────────────────────────────────────────────

cat("[4] Generating 9-mer peptides for HLA Class I...\n")

# Function to generate 9-mers containing mutation
generate_9mers <- function(peptide_21mer, mt_aa) {
  if (is.na(peptide_21mer) || is.na(mt_aa)) return(list())
  
  mut_pos <- 11  # Center of 21-mer
  peptides_9mer <- c()
  
  for (i in 1:13) {  # 13 possible 9-mers from 21-mer
    start <- i
    end <- i + 8
    peptide <- substr(peptide_21mer, start, end)
    
    # Only keep if mutation is within this 9-mer
    if (start <= mut_pos && mut_pos <= end) {
      peptides_9mer <- c(peptides_9mer, peptide)
    }
  }
  
  return(peptides_9mer)
}

# Generate 9-mers for each mutation
all_9mers <- data.frame()

for (i in 1:nrow(final_candidates)) {
  mut <- final_candidates[i, ]
  
  peptides <- generate_9mers(mut$peptide_21mer, mut$mt_aa)
  
  if (length(peptides) > 0) {
    for (pep in peptides) {
      all_9mers <- rbind(all_9mers, data.frame(
        gene = mut$Hugo_Symbol,
        patient = mut$Tumor_Sample_Barcode,
        mutation = mut$HGVSp_Short,
        peptide = pep,
        VAF = mut$VAF,
        peptide_length = 9,
        hla_class = "I",
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("  ✓ Generated", nrow(all_9mers), "9-mer peptides\n")

# Remove duplicates and stop codons
unique_peptides_9mer <- all_9mers %>%
  filter(!grepl("\\*", peptide)) %>%
  filter(!grepl("X", peptide)) %>%  # Remove peptides with X padding
  distinct(peptide, .keep_all = TRUE)

cat("  ✓ Unique clean 9-mers:", nrow(unique_peptides_9mer), "\n")

# Save
writeLines(unique_peptides_9mer$peptide, "results/netmhcpan_output/peptides_9mer.txt")
write.csv(all_9mers, "results/netmhcpan_output/peptides_9mer_metadata.csv", row.names = FALSE)

cat("  ✓ Saved: results/netmhcpan_output/peptides_9mer.txt\n\n")

# ─────────────────────────────────────────────────────────────────
# 5. Generate 15-mer Peptides (Class II)
# ─────────────────────────────────────────────────────────────────

cat("[5] Generating 15-mer peptides for HLA Class II...\n")

# Function to generate 15-mers containing mutation
generate_15mers <- function(peptide_21mer, mt_aa) {
  if (is.na(peptide_21mer) || is.na(mt_aa)) return(list())
  
  mut_pos <- 11  # Center of 21-mer
  peptides_15mer <- c()
  
  for (i in 1:7) {  # 7 possible 15-mers from 21-mer
    start <- i
    end <- i + 14
    peptide <- substr(peptide_21mer, start, end)
    
    # Only keep if mutation is within this 15-mer
    if (start <= mut_pos && mut_pos <= end) {
      peptides_15mer <- c(peptides_15mer, peptide)
    }
  }
  
  return(peptides_15mer)
}

# Generate 15-mers for each mutation
all_15mers <- data.frame()

for (i in 1:nrow(final_candidates)) {
  mut <- final_candidates[i, ]
  
  peptides <- generate_15mers(mut$peptide_21mer, mut$mt_aa)
  
  if (length(peptides) > 0) {
    for (pep in peptides) {
      all_15mers <- rbind(all_15mers, data.frame(
        gene = mut$Hugo_Symbol,
        patient = mut$Tumor_Sample_Barcode,
        mutation = mut$HGVSp_Short,
        peptide = pep,
        VAF = mut$VAF,
        peptide_length = 15,
        hla_class = "II",
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("  ✓ Generated", nrow(all_15mers), "15-mer peptides\n")

# Remove duplicates and stop codons
unique_peptides_15mer <- all_15mers %>%
  filter(!grepl("\\*", peptide)) %>%
  filter(!grepl("X", peptide)) %>%
  distinct(peptide, .keep_all = TRUE)

cat("  ✓ Unique clean 15-mers:", nrow(unique_peptides_15mer), "\n")

# Save
writeLines(unique_peptides_15mer$peptide, "results/netmhcpan_output/peptides_15mer.txt")
write.csv(all_15mers, "results/netmhcpan_output/peptides_15mer_metadata.csv", row.names = FALSE)

cat("  ✓ Saved: results/netmhcpan_output/peptides_15mer.txt\n\n")

# ─────────────────────────────────────────────────────────────────
# 6. Define HLA Alleles (Class I + Class II)
# ─────────────────────────────────────────────────────────────────

cat("[6] Defining HLA alleles...\n")

# HLA Class I alleles (common panel, ~70% coverage)
hla_class1 <- c(
  "HLA-A*02:01", "HLA-A*01:01", "HLA-A*24:02", "HLA-A*03:01",
  "HLA-B*07:02", "HLA-B*08:01", "HLA-B*44:02",
  "HLA-C*07:02", "HLA-C*07:01"
)

# HLA Class II alleles - DRB1 ONLY (DQB1 not available in NetMHCIIpan-4.3)
# Format for NetMHCIIpan: DRB1_0101 (underscore, no HLA- prefix)
hla_class2 <- c(
  "DRB1_0101", "DRB1_0301", "DRB1_0401",
  "DRB1_0701", "DRB1_1101", "DRB1_1301",
  "DRB1_1501", "DRB1_0901"
)

cat("  ✓ HLA Class I alleles:", length(hla_class1), "\n")
cat("  ✓ HLA Class II alleles (DRB1 only):", length(hla_class2), "\n")
cat("  ✓ Total HLA alleles:", length(hla_class1) + length(hla_class2), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 7. Run NetMHCpan (Class I)
# ─────────────────────────────────────────────────────────────────

cat("[7] Running NetMHCpan-4.2 predictions (Class I)...\n")

hla_class1_string <- paste(gsub("\\*", "", hla_class1), collapse = ",")

cat("  ✓ Peptides:", nrow(unique_peptides_9mer), "\n")
cat("  ✓ HLA alleles:", length(hla_class1), "\n")
cat("  ✓ Expected predictions:", nrow(unique_peptides_9mer) * length(hla_class1), "\n\n")

# Check NetMHCpan installation
netmhcpan_path <- "/Users/ranzo85p/netMHCpan-4.2/netMHCpan"

if (!file.exists(netmhcpan_path)) {
  stop("ERROR: NetMHCpan not found at: ", netmhcpan_path, "\n",
       "       Please update the path or install NetMHCpan-4.2")
}

cat("  ✓ NetMHCpan found:", netmhcpan_path, "\n")
cat("  ✓ Starting Class I predictions (this may take 10-15 minutes)...\n\n")

# Run NetMHCpan
netmhcpan_cmd <- sprintf(
  "%s -p %s -a %s -BA > %s 2>&1",
  netmhcpan_path,
  "results/netmhcpan_output/peptides_9mer.txt",
  hla_class1_string,
  "results/netmhcpan_output/predictions_class1.txt"
)

system(netmhcpan_cmd)

cat("\n  ✓ NetMHCpan complete!\n")
cat("  ✓ Output: results/netmhcpan_output/predictions_class1.txt\n\n")

# ─────────────────────────────────────────────────────────────────
# 8. Run NetMHCIIpan (Class II) - DRB1 ONLY
# ─────────────────────────────────────────────────────────────────

cat("[8] Running NetMHCIIpan-4.3 predictions (Class II - DRB1 only)...\n")

# Check NetMHCIIpan installation
netmhciipan_path <- "/Users/ranzo85p/netMHCIIpan-4.3/netMHCIIpan"

if (!file.exists(netmhciipan_path)) {
  cat("\n")
  cat("  ⚠ WARNING: NetMHCIIpan not found at:", netmhciipan_path, "\n")
  cat("  ⚠ Skipping Class II predictions\n")
  cat("  ℹ To enable Class II predictions:\n")
  cat("    1. Request NetMHCIIpan from: https://services.healthtech.dtu.dk/services/NetMHCIIpan-4.3/\n")
  cat("    2. Install and update the path above\n")
  cat("    3. Re-run this script\n\n")
  
  class2_available <- FALSE
  
} else {
  
  class2_available <- TRUE
  
  hla_class2_string <- paste(hla_class2, collapse = ",")
  
  cat("  ✓ Peptides:", nrow(unique_peptides_15mer), "\n")
  cat("  ✓ HLA alleles:", length(hla_class2), "(DRB1 only)\n")
  cat("  ✓ Expected predictions:", nrow(unique_peptides_15mer) * length(hla_class2), "\n\n")
  
  cat("  ✓ NetMHCIIpan found:", netmhciipan_path, "\n")
  cat("  ✓ Starting Class II predictions (this may take 15-20 minutes)...\n\n")
  
  # Run NetMHCIIpan with -inptype 1 (peptide mode - CRITICAL FIX!)
  netmhciipan_cmd <- sprintf(
    "%s -inptype 1 -f %s -a %s > %s 2>&1",
    netmhciipan_path,
    "results/netmhcpan_output/peptides_15mer.txt",
    hla_class2_string,
    "results/netmhcpan_output/predictions_class2.txt"
  )
  
  system(netmhciipan_cmd)
  
  cat("\n  ✓ NetMHCIIpan complete!\n")
  cat("  ✓ Output: results/netmhcpan_output/predictions_class2.txt\n\n")
}

# ─────────────────────────────────────────────────────────────────
# 9. Parse NetMHCpan Results (Class I)
# ─────────────────────────────────────────────────────────────────

cat("[9] Parsing NetMHCpan predictions (Class I)...\n")

lines <- readLines("results/netmhcpan_output/predictions_class1.txt")
data_lines <- lines[grepl("^\\s+[0-9]+\\s+HLA-", lines)]

cat("  ✓ Found", length(data_lines), "Class I prediction lines\n")

predictions_class1 <- data.frame()

for (line in data_lines) {
  fields <- unlist(strsplit(trimws(line), "\\s+"))
  
  if (length(fields) < 15) next
  
  predictions_class1 <- rbind(predictions_class1, data.frame(
    peptide = fields[3],
    hla = fields[2],
    affinity_nM = as.numeric(fields[15]),
    percentile = as.numeric(fields[13]),
    binding_level = if(length(fields) >= 17) fields[17] else "",
    hla_class = "I",
    peptide_length = 9,
    stringsAsFactors = FALSE
  ))
}

# Classify strong binders
predictions_class1$is_strong_binder <- 
  predictions_class1$affinity_nM < 500 & predictions_class1$percentile < 2

cat("  ✓ Parsed", nrow(predictions_class1), "Class I predictions\n")
cat("  ✓ Strong binders:", sum(predictions_class1$is_strong_binder), "\n")
cat("  ✓ Hit rate:", 
    round(100 * sum(predictions_class1$is_strong_binder) / nrow(predictions_class1), 1), 
    "%\n\n")

# ─────────────────────────────────────────────────────────────────
# 10. Parse NetMHCIIpan Results (Class II)
# ─────────────────────────────────────────────────────────────────

if (class2_available) {
  
  cat("[10] Parsing NetMHCIIpan predictions (Class II)...\n")
  
  lines <- readLines("results/netmhcpan_output/predictions_class2.txt")
  
  # NetMHCIIpan output format - look for data lines
  # Format: Pos  MHC  Peptide  Of  Core  Core_Rel  Inverted  Identity  Score_EL  %Rank_EL  Exp_Bind  BindLevel
  data_lines <- lines[grepl("^\\s+[0-9]+\\s+DRB1", lines)]
  
  cat("  ✓ Found", length(data_lines), "Class II prediction lines\n")
  
  predictions_class2 <- data.frame()
  
  for (line in data_lines) {
    fields <- unlist(strsplit(trimws(line), "\\s+"))
    
    # NetMHCIIpan-4.3 output has ~12 columns
    if (length(fields) < 10) next
    
    # Try to parse - columns: Pos MHC Peptide Of Core Core_Rel Inverted Identity Score_EL %Rank_EL Exp_Bind BindLevel
    tryCatch({
      
      # For EL mode (default): Score_EL is in column 9, %Rank_EL in column 10
      peptide_seq <- fields[3]
      mhc_allele <- fields[2]
      score_el <- as.numeric(fields[9])
      rank_el <- as.numeric(fields[10])
      
      # Convert rank to "affinity-like" score for consistency
      # Lower rank = stronger binding
      # Approximate affinity from rank (not real IC50, but for comparison)
      pseudo_affinity <- 50000 / (100 - rank_el + 1)  # Rough approximation
      
      predictions_class2 <- rbind(predictions_class2, data.frame(
        peptide = peptide_seq,
        hla = paste0("HLA-", mhc_allele),  # Add HLA- prefix
        affinity_nM = pseudo_affinity,
        percentile = rank_el,
        binding_level = "",
        hla_class = "II",
        peptide_length = 15,
        stringsAsFactors = FALSE
      ))
    }, error = function(e) {
      # Skip lines that can't be parsed
    })
  }
  
  # Classify strong binders (using %Rank threshold)
  # For EL predictions: %Rank < 2% = strong binder
  predictions_class2$is_strong_binder <- predictions_class2$percentile < 2
  
  cat("  ✓ Parsed", nrow(predictions_class2), "Class II predictions\n")
  cat("  ✓ Strong binders:", sum(predictions_class2$is_strong_binder), "\n")
  cat("  ✓ Hit rate:", 
      round(100 * sum(predictions_class2$is_strong_binder) / nrow(predictions_class2), 1), 
      "%\n\n")
  
} else {
  
  cat("[10] Skipping Class II parsing (NetMHCIIpan not available)\n\n")
  
  predictions_class2 <- data.frame()
}

# ─────────────────────────────────────────────────────────────────
# 11. Merge Class I and Class II Predictions
# ─────────────────────────────────────────────────────────────────

cat("[11] Merging Class I and Class II predictions...\n")

# Combine predictions
all_predictions <- bind_rows(predictions_class1, predictions_class2)

cat("  ✓ Total predictions:", nrow(all_predictions), "\n")
cat("    - Class I:", nrow(predictions_class1), "\n")
cat("    - Class II:", nrow(predictions_class2), "\n")
cat("  ✓ Total strong binders:", sum(all_predictions$is_strong_binder), "\n")
cat("    - Class I:", sum(predictions_class1$is_strong_binder), "\n")
cat("    - Class II:", sum(predictions_class2$is_strong_binder), "\n\n")

# ─────────────────────────────────────────────────────────────────
# 12. Merge with Metadata
# ─────────────────────────────────────────────────────────────────

cat("[12] Merging predictions with mutation metadata...\n")

# Combine 9-mer and 15-mer metadata
all_peptides_metadata <- bind_rows(all_9mers, all_15mers)

# Merge predictions with metadata
predictions_full <- all_predictions %>%
  left_join(all_peptides_metadata, 
            by = c("peptide", "hla_class"), 
            relationship = "many-to-many")

write.csv(predictions_full, 
          "data/processed/NetMHCpan_All_Predictions.csv", 
          row.names = FALSE)

cat("  ✓ Saved: data/processed/NetMHCpan_All_Predictions.csv\n")

# Filter strong binders
strong_binders <- predictions_full %>%
  filter(is_strong_binder == TRUE) %>%
  arrange(hla_class, affinity_nM)

write.csv(strong_binders, 
          "data/processed/Strong_Binders_Final.csv", 
          row.names = FALSE)

cat("  ✓ Saved: data/processed/Strong_Binders_Final.csv\n\n")

# ─────────────────────────────────────────────────────────────────
# 13. Summary by Gene and HLA Class
# ─────────────────────────────────────────────────────────────────

cat("[13] Summarizing results...\n\n")

# Summary by gene and HLA class
gene_class_summary <- strong_binders %>%
  group_by(gene, hla_class) %>%
  summarise(
    n_binders = n(),
    n_unique_peptides = n_distinct(peptide),
    n_patients = n_distinct(patient),
    best_affinity = min(affinity_nM, na.rm = TRUE),
    mean_affinity = round(mean(affinity_nM, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(gene, hla_class)

cat("Strong Binders by Gene and HLA Class:\n")
print(gene_class_summary, n = Inf)
cat("\n")

# Overall summary by HLA class
class_summary <- strong_binders %>%
  group_by(hla_class) %>%
  summarise(
    n_binders = n(),
    n_unique_peptides = n_distinct(peptide),
    n_patients = n_distinct(patient),
    n_genes = n_distinct(gene),
    .groups = "drop"
  )

cat("Overall Summary by HLA Class:\n")
print(class_summary)
cat("\n")

# Patient coverage
total_patients <- length(unique(strong_binders$patient))
cat("Patient Coverage:\n")
cat("  Total patients with strong binders:", total_patients, "\n")
cat("  Percentage of cohort:", 
    round(100 * total_patients / length(unique(final_candidates$Tumor_Sample_Barcode)), 1), 
    "%\n\n")

# ─────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat("NEOANTIGEN PREDICTION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Files created:\n")
cat("  • data/processed/protein_sequences.csv\n")
cat("  • results/netmhcpan_output/peptides_9mer.txt\n")
cat("  • results/netmhcpan_output/peptides_15mer.txt\n")
cat("  • results/netmhcpan_output/predictions_class1.txt\n")
if (class2_available) {
  cat("  • results/netmhcpan_output/predictions_class2.txt\n")
}
cat("  • data/processed/NetMHCpan_All_Predictions.csv\n")
cat("  • data/processed/Strong_Binders_Final.csv\n\n")

cat("Summary:\n")
cat("  Total predictions:", nrow(all_predictions), "\n")
cat("    - Class I:", nrow(predictions_class1), "\n")
cat("    - Class II:", nrow(predictions_class2), "\n")
cat("  Strong binders:", sum(all_predictions$is_strong_binder), "\n")
cat("    - Class I:", sum(predictions_class1$is_strong_binder), "\n")
cat("    - Class II:", sum(predictions_class2$is_strong_binder), "\n")
cat("  Genes:", paste(unique(gene_class_summary$gene), collapse = ", "), "\n")
cat("  Patients covered:", total_patients, "\n\n")

if (!class2_available) {
  cat("NOTE: Class II predictions were not run.\n")
  cat("      To enable, install NetMHCIIpan and update the path in this script.\n\n")
}

cat("Next step: Run scripts/06_survival_analysis.R\n\n")

cat("═══════════════════════════════════════════════════════════════\n\n")
