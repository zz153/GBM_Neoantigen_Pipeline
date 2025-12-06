#!/usr/bin/env Rscript
################################################################################
# 05_neoantigen_prediction.R
# 
# Predicts HLA-binding neoantigens using NetMHCpan-4.2
# 
# Steps:
#   1. Fetch protein sequences from UniProt
#   2. Generate 21-mer peptides (mutation centered)
#   3. Extract 9-mer peptides
#   4. Run NetMHCpan predictions
#   5. Parse and filter results
# 
# Inputs:
#   - data/processed/Evolution_Resistant_Candidates.csv
# 
# Outputs:
#   - data/processed/protein_sequences.csv
#   - results/netmhcpan_output/peptides.txt
#   - results/netmhcpan_output/peptides_metadata.csv
#   - results/netmhcpan_output/predictions.txt
#   - data/processed/NetMHCpan_All_Predictions.csv
#   - data/processed/Strong_Binders_Final.csv
################################################################################

cat("\n═══════════════════════════════════════════════════════════════\n")
cat("STEP 5: NEOANTIGEN PREDICTION (NetMHCpan)\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

# ─────────────────────────────────────────────────────────────────
# Setup
# ─────────────────────────────────────────────────────────────────

library(dplyr)
library(httr)
library(jsonlite)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)
dir.create("results/netmhcpan_output", recursive = TRUE, showWarnings = FALSE)

# ─────────────────────────────────────────────────────────────────
# 1. Load Clonal Mutations
# ─────────────────────────────────────────────────────────────────

cat("[1] Loading clonal mutations...\n")

final_candidates <- read.csv("data/processed/Evolution_Resistant_Candidates.csv")

cat("  ✓ Loaded", nrow(final_candidates), "clonal mutations\n")
cat("  ✓ Genes:", paste(unique(final_candidates$Hugo_Symbol), collapse = ", "), "\n")
cat("  ✓ Patients:", length(unique(final_candidates$patient_id)), "\n\n")

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
# 4. Generate 9-mer Peptides
# ─────────────────────────────────────────────────────────────────

cat("[4] Generating 9-mer peptides for NetMHCpan...\n")

# Function to generate 9-mers containing mutation
generate_9mers <- function(peptide_21mer, mt_aa) {
  if (is.na(peptide_21mer) || is.na(mt_aa)) return(list())
  
  mut_pos <- 11  # Center of 21-mer
  peptides_9mer <- c()
  
  for (i in 1:13) {
    start <- i
    end <- i + 8
    peptide <- substr(peptide_21mer, start, end)
    
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
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("  ✓ Generated", nrow(all_9mers), "9-mer peptides\n")

# Remove duplicates and stop codons
unique_peptides <- all_9mers %>%
  filter(!grepl("\\*", peptide)) %>%
  distinct(peptide, .keep_all = TRUE)

cat("  ✓ Unique clean peptides:", nrow(unique_peptides), "\n")

# Save
writeLines(unique_peptides$peptide, "results/netmhcpan_output/peptides.txt")
write.csv(all_9mers, "results/netmhcpan_output/peptides_metadata.csv", row.names = FALSE)

cat("  ✓ Saved: results/netmhcpan_output/peptides.txt\n")
cat("  ✓ Saved: results/netmhcpan_output/peptides_metadata.csv\n\n")

# ─────────────────────────────────────────────────────────────────
# 5. Run NetMHCpan
# ─────────────────────────────────────────────────────────────────

cat("[5] Running NetMHCpan-4.2 predictions...\n")

# HLA alleles (common panel, ~70% coverage)
hla_alleles <- c(
  "HLA-A*02:01", "HLA-A*01:01", "HLA-A*24:02", "HLA-A*03:01",
  "HLA-B*07:02", "HLA-B*08:01", "HLA-B*44:02",
  "HLA-C*07:02", "HLA-C*07:01"
)

hla_string <- paste(gsub("\\*", "", hla_alleles), collapse = ",")

cat("  ✓ Peptides:", nrow(unique_peptides), "\n")
cat("  ✓ HLA alleles:", length(hla_alleles), "\n")
cat("  ✓ Expected predictions:", nrow(unique_peptides) * length(hla_alleles), "\n\n")

# Check NetMHCpan installation
netmhcpan_path <- "/Users/ranzo85p/netMHCpan-4.2/netMHCpan"

if (!file.exists(netmhcpan_path)) {
  stop("ERROR: NetMHCpan not found at: ", netmhcpan_path, "\n",
       "       Please update the path or install NetMHCpan-4.2")
}

cat("  ✓ NetMHCpan found:", netmhcpan_path, "\n")
cat("  ✓ Starting predictions (this may take 10-15 minutes)...\n\n")

# Run NetMHCpan
netmhcpan_cmd <- sprintf(
  "%s -p %s -a %s -BA > %s 2>&1",
  netmhcpan_path,
  "results/netmhcpan_output/peptides.txt",
  hla_string,
  "results/netmhcpan_output/predictions.txt"
)

system(netmhcpan_cmd)

cat("\n  ✓ NetMHCpan complete!\n")
cat("  ✓ Output: results/netmhcpan_output/predictions.txt\n\n")

# ─────────────────────────────────────────────────────────────────
# 6. Parse NetMHCpan Results
# ─────────────────────────────────────────────────────────────────

cat("[6] Parsing NetMHCpan predictions...\n")

lines <- readLines("results/netmhcpan_output/predictions.txt")
data_lines <- lines[grepl("^\\s+[0-9]+\\s+HLA-", lines)]

cat("  ✓ Found", length(data_lines), "prediction lines\n")

predictions <- data.frame()

for (line in data_lines) {
  fields <- unlist(strsplit(trimws(line), "\\s+"))
  
  if (length(fields) < 15) next
  
  predictions <- rbind(predictions, data.frame(
    peptide = fields[3],
    hla = fields[2],
    affinity_nM = as.numeric(fields[15]),
    percentile = as.numeric(fields[13]),
    binding_level = if(length(fields) >= 17) fields[17] else "",
    stringsAsFactors = FALSE
  ))
}

# Classify strong binders
predictions$is_strong_binder <- predictions$affinity_nM < 500 & predictions$percentile < 2

cat("  ✓ Parsed", nrow(predictions), "predictions\n")
cat("  ✓ Strong binders (IC50 < 500 nM):", sum(predictions$is_strong_binder), "\n")
cat("  ✓ Percentage strong:", round(100 * sum(predictions$is_strong_binder) / nrow(predictions), 1), "%\n\n")

# ─────────────────────────────────────────────────────────────────
# 7. Merge with Metadata
# ─────────────────────────────────────────────────────────────────

cat("[7] Merging predictions with mutation metadata...\n")

predictions_full <- predictions %>%
  left_join(all_9mers, by = "peptide", relationship = "many-to-many")

write.csv(predictions_full, 
          "data/processed/NetMHCpan_All_Predictions.csv", 
          row.names = FALSE)

cat("  ✓ Saved: data/processed/NetMHCpan_All_Predictions.csv\n")

# Filter strong binders
strong_binders <- predictions_full %>%
  filter(is_strong_binder == TRUE) %>%
  arrange(affinity_nM)

write.csv(strong_binders, 
          "data/processed/Strong_Binders_Final.csv", 
          row.names = FALSE)

cat("  ✓ Saved: data/processed/Strong_Binders_Final.csv\n\n")

# ─────────────────────────────────────────────────────────────────
# 8. Summary by Gene
# ─────────────────────────────────────────────────────────────────

cat("[8] Summarizing results by gene...\n")

gene_summary <- strong_binders %>%
  group_by(gene) %>%
  summarise(
    n_binders = n(),
    n_unique_peptides = n_distinct(peptide),
    n_patients = n_distinct(patient),
    best_affinity = min(affinity_nM, na.rm = TRUE),
    mean_affinity = round(mean(affinity_nM, na.rm = TRUE), 2),
    .groups = "drop"
  ) %>%
  arrange(desc(n_binders))

cat("\nStrong Binders by Gene:\n")
print(gene_summary, row.names = FALSE)

cat("\n")

# ─────────────────────────────────────────────────────────────────
# Summary
# ─────────────────────────────────────────────────────────────────

cat("═══════════════════════════════════════════════════════════════\n")
cat("NEOANTIGEN PREDICTION COMPLETE\n")
cat("═══════════════════════════════════════════════════════════════\n\n")

cat("Files created:\n")
cat("  • data/processed/protein_sequences.csv\n")
cat("  • results/netmhcpan_output/peptides.txt\n")
cat("  • results/netmhcpan_output/peptides_metadata.csv\n")
cat("  • results/netmhcpan_output/predictions.txt\n")
cat("  • data/processed/NetMHCpan_All_Predictions.csv\n")
cat("  • data/processed/Strong_Binders_Final.csv\n\n")

cat("Summary:\n")
cat("  Total predictions:", nrow(predictions), "\n")
cat("  Strong binders:", sum(predictions$is_strong_binder), "\n")
cat("  Genes:", paste(unique(gene_summary$gene), collapse = ", "), "\n")
cat("  Patients covered:", length(unique(strong_binders$patient)), "\n\n")

cat("Next step: Run scripts/06_survival_analysis.R\n\n")

cat("═══════════════════════════════════════════════════════════════\n\n")