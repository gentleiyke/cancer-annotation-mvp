# Install packages and paths
# install.packages(c("jsonlite", "BiocManager", "dplyr", "purrr", "stringr"))
# BiocManager::install("Biostrings")

library(jsonlite)
library(Biostrings)
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)

# Helper: driver gene list
driver_genes <- c(
  "BRCA1", "BRCA2", "PALB2",
  "CDKN2A", "STK11", "KRAS",
  "EGFR", "TP53", "SMAD4"
)

# Create a tiering function
assign_tier <- function(mutation_classification,
                        afr_frequency,
                        gene) {
  # Normalise inputs
  cls <- ifelse(is.na(mutation_classification), "", tolower(mutation_classification))
  afr <- afr_frequency
  g   <- toupper(gene)
  
  # Define helper flags
  is_pathogenic   <- str_detect(cls, "pathogenic") & !str_detect(cls, "benign")
  is_benign       <- str_detect(cls, "benign")
  is_vus          <- str_detect(cls, "uncertain") | str_detect(cls, "vus")
  is_driver_gene  <- g %in% driver_genes
  
  afr_is_known    <- !is.na(afr)
  afr_rare        <- afr_is_known & afr < 0.01
  afr_common      <- afr_is_known & afr >= 0.01
  
  # Default output
  tier        <- "Tier 3 – Research interest"
  explanation <- c()
  
  # --- Tier 1: Actionable / high clinical relevance ---
  # Pathogenic/likely pathogenic, in driver gene, rare/non-African
  if (is_driver_gene && is_pathogenic && (afr_rare || !afr_is_known)) {
    tier <- "Tier 1 – Actionable / high clinical relevance"
    explanation <- c(explanation,
                     "Pathogenic/likely pathogenic variant in canonical cancer driver gene.",
                     if (afr_rare) "Rare in African populations (AFR < 1%)." else "No AFR frequency available (assumed rare).")
  }
  
  # --- Tier 2: Clinically relevant but less specific ---
  # Pathogenic in non-driver or is common-ish but still pathogenic
  else if (is_pathogenic) {
    tier <- "Tier 2 – Clinically relevant"
    explanation <- c(explanation,
                     "Pathogenic/likely pathogenic variant.",
                     if (is_driver_gene) "Driver gene, but AFR frequency suggests possible germline or founder effect."
                     else "Non-core driver gene; still clinically important.")
  }
  
  # --- Tier 3: Research interest ---
  # VUS or unclassified variant that is rare in AFR, in driver gene
  else if ((is_vus || (!is_pathogenic && !is_benign && cls != "")) &&
           is_driver_gene && (afr_rare || !afr_is_known)) {
    tier <- "Tier 3 – Research interest"
    explanation <- c(explanation,
                     "Variant of uncertain/other significance in driver gene.",
                     if (afr_rare) "Rare in African populations (AFR < 1%)." else "No AFR frequency available.")
  }
  
  # --- Tier 4: Likely passenger / benign ---
  # Benign / likely benign, or clearly common in AFR
  if (is_benign || afr_common || (!is_pathogenic && afr_common)) {
    tier <- "Tier 4 – Likely passenger / benign"
    explanation <- c(explanation,
                     if (is_benign) "Clinically classified as benign/likely benign." else NULL,
                     if (afr_common) "Common in African populations (AFR ≥ 1%)." else NULL)
  }
  
  # Combine explanation
  explanation_text <- if (length(explanation) == 0) NA_character_ else paste(unique(explanation), collapse = " ")
  
  list(
    tier = tier,
    tier_explanation = explanation_text
  )
}

# Adjust structure if different
DATA_ROOT <- "gene_panel"
FASTA_DIR <- file.path(DATA_ROOT, "fasta")
MUT_DIR   <- file.path(DATA_ROOT, "mutations")
AFR_DIR   <- file.path(DATA_ROOT, "afr_frequencies")
META_FILE <- file.path(DATA_ROOT, "metadata", "gene_metadata.json")

# Load FASTA sequences
# Return a named list: gene -> DNAStringSet (usually length 1)
load_fasta_panel <- function(fasta_dir = FASTA_DIR) {
  fasta_files <- list.files(fasta_dir, pattern = "\\.fasta$", full.names = TRUE)
  
  seqs <- lapply(fasta_files, function(f) {
    readDNAStringSet(f)
  })
  names(seqs) <- basename(fasta_files) |>
    str_replace("\\.fasta$", "")  # e.g. "BRCA1"
  
  seqs
}

# Load sequences
gene_fasta <- load_fasta_panel()

# Extract FASTA headers for each gene
gene_headers <- sapply(gene_fasta, function(x) {
  names(x)
})

# Show gene header mapping
tibble::tibble(
  gene = names(gene_headers),
  header = unname(gene_headers)
)

# Load single gene
tp53_seq <- as.character(gene_fasta[["TP53"]][[1]])
# View first 50 bases
substr(tp53_seq, 1, 50)  

# Load mutation JSONs into a data frame
load_mutations <- function(mut_dir = MUT_DIR) {
  mut_files <- list.files(mut_dir, pattern = "_mutations\\.json$", full.names = TRUE)
  
  mut_list <- lapply(mut_files, function(path) {
    x <- fromJSON(path)
    gene <- x$gene
    if (length(x$mutations) == 0) {
      return(NULL)
    }
    as_tibble(x$mutations) |>
      mutate(gene = gene, .before = 1)
  })
  
  bind_rows(mut_list)
}

mut_df <- load_mutations()
mut_df

# Load AFR frequency JSONs
load_afr <- function(afr_dir = AFR_DIR) {
  afr_files <- list.files(afr_dir, pattern = "_afr_frequencies\\.json$", full.names = TRUE)
  
  afr_list <- lapply(afr_files, function(path) {
    x <- fromJSON(path)
    gene <- x$gene
    if (length(x$frequencies) == 0) {
      return(NULL)
    }
    as_tibble(x$frequencies) |>
      mutate(gene = gene, .before = 1)
  })
  
  bind_rows(afr_list)
}

afr_df <- load_afr()
afr_df

# Load gene metadata
meta_list <- fromJSON(META_FILE)
# meta_list is a named list: BRCA1, BRCA2, ...
names(meta_list)

meta_df <- bind_rows(lapply(names(meta_list), function(g) {
  tibble(
    gene = g,
    chromosome = meta_list[[g]]$chromosome,
    gene_length = meta_list[[g]]$gene_length,
    function_summary = meta_list[[g]]$`function`,
    clinical_relevance = meta_list[[g]]$clinical_relevance,
    major_cancer_associations = paste(meta_list[[g]]$major_cancer_associations, collapse = ", "),
    hotspots = paste(meta_list[[g]]$hotspots, collapse = ", ")
  )
}))
meta_df

# Simple annotation function
annotate_variant <- function(gene, position, ref, alt,
                             mut_df, afr_df, meta_df) {
  position <- as.integer(position)
  ref <- toupper(ref)
  alt <- toupper(alt)
  
  # 1. Known mutation lookup
  mut_match <- mut_df |>
    filter(gene == !!gene,
           !is.na(position),
           position == !!position,
           (is.na(ref) | ref == !!ref),
           (is.na(alt) | alt == !!alt))
  
  if (nrow(mut_match) == 0) {
    mut_status <- "not_in_curated_list"
    mut_details <- NA
    mut_classification <- NA
  } else {
    mm <- mut_match[1, ]
    mut_status <- "matched_known_variant"
    mut_classification <- mm$classification
    mut_details <- paste0(
      "type=", mm$type,
      "; protein_change=", mm$protein_change,
      "; classification=", mm$classification
    )
  }
  
  # 2. AFR frequency lookup
  afr_match <- afr_df |>
    filter(gene == !!gene,
           !is.na(position),
           position == !!position,
           (is.na(ref) | ref == !!ref),
           (is.na(alt) | alt == !!alt))
  
  if (nrow(afr_match) == 0) {
    afr_freq <- NA_real_
    afr_note <- "no_AFR_record"
  } else {
    af <- afr_match$frequency[1]
    afr_freq <- af
    afr_note <- if (is.na(af)) {
      "AFR_frequency_missing"
    } else if (af == 0) {
      "AFR_frequency=0 (not observed in AFR dataset)"
    } else {
      paste0("AFR_frequency=", af)
    }
  }
  
  # 3. Gene metadata
  md <- meta_df |>
    filter(gene == !!gene) |>
    slice(1)
  
  # 4. Tiering
  tier_info <- assign_tier(
    mutation_classification = mut_classification,
    afr_frequency = afr_freq,
    gene = gene
  )
  
  tibble(
    gene = gene,
    position = position,
    ref = ref,
    alt = alt,
    mutation_status = mut_status,
    mutation_classification = mut_classification,
    mutation_details = mut_details,
    afr_frequency = afr_freq,
    afr_comment = afr_note,
    tier = tier_info$tier,
    tier_explanation = tier_info$tier_explanation,
    gene_chromosome = md$chromosome,
    gene_length = md$gene_length,
    gene_function = md$function_summary,
    gene_clinical_relevance = md$clinical_relevance,
    gene_cancer_associations = md$major_cancer_associations,
    gene_hotspots = md$hotspots
  )
}

# Run Example – OncoMatch-style demo
# We created a small example "variant list" to mimic what would come out
# of a variant calling pipeline in the full OncoMatch system.
variants <- tribble(
  ~gene,  ~position,  ~ref, ~alt,
  "BRCA1", 43045711, "C",   "G",
  "TP53",  7670669,  "C",   "A",
  "EGFR",  55191822, "T",   "G"
)

# Conceptually, this is the "annotation engine" step of OncoMatch where the
# Tumour/normal pipeline → VCF/variant table → Annotation engine
annotated <- variants |>
  pmap_dfr(~ annotate_variant(..1, ..2, ..3, ..4, mut_df, afr_df, meta_df))

annotated

# Single-variant deep dive (BRCA1)
# Here we manually call annotate_variant() for one BRCA1 variant.
# Conceptually, this is like:
#   "Given a specific mutation found in this patient's tumour,
#    tell me what it means in the context of our panel and African populations."
example <- annotate_variant(
  gene = "TP53",
  position = 7670669,
  ref = "C",
  alt = "A",
  mut_df = mut_df,
  afr_df = afr_df,
  meta_df = meta_df
)

# And we then select only the most clinically-relevant columns to show:
# - gene, position: which variant?
# - mutation_classification: how is it classified (pathogenic, benign, etc.)?
# - afr_frequency: how common is it in African populations?
# - tier: which "What Matters" bucket it falls into (Tier 1–4)?
# - tier_explanation: short text explaining why it got that tier.
example %>% select(gene, position, mutation_classification, afr_frequency, tier, tier_explanation)

