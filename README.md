# Cancer Annotation MVP
## Data

All reference data for the gene panel lives in:

`gene_panel/`

Contents:
- `fasta/` – GRCh38 reference sequences for:
  BRCA1, BRCA2, PALB2, CDKN2A, STK11, KRAS, EGFR, TP53, SMAD4.
- `mutations/` – Per-gene JSON files of curated pathogenic/benign variants, derived from ClinVar.
- `afr_frequencies/` – Per-gene JSON files with African allele frequencies (gnomAD AFR via GeneBe).
- `metadata/` – Gene-level metadata (chromosome, length, function, clinical relevance, hotspots).
- `docs/DATA_SOURCES_AND_METHODS.txt` – Full description of how the data was collected and processed.

See `gene_panel/docs/DATA_SOURCES_AND_METHODS.txt` for detailed data sources and methods.

**How to run the annotation demo:**

1. Clone this repo:

   ```bash
   git clone https://github.com/gentleiyke/cancer-annotation-mvp.git
   cd cancer-annotation-mvp
2. Open cancer-annotation-mvp.R in RStudio (or run it in R):
   - Ensure your working directory is the project root.
   - Run the script top to bottom.

