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
