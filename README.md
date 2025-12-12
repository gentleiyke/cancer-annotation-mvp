# Cancer Annotation MVP
Gene Panel with AFR frequencies
----------------------------------------------------------
Panel genes:
  BRCA1, BRCA2, PALB2, CDKN2A, STK11, KRAS, EGFR, TP53, SMAD4

Reference sequences (FASTA):
- FASTA files in the 'fasta/' directory were provided and are based on
  Homo sapiens gene sequences (GRCh38 build) with user-defined headers. No changes
  were made to the original sequences when preparing this package.

Known cancer mutation lists (mutations/*.json):
- Per-gene JSON files are derived from the curated table in 'known_cancer_mutation.json'.
- For each variant, we captured:
    * Gene
    * GRCh38 genomic position (parsed from the 'Genomic (GRCh38)' field)
    * Reference and alternate allele (parsed from the cDNA notation where possible)
    * Variant type (missense, nonsense, frameshift, synonymous, etc.)
    * Protein change (HGVS protein)
    * Clinical classification (e.g., pathogenic, likely benign)
    * Notes including source (ClinVar) and associated condition(s).
- Genes without entries in the original table have an empty 'mutations' list and a
  note indicating that no curated mutations were provided for that gene.

African population variant frequencies (afr_frequencies/*.json):
- Per-gene JSON files now include SNP frequencies parsed from 'african_snp_frequencies.json',
  which in turn was derived from gnomAD AFR via GeneBe.
- Each entry contains:
    * position: GRCh38 genomic coordinate parsed from 'genomic position' (e.g., 'chr17:43124023').
    * ref: reference allele from 'ref => alt'.
    * alt: alternate allele from 'ref => alt'.
    * frequency: AFR allele frequency as provided (0 indicates not observed in AFR in this dataset).
    * source: the original AFR dataset description (e.g., 'gnomAD AFR (via GeneBe)').
- Files are named '<GENE>_afr_frequencies.json' and are directly usable by the annotation pipeline.

Gene metadata (metadata/gene_metadata.json):
- For each gene, metadata includes chromosome, gene_length (calculated from the
  provided FASTA sequence), function summary, known hotspot regions or codons,
  clinical relevance, and major cancer associations.
- Descriptions are based on standard knowledge of these canonical cancer genes
  (e.g., BRCA1/2 homologous recombination repair, KRAS as a RAS/MAPK oncogene,
  TP53 as a master tumor suppressor, EGFR as a druggable receptor tyrosine kinase).

Intended use:
- These deliverables are designed for as an MVP of a cancer gene annotation
  engine (Oncomatch) focusing on a targeted gene panel. They provide:
    * Reference gene sequences (FASTA).
    * Curated pathogenic/benign variants (JSON per gene).
    * African population allele frequencies for common SNPs in each gene.
    * Gene-level metadata for reporting and interpretation.
- Together, these files can be used to test variant-calling, annotation and
  reporting components of the pipeline.

