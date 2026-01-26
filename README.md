### Cancer Annotation MVP
A minimal R-based pipeline for annotating cancer variants and assigning clinical-style tiers using curated gene panel data. The repository includes both a scripted annotation engine and a polished HTML report demonstrating how results can be summarised and reviewed.

![Report screenshot](report_oncomatch_demo.png)


| Input                                            | Output                                |
| ------------------------------------------------ | ------------------------------------- |
| Variant table (`gene`, `position`, `ref`, `alt`) | Annotated variants table with tiering |
| Curated gene panel (`gene_panel/`)               | Tier summary + deep-dive examples     |
| R Markdown report                                | HTML report with tables and plots     |

#### Repos Structure
```
├── cancer-annotation-mvp.R          # Core annotation logic
├── report_oncomatch_demo.Rmd        # Example HTML report
├── gene_panel/                      # Curated gene reference data
├── examples/
│   ├── input_variants.csv
│   ├── output_annotated_variants.csv
│   └── output_summary.json
└── README.md
```

See `gene_panel/docs/DATA_SOURCES_AND_METHODS.txt` for detailed data sources and methods

#### Inputs
#### Variant Input Table
The annotation engine expects a simple variant table with the following required columns:
| Column     | Description                |
| ---------- | -------------------------- |
| `gene`     | Gene symbol                |
| `position` | Genomic position (integer) |
| `ref`      | Reference allele           |
| `alt`      | Alternate allele           |

Example (examples/input_variants.csv)

#### Gene panel reference
Curated gene and variant-level metadata under:
```
gene_panel/
```
This data is used internally by the annotation functions to determine:
- driver status
- mutation classification
- allele frequency 
- tier assignment logic

#### Outputs
#### Annotated variants table
Running the annotation script produces a tibble with the following key columns:
| Column                    | Description                   |
| ------------------------- | ----------------------------- |
| `gene`                    | Gene symbol                   |
| `position`                | Genomic position              |
| `ref`                     | Reference allele              |
| `alt`                     | Alternate allele              |
| `mutation_status`         | Known / unknown pathogenicity |
| `mutation_classification` | Missense, nonsense, etc.      |
| `afr_frequency`           | Population allele frequency   |
| `tier`                    | Assigned clinical tier        |
| `tier_explanation`        | Human-readable rationale      |

#### Pipeline Overview
```
Variant table (gene, position, ref, alt)
                │
                ▼
     Annotation + tiering engine
          (cancer-annotation-mvp.R)
                │
                ▼
      Annotated variants table
                │
                ├── CSV / JSON outputs
                ▼
        HTML report
 (report_oncomatch_demo.Rmd)
```

#### How to run the annotation demo:
1. Clone this repo:

   ```bash
   git clone https://github.com/gentleiyke/cancer-annotation-mvp.git
   cd cancer-annotation-mvp
2. Open and run report_oncomatch_demo.Rmd 
