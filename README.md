# Phenotype Matching Tool

This repository contains tools to map phenotypes between FinnGen endpoints and PheCode/ICD-10 datasets using ICD-10 codes as an intermediary.

## Overview

The matching algorithm finds the best-matching target phenotype for each source phenotype by calculating ICD-10 code similarity:

1. **FinnGen ICD-10 codes** are defined as regular expressions matching sets of ICD-10 codes
2. **PheCodes** are mapped to specific ICD-10 codes via lookup tables
3. **Similarity score** is calculated using the Jaccard index: `|intersection| / |union|` of the ICD-10 code sets

The tool supports bidirectional matching:

- Map FinnGen endpoints → PheCodes/UKBB phenotypes
- Map PheCodes/UKBB phenotypes → FinnGen endpoints

The tool was designed originally to match UKBB Phecodes to FinnGen endpoints, but can be easily adapted to other phenotype definitions with similar structure, e.g. phecodeX: https://doi.org/10.1093/bioinformatics/btad655

## Features

- Bidirectional phenotype matching
- Handles FinnGen endpoint inclusion dependencies
- Configurable file separators and column names
- Optional filtering by number of mapped ICD-10 codes
- Returns top N alternative matches for each phenotype
- Supports conditional exclusion criteria for FinnGen endpoints

## Installation

### Requirements
- Python 3
- pandas
- numpy

### Install dependencies

```bash
pip3 install pandas numpy
```

## Usage

### Basic Syntax

```bash
python Scripts/phenomatch.py \
  --main-table [phecode|finngen] \
  --out OUTPUT_FILE.tsv \
  --phecode-source PHECODE_FILE \
  --fg-source FINNGEN_FILE \
  --map-source MAPPING_FILE \
  [additional options...]
```

### Required Arguments

| Argument | Description |
|----------|-------------|
| `--main-table` | Join direction: `phecode` (map FG→PheCodes) or `finngen` (map PheCodes→FG) |
| `--out` | Output filename |

#### PheCode/ICD-10 Data
| Argument | Description |
|----------|-------------|
| `--phecode-source` | Path to PheCode/ICD-10 file with phenotypes to match |
| `--pheno-pheno-col` | Column name for phenotype identifiers |
| `--pheno-type-col` | Column name for phenotype type (`phecode` or `icd10`) |

#### Mapping Data
| Argument | Description |
|----------|-------------|
| `--map-source` | Path to PheCode↔ICD-10 mapping file |
| `--map-pheno-col` | PheCode column name in mapping file |
| `--map-icd-col` | ICD-10 column name in mapping file |

#### FinnGen Data
| Argument | Description |
|----------|-------------|
| `--fg-source` | Path to FinnGen endpoints file |
| `--fg-pheno-col` | Phenotype column name in FinnGen file |
| `--fg-icd-col` | ICD-10 column name(s) in FinnGen file (space-separated if multiple) |
| `--fg-inc-col` | Column listing included FinnGen endpoints (pipe-separated) |

### Optional Arguments

| Argument | Default | Description |
|----------|---------|-------------|
| `--other-hits-n` | 5 | Number of alternative matches to include |
| `--pheno-sep` | `\t` | Separator for PheCode file |
| `--map-sep` | `\t` | Separator for mapping file |
| `--fg-sep` | `\t` | Separator for FinnGen file |
| `--map-filter` | 0 | Exclude PheCodes mapping to more than N ICD-10 codes (0=no filter) |
| `--icd-code-source` | - | Custom ICD-10 code list file (auto-extracted from mapping if not provided) |
| `--fg-cond-col` | - | Conditional columns in FinnGen file (excludes endpoints with non-empty values) |
| `--fg-icd-excl-col` | - | ICD-10 exclusion columns in FinnGen file |

### Example Commands

#### Map FinnGen endpoints to UKBB/PheCode phenotypes

Each output row contains a PheCode/UKBB endpoint with its best-matching FinnGen endpoint:

```bash
python Scripts/phenomatch.py \
  --main-table phecode \
  --phecode-source data/phenos_full_ukbb_gwas_exome_run.tsv \
  --fg-source data/finngen_R4_endpoints_utf8.tsv \
  --map-source data/Phecode_map_v1_2_icd10_beta.csv \
  --map-sep "," \
  --pheno-pheno-col pheno \
  --fg-pheno-col NAME \
  --map-pheno-col PHECODE \
  --fg-inc-col INCLUDE \
  --map-icd-col ICD10 \
  --fg-icd-col HD_ICD_10 \
  --pheno-type-col data_type \
  --out phecode_to_finngen_matches.tsv
```

#### Map UKBB/PheCode phenotypes to FinnGen endpoints

Each output row contains a FinnGen endpoint with its best-matching PheCode/UKBB endpoint:

```bash
python Scripts/phenomatch.py \
  --main-table finngen \
  --phecode-source data/phenos_full_ukbb_gwas_exome_run.tsv \
  --fg-source data/finngen_R4_endpoints_utf8.tsv \
  --map-source data/Phecode_map_v1_2_icd10_beta.csv \
  --map-sep "," \
  --pheno-pheno-col pheno \
  --fg-pheno-col NAME \
  --map-pheno-col PHECODE \
  --fg-inc-col INCLUDE \
  --map-icd-col ICD10 \
  --fg-icd-col HD_ICD_10 \
  --pheno-type-col data_type \
  --out finngen_to_phecode_matches.tsv
```

## Input Data Requirements

### PheCode/ICD-10 File
- Must contain a phenotype identifier column
- Must contain a phenotype type column with values: `phecode` or `icd10`
- Tab-separated by default (configurable)

### PheCode↔ICD-10 Mapping File
- Many-to-many mapping format (one mapping per line)
- Must contain PheCode column and ICD-10 column
- Tab-separated by default (configurable)

### FinnGen Endpoints File
- Must contain phenotype name column
- Must contain inclusion column (pipe-separated list of included phenotypes)
- Must contain one or more ICD-10 regex columns
- Tab-separated by default (configurable)
- Should be UTF-8 encoded

## Data Sources

- **PheCode-ICD10 mapping**: https://phewascatalog.org/phewas/_w_3494271dba2742b1a77cfc7926d714b9/data/Phecode_map_v1_2_icd9_icd10cm.csv.zip
- **FinnGen endpoints**: https://www.finngen.fi/sites/default/files/inline-files/FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx
- **UKBB phenotypes**: (https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8/edit#gid=903887429)

## Output Format

The output TSV file contains:
- Main phenotype identifier
- Best matching phenotype from the other dataset
- Similarity score (0-1)
- ICD-10 codes for both phenotypes
- Regular expressions used for matching
- Additional matching information
- Alternative matches (top N)

## Project Structure

```
.
├── Scripts/
│   ├── phenomatch.py      # Main matching script
│   ├── data_cleaning.py   # Data preprocessing functions
│   ├── join.py           # Matching algorithm implementation
│   ├── tree.py           # Dependency tree handling
│   ├── constants.py      # Configuration constants
│   └── progress.py       # Progress bar utilities
└── data/                 # Input data files
```

## Acknowledgements

- Tuomo Kiiskinen for the original matching algorithm
- Aki Havulinna & Tuomo Kiiskinen, clinical expert groups & others at FIMM & THL for the FinnGen phenotype definitions
