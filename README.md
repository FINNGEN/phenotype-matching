# Map phenotypes using ICD-10 as intermediary

This repository contains scripts to map FG endpoints to PheCodes.

The mapping is done by finding the target phenotype that best matches the source phenotype, by ICD10 similarity.

The similarity is calculated as follows:

FinnGen ICD10 codes are defined as a regular expression that matches some set of ICD-10 codes. As for Phecodes, the PheCode is enumerated for many ICD-10 codes.

The similarity is then calculated by counting, how many of a phecode's icd codes were matched by the FG regex. This is then divided by the total number of phecode icd codes.

## Data sources

Phecode-ICD10 mapping file: https://phewascatalog.org/phecodes_icd10  
FG endpoints: http://results.finngen.fi/static/FINNGEN_ENDPOINTS_DF4_V2.zip (needs to be converted to utf-8 if redonwloaded)  
UKBB endpoint listing: ??  

## Installation

Requires python 3 with packages `numpy` and `pandas`

Install using pip:

```
pip3 install pandas numpy
```

## Usage

Map FG endpoints to UKBB endpoints (each row contains UKBB/PheCode endpoint, and has the best matching FG endpoint on a column):
```
Scripts/phecode_to_fg.py --phecode-source data/phenos_full_ukbb_gwas_exome_run.tsv --fg-source data/finngen_R4_endpoints_utf8.tsv --map-source data/Phecode_map_v1_2_icd10_beta.csv  --pheno-col-phe pheno --pheno-col-fg NAME --pheno-col-map PHECODE --icd-col-map ICD10 --icd-col-fg HD_ICD_10 --phenotype-type-col data_type --out OUTPUT_FILE.tsv
```

Map UKBB endpoints to FG endpoints (each row contains one FG endpoint, and has the best matching UKBB endpoint on a column):
```
Scripts/fg_to_phecode.py --source1 data/finngen_R4_endpoints_utf8.tsv --source2 data/phecode_map_icd10.tsv --pheno-col-1 NAME --pheno-col-2 PHECODE --icd-col-1 HD_ICD_10 --icd-col-2 ICD10 --out OUTPUT_FILE.tsv
```
