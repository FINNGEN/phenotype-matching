# Map phenotypes using ICD-10 as intermediary

This repository contains scripts to map FG endpoints to PheCodes.

The mapping is done by finding the target phenotype that best matches the source phenotype, by ICD10 similarity.

The similarity is calculated as follows:

FinnGen ICD10 codes are defined as a regular expression that matches some set of ICD-10 codes. As for Phecodes, the PheCode is enumerated for many ICD-10 codes.

The FinnGen matching ICD10 codes are acquired by enumerating the ICD10 codes that match a FinnGen endpoint's ICD10 regex.

The similarity score is then calculated by dividing the size of the FG ICD10 codes' and Phecode ICD10 codes' intersection by the size of their union.

## Data sources

Phecode-ICD10 mapping file: https://phewascatalog.org/phecodes_icd10  
FG endpoints: http://results.finngen.fi/static/FINNGEN_ENDPOINTS_DF4_V2.zip (needs to be converted to utf-8 if redonwloaded)  
UKBB endpoint listing: ??  

## Input data

### Phecode/ICD10 endpoints
The file should contain a column for the Phecode/ICD10 code, and a column for its type, which can be either 'phecode' or 'icd10'. The user is responsible for filtering and transforming the input data. File separator is configurable.

### Phecode<->ICD10 mapping file
The file should be a many-to-many mapping for ICD10<->Phecode endpoints, with one mapping per line. The file should have a column for the phecodes, as well as one for the icd10s. File separator is configurable.

### FinnGen endpoint file
The file should have a column for phenotype name, a column for phenotype inclusion (a column that tells which other phenotypes are included for that phenotype), and column(s) for matching icd10 codes. File separator is configurable.

## Installation

Requires python 3 with packages `numpy` and `pandas`

Install using pip:

```
pip3 install pandas numpy
```

## Usage
Map FG endpoints to UKBB endpoints (each row contains UKBB/PheCode endpoint, and has the best matching FG endpoint on a column):
```
Scripts/phenomatch.py --main-table phecode --phecode-source data/phenos_full_ukbb_gwas_exome_run.tsv --fg-source data/finngen_R4_endpoints_utf8.tsv --map-source data/Phecode_map_v1_2_icd10_beta.csv --map-sep "," --pheno-pheno-col pheno --fg-pheno-col NAME --map-pheno-col PHECODE --fg-inc-col INCLUDE --map-icd-col ICD10 --fg-icd-col HD_ICD_10 --pheno-type-col data_type --out OUTPUT_FILE.tsv
```
Map UKBB endpoints to FG endpoints (each row contains one FG endpoint, and has the best matching UKBB endpoint on a column):
```
Scripts/phenomatch.py --main-table finngen --phecode-source data/phenos_full_ukbb_gwas_exome_run.tsv --fg-source data/finngen_R4_endpoints_utf8.tsv --map-source data/Phecode_map_v1_2_icd10_beta.csv --map-sep "," --pheno-pheno-col pheno --fg-pheno-col NAME --map-pheno-col PHECODE --fg-inc-col INCLUDE --map-icd-col ICD10 --fg-icd-col HD_ICD_10 --pheno-type-col data_type --out OUTPUT_FILE.tsv
```

## Acknowledgements
Tuomo Kiiskinen for the original matching algorithm  
Aki Havulinna & Tuomo Kiiskinen, clinical expert groups & others at FIMM & THL for the FinnGen phenotype definitions
