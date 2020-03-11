# Map phenotypes using ICD-10 as intermediary

This repository contains scripts to map FG endpoints to PheCodes.

The mapping is done by finding the target phenotype that best matches the source phenotype, by ICD10 similarity.

The similarity is calculated as follows:

FinnGen ICD10 codes are defined as a regular expression that matches some set of ICD-10 codes. As for Phecodes, the PheCode is enumerated for many ICD-10 codes.

The similarity is then calculated by counting, how many of a phecode's icd codes were matched by the FG regex. This is then divided by the total number of phecode icd codes.
