#! /usr/bin/env python3

import pandas as pd, numpy as np #typing: ignore
import argparse, re
from typing import AbstractSet, List, Dict, Optional, NamedTuple, Any
import itertools
from tree import * 
from data_cleaning import *
from join import match_endpoints, process_matches, write_matches 
from constants import FG_REGEX_COL, FG_MATCHING_ICD, ICD_MAP_COL

if __name__ == "__main__":

    parser=argparse.ArgumentParser("Map FG and UKBB/ICD10 codes to each other, either left or right join")
    parser.add_argument("--main-table",required=True, choices=["phecode","finngen"],help="The direction of the join goes from auxiliary data to main table. So 'phecode' would map FG endpoints to the phecode data.")
    parser.add_argument("--out",required=True, help="Output filename")
    #Phecode data
    phecode_parser = parser.add_argument_group("phecode")
    phecode_parser.add_argument("--phecode-source",required=True,help="Phecode/ICD10 file with phenotypes that are to be matched with FG data")
    phecode_parser.add_argument("--pheno-pheno-col",required=True,help="Phenotype column in phecode/ICD10 file")
    phecode_parser.add_argument("--pheno-type-col",required=True,help="Phenotype type column in PheCode/ICD10 file")
    phecode_parser.add_argument("--pheno-sep",default='\t',help="Phecode file separator")
    #Mapping data
    mapping_parser = parser.add_argument_group("mapping")
    mapping_parser.add_argument("--map-source",required=True,help="Phecode/ICD10 mapping file")
    mapping_parser.add_argument("--map-pheno-col",required=True,help="PheCode column in Phecode/ICD10 mapping")
    mapping_parser.add_argument("--map-icd-col",required=True,help="ICD10 column in Phe/ICD mapping")
    mapping_parser.add_argument("--map-sep",default='\t',help="mapping file separator")
    #FinnGen data
    fg_parser = parser.add_argument_group("finngen")
    fg_parser.add_argument("--fg-source",required=True,help="FinnGen file")
    fg_parser.add_argument("--fg-pheno-col",required=True,help="Phenotype column in FG file")
    fg_parser.add_argument("--fg-icd-col",required=True,nargs="+",help="ICD10 columns in FG file")
    fg_parser.add_argument("--fg-inc-col",required=True,help="The column which lists the FG endpoints included in an endpoint")
    fg_parser.add_argument("--fg-sep",default='\t',help="FinnGen file separator")

    args=parser.parse_args()

    join_direction = args.main_table

    print("Load data...",end="\r")
    pheno_data_ = pd.read_csv(args.phecode_source, sep = args.pheno_sep,usecols=[args.pheno_pheno_col, args.pheno_type_col])
    map_data_ = pd.read_csv(args.map_source, sep = args.map_sep, dtype={args.map_pheno_col:str},usecols=[args.map_pheno_col, args.map_icd_col])
    fg_data_ = pd.read_csv(args.fg_source, sep = args.fg_sep,na_values="NA",usecols=args.fg_icd_col+[args.fg_pheno_col, args.fg_inc_col])
    print("Load data... Done")

    print("Prepare data for joining...")
    #clean and transform data
    pheno_data = pheno_data_.copy().fillna("")
    map_data = clean_map_data(map_data_.copy(), args.map_icd_col)
    icd_codes = get_icd_codes(map_data.copy(),args.map_icd_col)
    phecode_data = prepare_phecode_data(pheno_data, map_data, args.pheno_pheno_col, args.pheno_type_col,args.map_pheno_col, args.map_icd_col)
    fg_data = prepare_fg_data(fg_data_.copy(), args.fg_icd_col, args.fg_inc_col, args.fg_pheno_col)
    fg_endpoints = create_fg_endpoints(fg_data,icd_codes,args.fg_pheno_col)
    phecode_endpoints = create_phecode_endpoints(phecode_data,args.pheno_pheno_col)
    print("Prepare data for joining... Done")

    print("Join data...",end="\r")
    if args.main_table == "finngen":
        matches = match_endpoints(fg_endpoints,phecode_endpoints)
    else:
        matches = match_endpoints(phecode_endpoints,fg_endpoints)
    processed_matches = process_matches(matches)
    print("Join data... Done")
    
    print("Write output to {}".format(args.out))
    write_matches(processed_matches,args.out)