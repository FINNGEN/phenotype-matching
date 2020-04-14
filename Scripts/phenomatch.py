#! /usr/bin/env python3

import pandas as pd, numpy as np
import argparse, re
from typing import AbstractSet, List, Dict, Optional
import itertools
from tree import * 

def load_phecode_data(args) -> pd.DataFrame:
    raise NotImplementedError("Loading of phecode data not implemented!")

def load_fg_data(args) -> pd.DataFrame:
    raise NotImplementedError("Loading of FG data not implemented!")

def join_data(fg_data,phecode_data,args) -> pd.DataFrame:
    raise NotImplementedError("Joining of data not implemented!")

if __name__ == "__main__":
    parser=argparse.ArgumentParser("Map FG and UKBB/ICD10 codes to each other, either left or right join")
    parser.add_argument("--phecode-source",required=True,help="Phecode/ICD10 file")
    parser.add_argument("--fg-source",required=True,help="FinnGen file")
    parser.add_argument("--map-source",required=True,help="Phecode/ICD10 mapping file")
    parser.add_argument("--pheno-col-phe",required=True,help="Phenotype column in phecode/ICD10 file")
    parser.add_argument("--pheno-col-fg",required=True,help="Phenotype column in FG file")
    parser.add_argument("--pheno-col-map",required=True,help="PheCode column in Phecode/ICD10 mapping")
    parser.add_argument("--icd-col-map",required=True,help="ICD10 column in Phe/ICD mapping")
    parser.add_argument("--icd-col-fg",required=True,nargs="+",help="ICD10 column in FG file")
    parser.add_argument("--phenotype-type-col",required=True,help="Phetype type column in PheCode/ICD10 file")
    parser.add_argument("--include-col-fg",required=True,help="The column which lists the FG endpoints included in an endpoint")
    parser.add_argument("--out",required=True, help="Output filename")

    args=parser.parse_args()
