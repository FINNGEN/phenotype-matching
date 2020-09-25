#!/usr/bin/env python3

import pandas as pd, numpy as np
from fg_to_phecode import *
import argparse, re
from typing import AbstractSet, List, Dict, Optional
import itertools
from tree import *

def build_dependency_tree(fg_df: pd.DataFrame, pheno: str, pheno_colname: str, icd_colname: str, include_colname: str, rec: bool=False, nodeset: Optional[AbstractSet]=None) -> Tree :
    """Build a tree from the include dependency chains, removing any cycles if necessary.
    """
    try:
        row = fg_df[fg_df[pheno_colname] == pheno ].iloc[0]
    except:
        return
    node_data = row[icd_colname]
    if not rec:
        nodeset=set()
    nodeset.add(pheno)
    subtree = Tree(pheno, node_data)
    if pd.isna(row[include_colname]):
        return subtree
    for c in row[include_colname].split("|"):
        if c not in nodeset:
            subnode = build_dependency_tree(fg_df,c,pheno_colname,icd_colname,include_colname,True,nodeset)
            if subnode != None:
                subtree.add_child(subnode)
            else:
                print("Warning: phenotype {} has included phenotype {} that does not exist.".format(pheno,c))

    return subtree


def get_matches(reg: str, lst: List[str]) -> List[str]:
    """Match list of strings to regex, returning those strings that match the regex.
    """
    retlist= [a for a in lst if bool(re.match(reg,a))]
    if not retlist:
        retlist=[reg] #in case the reg did not match any strings, it itself is a valid non-matching string. 'NA' would match with other 'NA'
    return retlist

def fg_combine_regexes(x: List[str]) -> str:
    """Combine regex expressions (with OR, not AND) into one regex.
    """
    reg_lst = []
    [reg_lst.append(tmp) for tmp in x if (((tmp not in reg_lst) and (tmp != "") )and (tmp != "$!$"))]
    return "|".join(reg_lst)

def solve_includes(fg_df: pd.DataFrame, pheno: str, pheno_colname: str, icd_colname: str, include_colname: str) -> str:
    """Solve regex column when multiple endpoints are included in endpoint. Handles cyclical cases. Does not handle missing phenotype names.
    """
    icds = "|".join(a for a in get_tree_nodes( build_dependency_tree(fg_df,pheno,pheno_colname,icd_colname,include_colname)).values() if a != "")
    return icds
