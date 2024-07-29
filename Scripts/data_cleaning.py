#! /usr/bin/env python3

import pandas as pd, numpy as np #typing: ignore
import argparse, re
from typing import AbstractSet, List, Dict, Optional, NamedTuple
import itertools
from tree import *
from join import Endpoint, EndpointMatch, Result
from constants import FG_REGEX_COL, FG_MATCHING_ICD, ICD_MAP_COL

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

def solve_includes(fg_df: pd.DataFrame, pheno: str, pheno_colname: str, icd_colname: str, include_colname: str) -> str:
    """Solve regex column when multiple endpoints are included in endpoint. Handles cyclical cases. Does not handle missing phenotype names.
    """
    icds = "|".join(a for a in get_tree_nodes( build_dependency_tree(fg_df,pheno,pheno_colname,icd_colname,include_colname)).values() if a != "")
    return icds

def get_icd_codes(map_data: pd.DataFrame, icd_column: str) -> List[str]:
    """Get List if ICD codes from ICD code map
    """
    return map_data[icd_column].unique()

def get_matches(reg: str, lst: List[str]) -> List[str]:
    """Match list of strings to regex, returning those strings that match the regex.
    """
    retlist= [a for a in lst if bool(re.match(reg,a))]
    if not retlist:
        return []
    return retlist

def create_fg_endpoints(fg_df: pd.DataFrame, icd_codes: List[str],fg_pheno_col)-> List[Endpoint]:
    """Create the finngen endpoint list
    """
    out=[]
    for t in fg_df.itertuples():
        out.append(
            Endpoint(
                getattr(t,fg_pheno_col),
                set(get_matches(getattr(t,FG_REGEX_COL),icd_codes ) ),
                getattr(t,FG_REGEX_COL)
            )
        )
    return out

def create_phecode_endpoints(phecode_df: pd.DataFrame, pheno_pheno_col: str) -> List[Endpoint]:
    """Create phecode endpoint list
    """
    out=[]
    for t in phecode_df.itertuples():
        out.append(
            Endpoint(
                getattr(t,pheno_pheno_col),
                set(getattr(t,ICD_MAP_COL).split(";")),
                ""
            )
        )
    return out

def clean_map_data(map_data: pd.DataFrame, map_icd_col: str) -> pd.DataFrame:
    """Clean up map data
    """
    map_data = map_data.fillna("")
    map_data[map_icd_col] =map_data[map_icd_col].apply(lambda x: str(x).replace(".",""))
    return map_data

def create_phecode_data(pheno_data: pd.DataFrame, map_data: pd.DataFrame, pheno_pheno_col: str, pheno_type_col: str, map_pheno_col: str, map_icd_col: str)-> pd.DataFrame:
    phecode_id = "phecode"
    icd_id = "icd10"

    #separate phecodes and icd10-codes

    #aggregate icd10 codes for phecodes from the map data. These will contain only valid icd10 codes, since they are the same as map data icd10 codes

    #For icd10 codes: get the matches by matching to the icd10 code name (no dots)

    pass

def prepare_phecode_data(pheno_data: pd.DataFrame, map_data: pd.DataFrame, pheno_pheno_col: str, pheno_type_col: str, map_pheno_col: str, map_icd_col: str) -> pd.DataFrame:
    """Data preprocessing for phecode data
    """
    # icd10 and phecode phenotype codes
    phecode_id = "phecode"
    icd_id = "icd10"

    #separate phecodes and icd10-codes
    phecode_data = pheno_data[pheno_data[pheno_type_col] == phecode_id].copy()
    icd_data = pheno_data[pheno_data[pheno_type_col] == icd_id].copy()

    #aggregate icd10 codes for phecodes from the map data. These will contain only valid icd10 codes, since they are the same as map data icd10 codes
    icd_codes_from_map=map_data.groupby(map_pheno_col).aggregate({map_icd_col:";".join}).reset_index()
    #join to phecode data
    phecode_data = phecode_data.merge(icd_codes_from_map[[map_pheno_col,map_icd_col]].rename(columns={map_icd_col:ICD_MAP_COL}),
        how="left",
        left_on=pheno_pheno_col,
        right_on=map_pheno_col,
        sort=False
        )
    phecode_data = phecode_data.drop(columns=[map_pheno_col])

    #For icd10 codes: get all matching ICD10 codes by matching the icd10 code to the list of icd codes in map file
    if not icd_data.empty:
        icd_codes = get_icd_codes(map_data,map_icd_col)
        icd_data[ICD_MAP_COL] = icd_data[pheno_pheno_col].apply(lambda x:str(x).replace(".",""))
        icd_data[ICD_MAP_COL] = icd_data[ICD_MAP_COL].apply(lambda x:";".join( list(set(get_matches(x,icd_codes))) ) )
        
        pheno_data = pd.concat([phecode_data,icd_data],sort=False).reset_index(drop=True)
    else:
        pheno_data = phecode_data.reset_index(drop=True)
    
    pheno_data = pheno_data.fillna("")

    return pheno_data

def fg_combine_regexes(x: List[str]) -> str:
    """Combine regex expressions (with OR, not AND) into one regex.
    """
    reg_lst = []
    [reg_lst.append(tmp) for tmp in x if (((tmp not in reg_lst) and (tmp != "") )and (tmp != "$!$"))]
    return "|".join(reg_lst)

def prepare_fg_data(fg_data: pd.DataFrame, fg_icd_col: List[str], fg_inc_col: str, fg_pheno_col: str) -> pd.DataFrame:
    """Data preprocessing for FinnGen data
    """
    fg_data=fg_data.dropna(subset = fg_icd_col +[fg_inc_col],how="all")
    
    #combine multiple regex columns into one
    fg_data[fg_icd_col] = fg_data[fg_icd_col].applymap(lambda x: str(x).replace(".","") if pd.notna(x) else "")
    fg_data["fg_icd_regex"] = fg_data[fg_icd_col].apply(fg_combine_regexes,axis=1)

    #add included phenotypes' regexes to phenotypes regexes

    fg_data[FG_REGEX_COL]=np.nan

    #index the FG endpoints with and without included columns. Those without will have their regexes unchanged,
    # while those that have them will have their regexes augmented by them.
    no_includes = pd.isna(fg_data[fg_inc_col])
    includes = ~pd.isna(fg_data[fg_inc_col])
    fg_data.loc[no_includes,FG_REGEX_COL] = fg_data.loc[no_includes,"fg_icd_regex"]
    for t in fg_data.loc[includes,:].itertuples():
        fg_data.loc[getattr(t,"Index"),FG_REGEX_COL] = solve_includes(fg_data,getattr(t,fg_pheno_col),fg_pheno_col,"fg_icd_regex",fg_inc_col)

    return fg_data