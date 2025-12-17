#! /usr/bin/env python3

import pandas as pd, numpy as np
import re
from typing import AbstractSet, List, Optional
from tree import *
from join import Endpoint
from constants import FG_REGEX_COL, ICD_MAP_COL

from progress import Progress_bar

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
    icd_list = list(set([a for a in get_tree_nodes( build_dependency_tree(fg_df,pheno,pheno_colname,icd_colname,include_colname)).values() if a != ""]))
    icd_list.sort()
    icds = "|".join(icd_list)
    return icds

def get_icd_codes(map_data: pd.DataFrame, icd_column: str) -> List[str]:
    """Get List if ICD codes from ICD code map
    """
    return map_data[icd_column].unique()

def get_icd_codes_from_file(icd_code_file: str) -> List[str]:
    """Get List if ICD codes from ICD code file
    """
    icd_codes = []
    with open(icd_code_file) as f:
        for line in f:
            icd_codes.append(line.strip())
    return icd_codes

def get_matches(reg: str, lst: List[str]) -> List[str]:
    """Match list of strings to regex, returning those strings that match the regex.
    """
    retlist= [a for a in lst if bool(re.match(reg,a))]
    if not retlist:
        return []
    return retlist

def format_regex_from_icd_codes(icd_codes: AbstractSet[str]) -> str:
    """Format a regex from a list of ICD codes
    """
    cl = sorted(list(icd_codes), key = len)
    matched = []
    for c_a in cl:
        for c_b in cl:
            if c_b in matched:
                continue
            if c_a == c_b:
                continue
            if re.match(c_a, c_b):
                matched.append(c_b)
    regex = [c for c in cl if c not in matched]
    return "|".join(regex)

def format_regex_from_icd_string(icd_string: str) -> str:
    """Format a regex from a string of ICD codes
    """
    return format_regex_from_icd_codes(set(tokenize_icd_string(icd_string)))

def tokenize_icd_string(icd_string: str) -> List[str]:
    """Tokenize an ICD string.
    Very crude way to split most of the regexes in FG definitions.
    """
    tokens = []
    token = ""
    bracket_start_seen = False
    for c in icd_string:
        if c == "[":
            bracket_start_seen = True
        if c == "]":
            bracket_start_seen = False
        if c == "|" and not bracket_start_seen:
            tokens.append(token)
            token = ""
        else:
            token += c
    tokens.append(token)
    return tokens

def create_fg_endpoints(fg_df: pd.DataFrame, icd_codes: List[str], fg_pheno_col: str, exclude_cols: List[str])-> List[Endpoint]:
    """Create the finngen endpoint list
    """
    out=[]
    with Progress_bar as p:
        for t in p.track(fg_df.itertuples(), total=len(fg_df), description="Creating FinnGen endpoints..."):
            #get included icd codes
            incl_icd_codes = set(get_matches(getattr(t, FG_REGEX_COL), icd_codes))
            if exclude_cols:
                #get excluded icd codes
                excl_icd_codes = set()
                for c in exclude_cols:
                    if pd.notna(getattr(t,c)):
                        excl_icd_codes.update(set(get_matches(getattr(t, c), icd_codes)))
                #remove excluded icd codes from included icd codes
                if excl_icd_codes:
                    incl_icd_codes = incl_icd_codes - excl_icd_codes
                    if not incl_icd_codes:
                        print(f"Warning: FinnGen endpoint {getattr(t, fg_pheno_col)} has all ICD-10 codes excluded.")
                        continue
            out.append(
                Endpoint(
                    getattr(t, fg_pheno_col),
                    incl_icd_codes,
                    getattr(t, FG_REGEX_COL),
                    True if excl_icd_codes else False
                )
            )
    return out

def create_phecode_endpoints(phecode_df: pd.DataFrame, icd_codes: List[str], pheno_pheno_col: str) -> List[Endpoint]:
    """Create phecode endpoint list
    """
    out=[]
    with Progress_bar as p:
        for t in p.track(phecode_df.itertuples(), total=len(phecode_df), description="Creating phecode endpoints..."):
            icd_code_regex = format_regex_from_icd_codes(set(getattr(t,ICD_MAP_COL).split(";")))
            incl_icd_codes = set(get_matches(icd_code_regex, icd_codes))
            out.append(
                Endpoint(
                    getattr(t,pheno_pheno_col),
                    incl_icd_codes,
                    icd_code_regex,
                    False
                )
            )
    return out

def clean_map_data(map_data: pd.DataFrame, map_icd_col: str) -> pd.DataFrame:
    """Clean up map data
    """
    map_data = map_data.fillna("")
    map_data[map_icd_col] = map_data[map_icd_col].apply(lambda x: str(x).replace(".",""))
    map_data = map_data.drop_duplicates()
    return map_data

def create_phecode_data(pheno_data: pd.DataFrame, map_data: pd.DataFrame, pheno_pheno_col: str, pheno_type_col: str, map_pheno_col: str, map_icd_col: str)-> pd.DataFrame:
    phecode_id = "phecode"
    icd_id = "icd10"

    #separate phecodes and icd10-codes

    #aggregate icd10 codes for phecodes from the map data. These will contain only valid icd10 codes, since they are the same as map data icd10 codes

    #For icd10 codes: get the matches by matching to the icd10 code name (no dots)

    pass

def prepare_phecode_data(pheno_data: pd.DataFrame, map_data: pd.DataFrame, pheno_pheno_col: str, pheno_type_col: str, map_pheno_col: str, map_icd_col: str, filter: int) -> pd.DataFrame:
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
    if map_pheno_col != pheno_pheno_col:
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

    # Filter out phecodes mapping to code counts greater than filter threshold
    if filter > 0:
        pheno_data["code_count"] = pheno_data[ICD_MAP_COL].apply(lambda x: len(x.split(";")) if x else 0)
        N = len(pheno_data)
        pheno_data = pheno_data[pheno_data.code_count < filter]
        if len(pheno_data) < N:
            print(f"Removed {N-len(pheno_data)} phecodes mapping to more than {filter} ICD-10 codes.")

    return pheno_data

def fg_combine_regexes(x: List[str]) -> str:
    """Combine regex expressions (with OR, not AND) into one regex.
    """
    reg_lst = []
    [reg_lst.append(tmp) for tmp in x if tmp not in reg_lst and tmp != "" and tmp != "$!$" and tmp != "ANY"]
    return "|".join(reg_lst)

def prepare_fg_data(fg_data: pd.DataFrame, fg_icd_col: List[str], fg_inc_col: str, fg_pheno_col: str, fg_cond_col: List[str]) -> pd.DataFrame:
    """Data preprocessing for FinnGen data
    """
    N = len(fg_data)
    fg_data = fg_data.dropna(subset = fg_icd_col + [fg_inc_col], how = "all")
    if len(fg_data) < N:
        print(f"Removed {N-len(fg_data)} phenotypes with missing all ICD codes and included phenotypes.")

    N = len(fg_data)
    if fg_cond_col:
        fg_data = fg_data.loc[~fg_data.index.isin(fg_data.dropna(subset = fg_cond_col, how = "all").index)]
        if len(fg_data) < N:
            print(f"Removed {N-len(fg_data)} phenotypes with conditions.")
        N = len(fg_data)

    #remove phenotypes with '%' (mode) in any code column
    for c in fg_icd_col:
        fg_data = fg_data.loc[~fg_data[c].str.contains("%",na=False)]
    if len(fg_data) < N:
        print(f"Removed {N-len(fg_data)} phenotypes with '%' (mode) in ICD code.")

    #remove dots from ICD codes
    fg_data[fg_icd_col] = fg_data[fg_icd_col].applymap(lambda x: str(x).replace(".","") if pd.notna(x) else "")

    #combine multiple regex columns into one
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

    #format regexes
    fg_data[FG_REGEX_COL] = fg_data[FG_REGEX_COL].apply(format_regex_from_icd_string)

    return fg_data