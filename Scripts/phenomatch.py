#! /usr/bin/env python3

import pandas as pd, numpy as np #typing: ignore
import argparse, re
from typing import AbstractSet, List, Dict, Optional
import itertools
from tree import * 
from phecode_to_fg import phenotype_data_filtering, get_matches, fg_combine_regexes, solve_includes #typing: ignore
from fg_to_phecode import union_similarity

#constants
FG_REGEX_COL="fg_regex_column"
FG_MATCHING_ICD="matching_icd"
ICD_MAP_COL="icd_map_col"

def clean_map_data(map_data: pd.DataFrame, map_icd_col: str) -> pd.DataFrame:
    map_data = map_data.fillna("")
    map_data[map_icd_col] =map_data[map_icd_col].apply(lambda x: str(x).replace(".",""))
    return map_data

def get_icd_codes(map_data,icd_column) -> List[str]:
    """Get List if ICD codes from ICD code map
    """
    return map_data[icd_column].unique()

def prepare_phecode_data(pheno_data: pd.DataFrame, map_data: pd.DataFrame, pheno_pheno_col: str, pheno_type_col: str, map_pheno_col: str, map_icd_col: str) -> pd.DataFrame:

    # icd10 and phecode phenotype codes
    phecode_id = "phecode"
    icd_id = "icd_all"

    #separate icd10 and phecode entries in pheno_data
    phecode_data = pheno_data[pheno_data[pheno_type_col] == phecode_id].copy()
    icd_data = pheno_data[pheno_data[pheno_type_col] == icd_id].copy()

    icd_codes_from_map=map_data.groupby(map_pheno_col).aggregate({map_icd_col:";".join}).reset_index()
    phecode_data = phecode_data.merge(icd_codes_from_map[[map_pheno_col,map_icd_col]].rename(columns={map_icd_col:ICD_MAP_COL}),
        how="left",
        left_on=pheno_pheno_col,
        right_on=map_pheno_col,
        sort=False
        )
    phecode_data = phecode_data.drop(columns=[map_pheno_col])
    icd_codes = get_icd_codes(map_data,map_icd_col)
    icd_data[ICD_MAP_COL] = icd_data[pheno_pheno_col].apply(lambda x:str(x).replace(".",""))
    icd_data[ICD_MAP_COL] = icd_data[ICD_MAP_COL].apply(lambda x:";".join( get_matches(x,icd_codes) ) )
    
    pheno_data = pd.concat([phecode_data,icd_data],sort=False).reset_index(drop=True)
    pheno_data = pheno_data.fillna("")

    return pheno_data

def prepare_fg_data(fg_data: pd.DataFrame, fg_icd_col: List[str], fg_inc_col: str, fg_pheno_col: str) -> pd.DataFrame:
    
    fg_cols = [fg_pheno_col, fg_inc_col]+fg_icd_col
    fg_data=fg_data[fg_cols]
    fg_data=fg_data.dropna(how="all")
    
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

def join_data(phecode_data,fg_data,icd_codes,join_direction,fg_pheno_col, pheno_pheno_col) -> pd.DataFrame:
    #if join direction is right
    if join_direction == "left":
        ##add icd codes for fg
        fg_data[FG_MATCHING_ICD] = fg_data[FG_REGEX_COL].apply(lambda x:";".join(get_matches(x,icd_codes)) if pd.notna(x) and x!="" else "NAN")
        ##create fg dicts
        fg_dict={}
        fg_regex_dict={}
        for t in fg_data.itertuples():
            fg_dict[getattr(t, fg_pheno_col)] = getattr(t, FG_MATCHING_ICD)
            fg_regex_dict[getattr(t,fg_pheno_col)] = getattr(t,FG_REGEX_COL)
        ##calculate best matching FG pheno for phenotype
        output_records: List[Dict[str,Any]]=[]
        for t in phecode_data.itertuples():
            record={}
            record["phenotype"] = getattr(t,pheno_pheno_col)
            record["phenotype_icd10_matches"] = getattr(t,ICD_MAP_COL)
            icd_matches = set(getattr(t,ICD_MAP_COL).split(";"))
            fg_match_list = []
            fg_score_list = []
            for phenoname,icds in fg_dict.items():
                similarity_score=union_similarity(icd_matches,icds.split(";"))
                if similarity_score > 0:
                    fg_match_list.append(phenoname)
                    fg_score_list.append(similarity_score)
            if len(fg_score_list)>0:
                best_score_idx = np.argmax(fg_score_list)
                best_score = fg_score_list[best_score_idx]
                best_pheno = fg_match_list[best_score_idx]
                best_matches = fg_dict[best_pheno]
                best_reg = fg_regex_dict[best_pheno]
                listed_phenos=";".join( "{}|{:.3g}".format(a,b) for (a,b) in itertools.zip_longest(fg_match_list,fg_score_list ) )
            else:
                best_pheno="NA"
                best_score="NA"
                listed_phenos="NA"
                best_matches = "NA"
                best_reg = "NA"
            record["best_fg_phenotype"] = best_pheno
            record["fg_icd10"] = best_matches
            record["best_fg_score"] = best_score
            record["all_phenos"] = listed_phenos
            record["fg_regex"] = best_reg
            output_records.append(record)
        output_df=pd.DataFrame(output_records)
        output_df = phecode_data.merge(output_df,how="left",left_on=pheno_pheno_col,right_on="phenotype",sort=False)
        output_df = output_df.drop(columns=["phenotype",ICD_MAP_COL])
        return output_df

    #if join directions is left
    elif join_direction =="right":
        raise NotImplementedError("Joining of data in left direction not yet implemented!")
    else:
        raise NotImplementedError("Joining in other direction than right or left not implemented!!")
if __name__ == "__main__":
    parser=argparse.ArgumentParser("Map FG and UKBB/ICD10 codes to each other, either left or right join")
    parser.add_argument("--main-table",required=True, choices=["phecode","finngen"],help="The direction of the join goes from auxiliary data to main table. So 'phecode' would map FG endpoints to the phecode data.")
    parser.add_argument("--out",required=True, help="Output filename")

    phecode_parser = parser.add_argument_group("phecode")
    phecode_parser.add_argument("--phecode-source",required=True,help="Phecode/ICD10 file")
    phecode_parser.add_argument("--pheno-pheno-col",required=True,help="Phenotype column in phecode/ICD10 file")
    phecode_parser.add_argument("--pheno-type-col",required=True,help="Phenotype type column in PheCode/ICD10 file")
    phecode_parser.add_argument("--map-source",required=True,help="Phecode/ICD10 mapping file")
    phecode_parser.add_argument("--map-pheno-col",required=True,help="PheCode column in Phecode/ICD10 mapping")
    phecode_parser.add_argument("--map-icd-col",required=True,help="ICD10 column in Phe/ICD mapping")

    fg_parser = parser.add_argument_group("finngen")
    fg_parser.add_argument("--fg-source",required=True,help="FinnGen file")
    fg_parser.add_argument("--fg-pheno-col",required=True,help="Phenotype column in FG file")
    fg_parser.add_argument("--fg-icd-col",required=True,nargs="+",help="ICD10 column in FG file")
    fg_parser.add_argument("--fg-inc-col",required=True,help="The column which lists the FG endpoints included in an endpoint")

    args=parser.parse_args()

    join_direction = "left" if args.main_table == "phecode"  else "right"

    print("Load data...",end="\r")
    pheno_data_ = pd.read_csv(args.phecode_source, sep = '\t')
    map_data_ = pd.read_csv(args.map_source, sep = ',', dtype={args.map_pheno_col:str})
    fg_data_ = pd.read_csv(args.fg_source, sep = '\t')
    print("Load data... Done")

    #clean phecode data
    pheno_data = pheno_data_.copy().fillna("")
    pheno_data=phenotype_data_filtering(pheno_data)
    map_data = clean_map_data(map_data_.copy(), args.map_icd_col)

    print("Prepare data for joining...")
    icd_codes = get_icd_codes(map_data.copy(),args.map_icd_col)
    phecode_data = prepare_phecode_data(pheno_data, map_data, args.pheno_pheno_col, args.pheno_type_col,args.map_pheno_col, args.map_icd_col)
    fg_data = prepare_fg_data(fg_data_.copy(), args.fg_icd_col, args.fg_inc_col, args.fg_pheno_col)
    print("Prepare data for joining... Done")

    print("Join data...",end="\r")
    joined_data = join_data(phecode_data,fg_data,icd_codes,join_direction,args.fg_pheno_col,args.pheno_pheno_col)
    print("Join data... Done")
    
    print("Write output to {}".format(args.out))
    joined_data.to_csv(args.out, sep='\t',index=False)