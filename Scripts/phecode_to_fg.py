#! /usr/bin/env python3

import pandas as pd, numpy as np
from fg_to_phecode import *
import argparse, re
from typing import AbstractSet, List, Dict
import itertools

def phenotype_data_filtering(df: pd.DataFrame) -> pd.DataFrame :
    #filter the input data so it only contains ice10s and phecodes
    df= df.loc[df["n_cases_both_sexes"]>100,:]     
    df= df[df["coding"]!="icd9"]  
    df= df[~df["pheno"].isin(['22601', '22617', '20024', '41230', '41210'])]   
    df= df[df["pop"]=="EUR"] 
    df = df[df["data_type"].isin(["icd_all","phecode"])]
    return df

def fg_matches(reg: str, lst: List[str]) -> List[str]:
    try:
        return [a for a in lst if bool(re.match(reg,a))]
    except:
        print(reg, lst)
        raise

def map_fg_for_phenotypes(args) -> pd.DataFrame:
    #load phecode file, already filtered down to icd10/phecodes
    pheno_data = pd.read_csv(args.phecode_source,sep="\t",dtype={args.pheno_col_phe:str})
    map_data = pd.read_csv(args.map_source,sep=",",dtype={args.pheno_col_map:str})
    pheno_data = pheno_data.fillna("")
    pheno_data = phenotype_data_filtering(pheno_data) #filter the data to only phecodes and icd10
    map_data = map_data.fillna("")
    map_data[args.icd_col_map] =map_data[args.icd_col_map].apply(lambda x: str(x).replace(".","")) 
    #divide dataset into phecodes and icd10s
    phecode_id = "phecode" #from datafile, exomes_full_run etc
    icd_id = "icd_all" #same as above
    phecode_data = pheno_data[pheno_data[args.phenotype_type_col] == phecode_id].copy()
    icd_data = pheno_data[pheno_data[args.phenotype_type_col] == icd_id].copy()
    #for icd10, copy icd10 ,column to mapping column
    icd_data["icd10_map_col"] = icd_data[args.pheno_col_phe].apply(lambda x:str(x).replace(".",""))
    #for phecodes, use phecode mapping to get matching list of phecodes
    #aggregate map data per phecode
    icd_codes_from_map=map_data.groupby(args.pheno_col_map).aggregate({args.icd_col_map:";".join}).reset_index() #colname from input data
    #join map data with icd_data
    phecode_data = phecode_data.merge(icd_codes_from_map[[args.pheno_col_map,args.icd_col_map]].rename(columns={args.icd_col_map:"icd10_map_col"}),
        how="left",
        left_on=args.pheno_col_phe,
        right_on=args.pheno_col_map,
        sort=False
        )
    phecode_data = phecode_data.drop(columns=[args.pheno_col_map])
    #concat those two DFs
    pheno_data = pd.concat([phecode_data,icd_data],sort=False).reset_index(drop=True)
    pheno_data = pheno_data.dropna(subset=["icd10_map_col"])
    #for each phenotype, find closest FG phenotype
    #read in fg data
    fg_data = pd.read_csv(args.fg_source,sep="\t")
    #constrain fg to be only pheno name, icd10
    fg_data = fg_data[[args.pheno_col_fg, args.icd_col_fg]]
    fg_data = fg_data.dropna()
    fg_data[args.icd_col_fg] = fg_data[args.icd_col_fg].apply(lambda x: str(x).replace(".",""))
    #for each phenotype, get the closest fg match.
    print("Get ICD codes for FG...")
    icd_codes = map_data[args.icd_col_map].unique()
    fg_data["matching_ICD"] = fg_data[args.icd_col_fg].apply(lambda x: ";".join(fg_matches(x,icd_codes)))
    fg_dict={}
    fg_regex_dict={}
    print("Create FG dict...")
    for t in fg_data.itertuples():
        fg_dict[getattr(t, args.pheno_col_fg)] = t.matching_ICD
        fg_regex_dict[getattr(t, args.pheno_col_fg)] = getattr(t,args.icd_col_fg)
    #first, get all fg matches for a single . Then, rank them according to similarity score.
    print("start to match PheCodes to FG...")
    output=[]
    am_rows = pheno_data.shape[0]
    i = 0
    for t in pheno_data.itertuples():
        record={}
        record["phenotype"] = getattr(t,args.pheno_col_phe)
        record["phenotype_icd10_matches"] = getattr(t,"icd10_map_col")
        icd_matches = set(getattr(t,"icd10_map_col").split(";"))
        fg_match_list = []
        fg_score_list = []
        for key,value in fg_dict.items():
            similarity_score=union_similarity(icd_matches,value.split(";"))
            if similarity_score > 0:
                fg_match_list.append(key)
                fg_score_list.append(similarity_score)
        #get best score
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
        output.append(record)
        
        if i % (am_rows//10) == 0:
            print("Progress: {:3.1f}%".format(i/am_rows*100))
        i+=1
    #output
    output_df = pd.DataFrame(output)
    output_df = pheno_data.merge(output_df, how="left", left_on=args.pheno_col_phe, right_on= "phenotype",sort=False)
    output_df=output_df.drop(columns=["phenotype","icd10_map_col"])
    return output_df

if __name__ == "__main__":
    parser=argparse.ArgumentParser("Map Phecodes/ICD10 to Finngen (best FG phenotype for each ICD10/PheCode) using ICD10 as intermediary")
    parser.add_argument("--phecode-source",required=True,help="Phecode/ICD10 file")
    parser.add_argument("--fg-source",required=True,help="FinnGen file")
    parser.add_argument("--map-source",required=True,help="Phecode/ICD10 mapping file")
    parser.add_argument("--pheno-col-phe",required=True,help="Phenotype column in phecode/ICD10 file")
    parser.add_argument("--pheno-col-fg",required=True,help="Phenotype column in FG file")
    parser.add_argument("--pheno-col-map",required=True,help="PheCode column in Phecode/ICD10 mapping")
    parser.add_argument("--icd-col-map",required=True,help="ICD10 column in Phe/ICD mapping")
    parser.add_argument("--icd-col-fg",required=True,help="ICD10 column in FG file")
    parser.add_argument("--phenotype-type-col",required=True,help="Phetype type column in PheCode/ICD10 file")
    parser.add_argument("--out",required=True, help="Output filename")
    args = parser.parse_args()
    output=map_fg_for_phenotypes(args)
    output.to_csv(args.out,sep="\t",index=False)

