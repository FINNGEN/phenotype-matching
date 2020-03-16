#!/usr/bin/python3
import argparse
import pandas as pd, numpy as np
from expand_codes import *
import re
from typing import List, Dict, AbstractSet, Any

def union_similarity(fg_set: AbstractSet[str], pheno_list: List[str]):
    intersection = fg_set.intersection(pheno_list)
    union = fg_set.union(pheno_list)
    return len(intersection)/len(union)


def main(args):
    """Find matching PheCodes for FinnGen phenotypes
    """
    #load files into dataframes etc
    source_1 = pd.read_csv(args.source1, sep="\t",dtype=str)
    source_2 = pd.read_csv(args.source2, sep="\t",dtype=str)
    #restrict to relevant columns
    df_1 = source_1[[args.pheno_col_1,args.icd_col_1]].copy()
    df_1=df_1.dropna()
    df_2 = source_2[[args.pheno_col_2,args.icd_col_2]].copy()
    df_2=df_2.dropna()
    df2_orig=df_2.copy()
    #FG has regex, the other has separate codes
    #aggregate second based on phecodes
    df_2 = df_2.groupby([args.pheno_col_2]).aggregate({args.icd_col_2:";".join}).reset_index()
    phecode_dict={}
    df_2_records = df_2.to_dict("records")
    for entry in df_2_records:
        phecode_dict[entry[args.pheno_col_2]] = entry[args.icd_col_2].replace(".","").split(";")
    #for each fg phenotype, calculate scores for each phecode phenotype
    fg_dict={}
    df_1_records=df_1.to_dict("records")
    for entry in df_1_records:
        fg_dict[entry[args.pheno_col_1]] = entry[args.icd_col_1]
    out=[]
    all_icd10 = df2_orig[args.icd_col_2].apply(lambda x: x.replace(".","")).unique()
    phecodes={}
    phecodes["key"]=[]
    phecodes["value"]=[]
    for (key,value) in phecode_dict.items():
        phecodes["key"].append(key)
        phecodes["value"].append(value)
    cnt=0
    for (pheno,reg) in fg_dict.items():
        #calculate normalizing constant for similarity v2
        fg_set = set((a for a in all_icd10 if re.match(reg,a)))# amount of ICD10 matching to the regular expression in the ICD10-PheCode mapping
        #that is the closest approximation to enumerating FG matches in the ICD10 that I can do without enumerating them
        entry={}
        entry["phenotype"]=pheno
        phecode_list=[]
        similarity_score = []
        for value in phecodes["value"]:
            #union similarity
            sim_2 = union_similarity(fg_set,value)
            similarity_score.append(sim_2)
        idx = np.argmax(similarity_score)

        entry["fg_icd10"]=reg
        if similarity_score[idx]>0:
            entry["phecode"]=phecodes["key"][idx]
            entry["phecode_icd10"]="|".join(phecodes["value"][idx])
            entry["similarity"]=similarity_score[idx]
        else:
            entry["phecode"]=np.nan
            entry["phecode_icd10"]=np.nan
            entry["similarity"]=np.nan
        out.append(entry)
        if cnt%100 == 0:
            print("{:5.1f}%".format(cnt/len(fg_dict.keys())*100  ) )
        cnt+=1
    #construct output file
    out_df=pd.DataFrame(out)
    out_df.to_csv(args.out,sep="\t",index=False)


if __name__ == "__main__":
    parser=argparse.ArgumentParser("Map phenotypes using ICD-10 as intermediary")
    parser.add_argument("--source1",required=True,help="Phenotype source 1")
    parser.add_argument("--source2",required=True,help="Phenotype source 2")
    parser.add_argument("--pheno-col-1",required=True,help="Phenotype column in source 1")
    parser.add_argument("--pheno-col-2",required=True,help="Phenotype column in source 2")
    parser.add_argument("--icd-col-1",required=True,help="ICD10 column in source 1")
    parser.add_argument("--icd-col-2",required=True,help="ICD10 column in source 2")
    parser.add_argument("--out",required=True, help="Output filename")
    args = parser.parse_args()
    main(args)