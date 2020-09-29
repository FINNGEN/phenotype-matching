from typing import AbstractSet, List, Dict, Optional, NamedTuple, Any
import pandas as pd
import itertools

class Endpoint(NamedTuple):
    """Endpoint with a name and matching ICD10 codes
    """
    name: str
    matches: AbstractSet[str]
    regex: str

class EndpointMatch(NamedTuple):
    """Match between endpoints
    """
    endpoint_1: Endpoint
    endpoint_2: Endpoint
    score: float

class Result(NamedTuple):
    endpoint_1: str
    endpoint_2: str
    score: float
    matches_1: str
    matches_2: str
    regex_1: str
    regex_2: str
    other_hits: str

def union_similarity(match_set: AbstractSet[str], match_list: AbstractSet[str]):
    intersection = match_set.intersection(match_list)
    union = match_set.union(match_list)
    return len(intersection)/len(union)

def match_endpoints(endpoints_1: List[Endpoint], endpoints_2: List[Endpoint]) -> List[EndpointMatch]:
    """Join endpoints by calculating union similarity and output all matches over those endpoints 
    """
    nomatch = Endpoint("NO MATCH",{},"")
    matches = []
    for end_1 in endpoints_1:
        if end_1.matches:
            for end_2 in endpoints_2:
                if end_2.matches:
                    score = union_similarity(end_1.matches, end_2.matches)
                    if score>0.0:
                        matches.append(EndpointMatch(
                            end_1,
                            end_2,
                            score
                        ))
        else:
            matches.append(EndpointMatch(
                        end_1,
                        nomatch,
                        0.0
                    ))
    return matches

def process_matches(matches: List[EndpointMatch]) -> List[Result] :
    """Aggregate the matches on endpoint_1, resulting in a list of best matches for endpoint 1
    """
    out=[]
    #sort data
    sort_f = lambda match: match.endpoint_1.name
    sorted_matches = sorted(matches, key = sort_f)
    #group by endpoint 1
    for key, group in itertools.groupby(sorted_matches, key = sort_f ):
        #sort by score
        endpoint_matches = sorted(list(group), key = lambda match: match.score,reverse=True)
        #best match
        best_match = endpoint_matches[0]
        #others
        other_matches = ";".join([f"{a.endpoint_2.name}|{a.score:.3g}" for a in endpoint_matches[1:] ])

        out.append(Result(
            best_match.endpoint_1.name,
            best_match.endpoint_2.name,
            best_match.score,
            ";".join(best_match.endpoint_1.matches),
            ";".join(best_match.endpoint_2.matches),
            best_match.endpoint_1.regex,
            best_match.endpoint_2.regex,
            other_matches
        ))
    return out

def write_matches(matches: List[Result],out: str):
    """Write matches to tsv
    Currently using pandas dataframe flattening as the parser
    """
    data = pd.DataFrame(matches)
    data.to_csv(out,sep="\t",index=False)