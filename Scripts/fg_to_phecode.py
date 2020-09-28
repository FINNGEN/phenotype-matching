#!/usr/bin/env python3
import argparse
import pandas as pd, numpy as np
from expand_codes import *
import re
from typing import List, Dict, AbstractSet, Any

def union_similarity(match_set: AbstractSet[str], match_list: List[str]):
    intersection = match_set.intersection(match_list)
    union = match_set.union(match_list)
    return len(intersection)/len(union)
