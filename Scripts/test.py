import unittest

import pandas as pd, numpy as np
from typing import AbstractSet, List, Dict, Optional, NamedTuple, Any
import itertools
from tree import * 
from data_cleaning import *
from join import match_endpoints, process_matches, write_matches 
from constants import FG_REGEX_COL, FG_MATCHING_ICD, ICD_MAP_COL


class TestConstants(unittest.TestCase):
    def test_constants_not_empty(self):
        """
        Test that the imported constants are not empty
        """
        self.assertIsNotNone(FG_REGEX_COL)
        self.assertIsNotNone(FG_MATCHING_ICD)
        self.assertIsNotNone(ICD_MAP_COL)

if __name__ == '__main__':
    unittest.main()