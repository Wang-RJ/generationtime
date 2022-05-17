#!/usr/bin/env python3

"""
This script reads a dataframe and a bed file with introgressed sites 
coordinates, and adds these coordiantes to the dataframe.
"""

import numpy as np
import pandas as pd
import sys

in_df = sys.argv[1]     # input dataframe
in_mask = sys.argv[2]   # bed file
out_df = sys.argv[3]    # output dataframe

header_names = ['Chromosome', 'Position', 'AlleleRef', 'AlleleAlt', 'AlleleAnc', 'DataSource', 'AgeMedian_Jnt', 'Ancestral', 'Derived', 'RecRate', 'ID', 'rs', 'Ref', 'Alt', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS', 'sameAlt']

data = pd.read_csv(in_df, sep='\t', header=None, names=header_names)
mask = pd.read_csv(in_mask, sep='\t', header=None, names=['chr', 'pos'])


mask_col = data['Position'].isin(mask['pos'])
data['Masked'] = mask_col

data.to_csv(out_df+'.tsv', sep='\t', index=False)
