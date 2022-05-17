#!/usr/bin/env python

"""
Add CI columns from atlas files to opur dataframe.
"""

import numpy as np
import pandas as pd
import sys

in_df = sys.argv[1]     # input dataframe
in_atlas = sys.argv[2]  # atlas file
out_df = sys.argv[3]    # output dataframe

data = pd.read_csv(in_df, sep='\t', index_col=False, skipinitialspace=True)
atlas = pd.read_csv(in_atlas, comment='#', skipinitialspace=True, usecols=['Position', 'AgeCI95Lower_Jnt', 'AgeCI95Upper_Jnt'], index_col=False)

data2 = pd.merge(data, atlas, on='Position', how='inner')

data3 = data2.drop_duplicates(subset='Position', keep='first', inplace=False, ignore_index=True)

data3.to_csv(out_df+'.tsv', sep='\t', index=False)
