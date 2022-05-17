#!/usr/bin/env python

"""
This script join introgression annotations across all individuals and
all populations found in the Steinrucken et al, 2018 paper, and outputs
a bed file with combined region annotations. 

The input is a directory containing all files. 
"""

import os
import sys
import pandas as pd

indir = sys.argv[1]

files = os.listdir(indir)

biglist = []
for file in files:
    if os.path.isfile(os.path.join(indir, file)):
        #with open(os.path.join(indir, file), 'r') as f:
        df = pd.read_csv(os.path.join(indir, file), sep='\t', names=['chr', 'start', 'end'])
        df['length'] = df['end'] - df['start']
        sum_diff = sum(df['length'])
        print(f"{file[:7]}\t{sum_diff}")


