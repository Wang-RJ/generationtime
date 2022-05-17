#!/usr/bin/env python

"""
This script reads a dataframe and a bed file with introgression
coordinates, and outputs a new bed file containing the introgressed 
sites to be removed later.
"""

import numpy as np
import pandas as pd
import sys

chrom = int(sys.argv[1])    # chromosome number
in_df = sys.argv[2]     # dataframe
in_bed = sys.argv[3]    # coordinate file

data = pd.read_csv(in_df, sep='\t', usecols=[0, 1])
mask = pd.read_csv(in_bed, sep='\t')

idx1 = data['Chromosome'] == chrom
idx2 = mask['chrom'] == chrom

pointer = 0
for i, row in mask[idx2].iterrows():
  for p, pos in data[idx1].iterrows():
    if np.logical_and(pos['Position'] >= row['start'], pos['Position'] <= row['end']):
      print('%s\t%d' % (chrom, pos['Position']))
    elif pos['Position'] > row['end']:
      break
    else:
      continue


