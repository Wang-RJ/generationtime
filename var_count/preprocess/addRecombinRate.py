#!/usr/bin/env python

"""
This script read a recombination map file and add the a column for 
rates to our dataframe.
"""

import numpy as np
import pandas as pd
import sys

recmap = sys.argv[1]    # recombination map file
in_df = sys.argv[2]     # dataframe

genmap = pd.read_csv(recmap, sep='\t', index_col=False)

data = pd.read_csv(in_df, sep='\t', index_col=False)

frompos = []
prev = 1
for i, pos in enumerate(genmap['Position(bp)']):
    if i == 0:
        frompos.append(prev)
        prev = pos
    else:
        frompos.append(prev+1)
        prev = pos

genmap['FromPos'] = frompos

Rec_Rate = []
idx = 0
rcount = 0
for i, row in genmap.iterrows():
    for i, pos in enumerate(data['Position'][idx:]):
        if pos >= row[4] and pos <= row[1]:
            rcount += 1
            Rec_Rate.append(row[2])
        else:
            idx = rcount
            break

rate = pd.DataFrame(data=np.array(Rec_Rate), columns=['Rate'])
df = pd.concat([data, rate], ignore_index=True, axis=1)

df.to_csv('chr1-rec.txt', sep='\t', index=False, header=True)


