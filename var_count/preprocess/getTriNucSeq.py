#!/usr/bin/env python3

"""
This script read the atlas and ref chromosome files and add columns
for tri-nucleotide mutation spectrum to its output.
"""

import pandas as pd
from pyfaidx import Fasta
import sys

fn = sys.argv[1]    # atlas file
ref = sys.argv[2]   # refseq file
out = sys.argv[3]   # name for output file

df = pd.read_csv(fn, usecols=[1,2,3,4,5,6,23], comment='#', skipinitialspace=True, index_col=False)
df = df[(df['AlleleAnc'] != '.')]
#df.drop_duplicates(subset=['Position'], keep='last', inplace=True)

chrom = Fasta(ref)
name = [i for i in chrom.keys()]


Ancestral, Derived = [], []
for i, row in df.iterrows():
    pos = row[1]-1
    prev_nuc = chrom[name[0]][pos-1].seq
    next_nuc = chrom[name[0]][pos+1].seq
    #nuc = chrom[name[0]][pos].seq
    #Ancestral.append(prev_nuc + nuc + next_nuc)
    Ancestral.append(prev_nuc.upper() + row[4] + next_nuc.upper())
    if row[4] == row[2]:
        Derived.append(prev_nuc.upper() + row[3] + next_nuc.upper())
#    elif row[4] != row[2]:
    else:
        Derived.append(prev_nuc.upper() + row[2] + next_nuc.upper())
#    else:
#        raise Exception('I do not know which base in anc/der')


df['Ancestral'] = Ancestral
df['Derived'] = Derived

noTCC_index = df[(df['Ancestral'] == 'TCC') & (df['Derived'] == 'TTC')].index
df.drop(noTCC_index, inplace = True)

noAGG_index = df[(df['Ancestral'] == 'AGG') & (df['Derived'] == 'AAG')].index
df.drop(noAGG_index, inplace = True)

df.to_csv(out, sep='\t', index=False)

