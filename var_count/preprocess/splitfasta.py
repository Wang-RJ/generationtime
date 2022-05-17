#!/usr/bin/env python3

"""
This script extract per chromosome sequence from the GRCh37 
reference genome downloaded from:
https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.25/

It takes the ref genome file as its only argument.
"""

from Bio import SeqIO
import sys

RefSeq = ['NC_000001.10', 'NC_000002.11'
'NC_000003.11', 'NC_000004.11', 'NC_000005.9', 'NC_000006.11',
'NC_000007.13', 'NC_000008.10', 'NC_000009.11', 'NC_000010.10',
'NC_000011.9', 'NC_000012.11', 'NC_000013.10', 'NC_000014.8',
'NC_000015.9', 'NC_000016.9', 'NC_000017.10', 'NC_000018.9',
'NC_000019.9', 'NC_000020.10', 'NC_000021.8', 'NC_000022.10']

#wantedSeqs = sys.argv[2]
seqs = SeqIO.parse(open(sys.argv[1]), 'fasta')
for name in RefSeq:
    SeqIO.write((seq for seq in seqs if seq.id in name), sys.stdout, 'fasta')

