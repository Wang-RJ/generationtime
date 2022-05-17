#!/usr/bin/env python3

"""
This script extract population-specific annotations from a genome 
variation format (GVF) file.
"""

from argparse import ArgumentParser
import numpy as np
import pandas as pd
import urllib.request, urllib.parse, urllib.error


def parseGFFAttributes(attributeString):
    """Parse the GFF3 attribute column and return a dict"""
    if attributeString == ".": return {}
    ret = {}
    for attribute in attributeString.split(';'):
        key, value = attribute.split('=')
        ret[urllib.parse.unquote(key)] = urllib.parse.unquote(value)
    return ret


def parseGVF(row):
    """
    A simple GVF format parser, modified from a GFF3 format parser.
    Return a dictionary per row for our dataframe.
    """
    
    populations = ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    
    record = {}
    
    parts = row.strip().split('\t')

    record['Chr'] = parts[0]
    record['Pos'] = parts[3]
    
    attr = parseGFFAttributes(parts[8])
    record['ID'] = attr['ID']
    dbSNP = attr['Dbxref'].split(':')
    record['rs'] = dbSNP[1]
    #maf = attr['global_minor_allele_frequency']
    record['Ref'] = attr['Reference_seq']
    record['Alt'] = attr['Variant_seq']
    
    for pop in populations:
        if pop in attr.keys(): 
            record[pop] = attr[pop] 
        else: 
            record[pop] = 0.0
            
    return record


if __name__ == '__main__':
    parser = ArgumentParser(description='Parse a file in GVF format and create a dataframe')
    parser.add_argument('-i', '--infile', help='Input a file in GVF format')
    parser.add_argument('-o', '--outfile', help='Outputfile name and path')
    args = parser.parse_args()

    fn = args.infile
    out = args.outfile

    rows_list = []
    with open(fn) as infile:
        for row in infile:
            if row.startswith('#'):
                continue
            else:    
                rows_list.append(parseGVF(row))

    names = ['Chr', 'ID', 'rs', 'Pos', 'Ref', 'Alt', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS']
    df = pd.DataFrame(rows_list, columns=names) 

    df.to_csv(out, sep='\t', index=False)


