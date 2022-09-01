#!/usr/bin/env python

"""
This script reads a generation time dataframe per chromosome and a bed 
file with replication timings. It outputs a single column with rep
timings to be added later to put dataframe.
"""

import sys
import numpy as np
import pandas as pd

# reading only the "Position" from our df
sites = pd.read_csv(sys.argv[1], sep='\t', usecols=[1], header=0, names=['Position']).values
bed = pd.read_csv(sys.argv[2], sep='\t', usecols=[2, 3], names=['point', 'time'])


times = []
prev_pt = 1     # starting from position 1
prev_time = bed['time'][0]  # position 1 takes first rep time
idx = 0
count = 0
for r, row in bed.iterrows():
	for s, site in enumerate(sites[idx:]):
		if row['time'] == None or prev_time == None:
			times.append('NaN')
		else:	
			diff = (int(row['point']) - int(prev_pt)) / 2
			reg1 = int(prev_pt) + int(diff)
			if int(site[0]) > int(prev_pt) and int(site[0]) <= int(reg1):
				times.append(prev_time)
				count += 1
			elif int(site[0]) > int(reg1) and int(site[0]) <= int(row['point']):
				times.append(row['time'])
				count += 1
			elif int(site[0]) > int(row['point']):
				prev_pt = int(row['point'])
				prev_time = row['time']
				idx = count + 1
				break
			else:
				raise Exception('something is wrong!')


# we want to fill the last few positions with rates
while len(times) <= len(sites):
	times.append(row['time'])

# print result
for i, j in zip(sites, times):
	print(f"{i}\t{j}")


