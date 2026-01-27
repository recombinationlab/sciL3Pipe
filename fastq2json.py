#!/usr/bin/env python3

import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict

'''
Modified from https://github.com/crazyhottommy/pyflow-RNAseq/blob/master/fastq2json.py

'''

parser = argparse.ArgumentParser()
parser.add_argument('--fastq_dir', nargs='+',
					help='Required. The FULL path to the fastq folder(s)')
parser.add_argument('--split', action='store_true',
					help='If input fastqs have been split into many parts with seqkit split2')
parser.add_argument('-o', '--out_file', action='store', 
					help='Optional. Name of output file, otherwise default is samples.json.')
args = parser.parse_args()

assert args.fastq_dir is not None, "Please provide the path to the fastq folder"


## default dictionary is quite useful!

FILES = defaultdict(lambda: defaultdict(list))

## build the dictionary with full path for each fastq.gz file
for folder in args.fastq_dir:
	for root, dirs, files in os.walk(folder):
		for f in files:
			if f.endswith("fastq.gz") or f.endswith("fq.gz") or f.endswith("fastq") or f.endswith("fq"):
				full_path = join(root, f)
				if '_I1_' in f:
					continue
				if '_00' in f:
					if 'L00' in f:
						if args.split:
							PE = re.search(r"(.+)_(L[0-9]{3})_(R[12])_00[0-9].part_([0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
							SE = re.search(r"(.+)_(L[0-9]{3})_00[0-9].part_([0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
							split_g = 4
						else:
							PE = re.search(r"(.+)_(L[0-9]{3})_(R[12])_00[0-9].(fastq.gz|fq.gz|fastq|fq)", f)
							SE = re.search(r"(.+)_(L[0-9]{3})_00[0-9].(fastq.gz|fq.gz|fastq|fq)", f)
						reads_g = 3
					else:
						if args.split:
							PE = re.search(r"(.+)_(R[12])_00[0-9].part_([0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
							SE = re.search(r"(.+)_00[0-9].part_([0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
							split_g = 3
						else:
							PE = re.search(r"(.+)_(R[12])_00[0-9].(fastq.gz|fq.gz|fastq|fq)", f)
							SE = re.search(r"(.+)_00[0-9].(fastq.gz|fq.gz|fastq|fq)", f)
						reads_g = 2
				else:	
					if 'L00' in f:
						if args.split:
							PE = re.search(r"(.+)_(L[0-9]{3})_(R[12]).part_([0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
							SE = re.search(r"(.+)_(L[0-9]{3}).part_([0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
							split_g = 4
						else:
							PE = re.search(r"(.+)_(L[0-9]{3})_(R[12]).(fastq.gz|fq.gz|fastq|fq)", f)
							SE = re.search(r"(.+)_(L[0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
						reads_g = 3
					else:
						if args.split:
							PE = re.search(r"(.+)_(R[12]).part_([0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
							SE = re.search(r"(.+).part_([0-9]{3}).(fastq.gz|fq.gz|fastq|fq)", f)
							split_g = 3
						else:
							PE = re.search(r"(.+)_(R[12]).(fastq.gz|fq.gz|fastq|fq)", f)
							SE = re.search(r"(.+).(fastq.gz|fq.gz|fastq|fq)", f)
						reads_g = 2
				if PE:
					if args.split:
						sample = PE.group(1).split('_')[0] + '_' + PE.group(split_g)
					else:
						sample = PE.group(1).split('_')[0]
						
					reads = PE.group(reads_g)  
					FILES[sample][reads].append(full_path)
				else:
					sample = SE.group(1).split('_')[0]
					FILES[sample]['R1'].append(full_path)
				
#Make sure file from different lanes are in correct order
FILES_sorted = defaultdict(lambda: defaultdict(list))

for sample in FILES.keys():
		for read in FILES[sample]:
			FILES_sorted[sample][read] = sorted(FILES[sample][read])


print()
print ("total {} unique samples will be processed".format(len(FILES.keys())))
print ("------------------------------------------")
for sample in FILES.keys():
	for read in FILES[sample]:
		print ("{sample} {read} has {n} fastq".format(sample = sample, read = read, n = len(FILES[sample][read])))
print ("------------------------------------------")
print("check the samples.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES_sorted, indent = 4, sort_keys=True)
if args.out_file == None:
	open('samples.json', 'w').writelines(js)
else:
	open(args.out_file, 'w').writelines(js)

