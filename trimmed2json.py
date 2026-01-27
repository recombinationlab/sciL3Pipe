#!/usr/bin/env python3

import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument('--fastq_dir', nargs='+',
					help='Required. The FULL path to the fastq folder(s)')
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
			if f.endswith(".trimmed.attached.fastq.gz"):
				full_path = join(root, f)
				B = re.search(r"(.+)_([A-Z]{6})", f)
				sample = B.group(1)
				sss = B.group(1) + "_" + B.group(2)
					 
				FILES[sample][sss].append(full_path)
				
				
#Make sure file are in correct order
FILES_sorted = defaultdict(lambda: defaultdict(list))

for sample in FILES.keys():
		for sss in FILES[sample]:
			FILES_sorted[sample][sss] = sorted(FILES[sample][sss])


print()
print ("total {} unique samples will be processed".format(len(FILES.keys())))
print ("------------------------------------------")
for sample in FILES.keys():
	for read in FILES[sample]:
		print ("{sample} {read} has {n} fastq".format(sample = sample, read = read, n = len(FILES[sample][read])))
print ("------------------------------------------")
print("check the trimmed.json file for fastqs belong to each sample")
print()

js = json.dumps(FILES_sorted, indent = 4, sort_keys=True)
if args.out_file == None:
	open('trimmed.json', 'w').writelines(js)
else:
	open(args.out_file, 'w').writelines(js)

