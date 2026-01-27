#!/usr/bin/env python3

import json
import os
import re
from os.path import join
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument('--bam_dir', nargs='+', required=True,
					help='Required. the FULL path to the SSS BAM folder(s)')
parser.add_argument("--ends", action='store', default='srt.bam',
					help='File ending to limit fetching to a subset of files.')
parser.add_argument("--out_file", action='store', 
					help='Output json file path')
args = parser.parse_args()

assert args.bam_dir is not None, "please provide the path to the SSS BAM folder"


## default dictionary is quite useful!

FILES = defaultdict(lambda: defaultdict(list))

## build the dictionary with full path for each BAM file
for folder in args.bam_dir:
	for root, dirs, files in os.walk(folder):
		for f in files:
			if f.endswith(args.ends):
				full_path = join(root, f)
				# if f.count('_') == 2:
				# 	B = re.search(r"(.+)_([0-9]{3})_([A-Z]{6})", f)
				# 	sample = B.group(1)
				# 	sss = B.group(1) + "_" + B.group(3)  
				# else:
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
	for sss in FILES[sample]:
		print ("{sample} {sss} has {n} BAMs".format(sample = sample, sss = sss, n = len(FILES[sample][sss])))
print ("------------------------------------------")
print("check the sss_samples.json file for BAMs belong to each sample")
print()

js = json.dumps(FILES_sorted, indent = 4, sort_keys=True)

if args.out_file != None:
	open(args.out_file, 'w').writelines(js)
else:
	open('sss_samples.json', 'w').writelines(js)

