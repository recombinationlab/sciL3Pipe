

import argparse
import glob
import os
from collections import defaultdict
import re
import pandas as pd

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(description = print(__doc__))
    parser.add_argument('-i', '--log_dir', action = 'store', metavar = 'LOG', required=True,
                        help = 'Path to log directory')
    parser.add_argument('-o', '--output_summary', action = 'store', metavar = 'FILE', required=True, 
                        help='Output file path')
    parser.add_argument('-f', '--format', action = 'store', metavar = 'FORMAT', required=False, default='csv', choices=['csv','html'],
                        help='Output file format')

    args = parser.parse_args(args)


    return args



def main():

    args = parse_arguments()

    # args = parse_arguments('-i /u/project/yeastyin/chovanec/projects/nextseq250312/logs/ ' \
    #                        '-o /u/project/yeastyin/chovanec/projects/nextseq250312/logs/sss_log_summary.tsv '.split())

    df = parse_log(args)

    if args.format == 'csv':
        df.to_csv(args.output_summary, sep='\t', index=False)
    elif args.format == 'html':
        df.to_html(args.output_summary, index=False, border=0, classes="table table-striped")
    


def parse_log(args):

    pattern = r'barcode (\w+) read pair count: (\d+)'
    sss_log_dict = defaultdict(lambda: defaultdict(int))
    sss_logs = glob.glob(args.log_dir + "/*_sss.out")

    for f in sss_logs:
        lib = os.path.basename(f).split('_')[0]
        with open(f, 'r') as f_in:
            for line in f_in:
                match = re.search(pattern, line)
                # Extract the matched groups
                if match:
                    barcode = match.group(1)
                    count = int(match.group(2))
                    sss_log_dict[lib][barcode] += count


    # Get cell count per SSS barcode based on SSS split log
    sss_split_counts = defaultdict(int)
    for library, sss_counts in sss_log_dict.items():
        for sss, count in sss_counts.items():
            sss_gc_logs = glob.glob(args.log_dir + f"/{library}/{library}_{sss}*.genome_coverage.log")
            if len(sss_gc_logs):
                sss_split_counts[library + '_' + sss] = len(sss_gc_logs)


    # Flatten nested dicts into a single dict
    rows = []
    for library, sss_counts in sss_log_dict.items():
        for sss, count in sss_counts.items():
            cell_count = sss_split_counts.get(library + '_' + sss, 0)
            rows.append({'library': library, 'SSS': sss, 'count': count, '# of cells >10k': cell_count})

    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(rows)

    return df


if __name__ == '__main__':
    main()
