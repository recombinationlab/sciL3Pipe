

import argparse
import glob
import os
from collections import defaultdict
import re
import pandas as pd
from functools import reduce

def parse_dict(arg):
    '''
    Docstring for parse_dict

    Custom type function to convert key:val,key:val to a dictionary.
    
    :param arg: Arguments passed from argparse
    '''
    result = {}
    try:
        # Split by comma to get individual pairs
        for item in arg.split(','):
            key, value = item.split(':')
            result[key.strip()] = int(value.strip())
        return result
    except ValueError:
        raise argparse.ArgumentTypeError("Cutoff must be in format 'key1:value1,key2:value2' with integer values.")


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(description = print(__doc__))
    parser.add_argument('-i', '--log_dir', action = 'store', metavar = 'LOG', required=True,
                        help = 'Path to log directory')
    parser.add_argument('-o', '--output_summary', action = 'store', metavar = 'FILE', required=True, 
                        help='Output file path')
    parser.add_argument('-f', '--format', action = 'store', metavar = 'FORMAT', required=False, default='csv', choices=['csv','html'],
                        help='Output file format')
    parser.add_argument('--cutoff', 
                        type=parse_dict, 
                        default={'hs37d': 10000},
                        help='Key-value pairs for cutoffs (e.g., hg19:10000,cerv:1000)')

    args = parser.parse_args(args)

    return args

# Use the custom function as the 'type'



def main():

    args = parse_arguments()

    # args = parse_arguments('-i /u/project/yeastyin/chovanec/projects/novaseq251219/logs/ ' \
    #                        '-o /u/project/yeastyin/chovanec/projects/novaseq251219/logs/removed_read_pairs_log_summary.tsv ' \
    #                         '--cutoff hg19:10000,sacCer3:1000,GRCm38:10000'.split())

    df_1 = parse_sss_log(args)
    df_2 = parse_tn5_lig_log(args)
    df_3 = parse_split_bam_log(args)


    # List your dataframes
    dfs = [df_1, df_2, df_3]

    # Merge them all on the 'library' column
    df_combined = reduce(lambda left, right: pd.merge(left, right, on='library', how='outer'), dfs)
  
    if args.format == 'csv':
        df_combined.to_csv(args.output_summary, sep='\t', index=False)
    elif args.format == 'html':
        df_combined.to_html(args.output_summary, index=False, border=0, classes="table table-striped")
    


def parse_sss_log(args):

    '''
    Docstring for parse_sss_log
    
    Log example:
        Total read pairs processed: 10000000
        Read pairs without SSS barcode: 419186
        Read pairs without RT primer: 697096
        R1-R2 orientation read pairs written out: 4749428
        R2-R1 orientation read pairs written out: 4134290
        SSS with barcode TAATGC read pair count: 3815329
        SSS with barcode GTCTAT read pair count: 2452331
        SSS with barcode TTCGAG read pair count: 1425392
        SSS with barcode GTGCAG read pair count: 169588
        SSS with barcode CCAGTG read pair count: 289567
        SSS with barcode ACTTCG read pair count: 665385
        SSS with barcode TGAACG read pair count: 66126
        Finished

    :param args: Arguments passed from argparse
    '''

    pattern_1 = r'Total read pairs processed: (\d+)'
    pattern_2 = r'Read pairs without SSS barcode: (\d+)'
    pattern_3 = r'Read pairs without RT primer: (\d+)'

    sss_log_dict = defaultdict(lambda: defaultdict(int))
    sss_logs = glob.glob(args.log_dir + "/*_sss.out")

    for f in sss_logs:
        lib = os.path.basename(f).split('_')[0]
        with open(f, 'r') as f_in:
            for line in f_in:
                match_1 = re.search(pattern_1, line.strip())
                match_2 = re.search(pattern_2, line.strip())  
                match_3 = re.search(pattern_3, line.strip())    
                # Extract the matched groups
                if match_1:
                    total_read_pairs = int(match_1.group(1))
                    sss_log_dict[lib]['Total read pairs processed SSS'] += total_read_pairs
                elif match_2:
                    without_sss = int(match_2.group(1))
                    sss_log_dict[lib]['Read pairs without SSS barcode'] += without_sss
                elif match_3:
                    without_rt = int(match_3.group(1))
                    sss_log_dict[lib]['Read pairs without RT primer'] += without_rt

   
    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(sss_log_dict).T.reset_index().rename(columns={'index': 'library'})

    return df


def parse_tn5_lig_log(args):

    '''
    Docstring for parse_tn5_lig_log
    
    Total read pairs: 20755337
    Total read pairs written: 19404777
    Read pairs without Tn5 or ligation barcode: 237183


    :param args: Arguments passed from argparse
    '''

    pattern_1 = r'Total read pairs: (\d+)'
    pattern_2 = r'Total read pairs written: (\d+)'
    pattern_3 = r'Read pairs without Tn5 or ligation barcode: (\d+)'

    tl_log_dict = defaultdict(lambda: defaultdict(int))
    tl_logs = glob.glob(args.log_dir + "/*/*tn5_ligation_barcodes.log")

    for f in tl_logs:
        lib = os.path.basename(f).split('_')[0]
        with open(f, 'r') as f_in:
            for line in f_in:
                match_1 = re.search(pattern_1, line.strip())
                match_2 = re.search(pattern_2, line.strip())  
                match_3 = re.search(pattern_3, line.strip())    
                # Extract the matched groups
                if match_1:
                    total_read_pairs = int(match_1.group(1))
                    tl_log_dict[lib]['Total read pairs Tn5 Ligation'] += total_read_pairs
                elif match_2:
                    total_read_pairs_out = int(match_2.group(1))
                    tl_log_dict[lib]['Total read pairs written Tn5 Ligation'] += total_read_pairs_out
                elif match_3:
                    removed_rp = int(match_3.group(1))
                    tl_log_dict[lib]['Read pairs without Tn5 or ligation barcode'] += removed_rp

   
    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(tl_log_dict).T.reset_index().rename(columns={'index': 'library'})

    return df




def parse_split_bam_log(args):

    '''
    Docstring for parse_split_bam_log
    
    /u/project/yeastyin/chovanec/projects/novaseq251219/split_bam_bwa_hg19/yi492/yi492_AATTGA

    TGATATTGGCCACCT, 368
    TGATATTGGCTATTA, 64
    TGATATTGTATGTTG, 118
    TGATATTGTTAACCG, 194
    TGATATTGGGGTACA, 118
    TGATATTGTTCTATT, 162
    TGATATTGGCTTGAT, 27568
    TGATATTGTCCAAGT, 120
    TGATATTGAGAACTA, 88358
    TGATATTGTCGGTTT, 166
    TGATATTGAGAGACT, 714
    TGATATTGCTGTTAA, 150
    TGATATTGAGGTAGT, 53650
    TGATATTGGGCTTGG, 342
    TGATATTGCACTGTT, 43314

    :param args: Arguments passed from argparse
    '''

    split_bam_log_dict = defaultdict(lambda: defaultdict(int))
    logs = glob.glob(args.log_dir + "../split_bam*/*/*/*.txt")
    barcode_pattern = re.compile(r"_[ACTG]{6}\.txt$")

    split_bam_logs = [f for f in logs if not f.endswith(".sample_list.txt") and barcode_pattern.search(f)]

    for f in split_bam_logs:
        lib = os.path.basename(f).split('_')[0]
        # select cutoff from dict
        for k,v in args.cutoff.items():
            if k in f:
                cutoff = v

        with open(f, 'r') as f_in:
            for line in f_in:
                barcode, count = line.strip().split(',')

                if int(count) >= int(cutoff):
                    split_bam_log_dict[lib]['Read pair counts of cells passing cutoff'] += int(count)
                else:
                    split_bam_log_dict[lib]['Read pair counts of cells failing cutoff'] += int(count)
   
    # Convert the list of dictionaries to a DataFrame
    df = pd.DataFrame(split_bam_log_dict).T.reset_index().rename(columns={'index': 'library'})

    return df




if __name__ == '__main__':
    main()
