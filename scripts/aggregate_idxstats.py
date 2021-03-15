import pysam
import argparse
import glob
import os
from collections import defaultdict
import re
from natsort import index_natsorted
import pandas as pd

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description =
            'Validate Paired-end BAM file, making sure it has proper pairs, ' +
            'used to QC other scripts.')
    parser.add_argument('-d', '--dir', action = 'store', metavar = 'FILE',
                        nargs='+', help = 'Input BAM directory. Allows multiple entries.')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Summary file output')
    parser.add_argument('-f', '--format', action = 'store', metavar = 'STR',
                        choices=['UCSC', 'ensembl', 'none'], default='none',
                        help = 'UCSC (mouse) or ensembl (human) format to extract')
    parser.add_argument('-r', '--regex', action = 'store', metavar = 'STR',
                        default='*.bam',
                        help = 'Regex to specify which BAM files to process in folder, ' +
                               'default: "*.bam"')
        
    args = parser.parse_args(args)

    return args


def main():

    args = parse_arguments()

    # args = parse_arguments('-d /mnt/data/sci-l3/test/split_bam_bwa/yi292/yi292_GTGCAG ' \
    #                        '/mnt/data/sci-l3/test/split_bam_bwa/yi292/yi292_CACGTG ' \
    #                        '-o /mnt/data/sci-l3/test/idxstats_summary.tsv '\
    #                        '--format none'.split())

    summary = aggregate_idxstats(args)

    summary.to_csv(args.output, index=True, sep='\t')

    # test = get_idxstats('/mnt/data/sci-l3/test/split_bam_bwa/yi292/yi292_GTGCAG/yi292_GTGCAG.CCGCGCTTGGCCAAA.bam', args.format)

def aggregate_idxstats(args):
    '''
    Loop over a directory of BAM files and aggregate idxstats into a pandas dataframe
    '''
    idx_summary = defaultdict()
    col_names = []

    for directory in args.dir:
        for f in glob.glob(os.path.join(directory, args.regex)):
            bname = os.path.basename(f)
            # print(f'Getting stats for: {bname}')

            col_nam, row = get_idxstats(f, args.format)
            idx_summary[bname] = row

            if len(col_names) == 0:
                col_names = col_nam
            else:
                assert col_names == col_nam, 'Column names do not match'

    df = pd.DataFrame.from_dict(idx_summary,orient='index',
                                columns=col_names)

    return df


def get_idxstats(file_path, assembly_format):
    '''
    Run samtools idxstats and parse

    reference sequence name, 
    sequence length, 
    # mapped read-segments and 
    # unmapped read-segments
    '''
    idxstats = pysam.idxstats(file_path)

    chroms = []
    mapped_reads = []
    # Exclude alt chromosomes (anything that doesn't end in .1 .2 _ etc.)
    
    for line in idxstats.rstrip().split('\n'):
        chrom, _, mapped, _ = line.split('\t')

        # filter alt chromosomes
        if all([assembly_format == 'UCSC',
                chrom.startswith('chr'),
                not bool(re.search('_|\.\d$', chrom))]):
            chroms.append(chrom)
            mapped_reads.append(mapped)
        elif all([assembly_format == 'ensembl',
                  not line.startswith('chr'),
                  not bool(re.search('_|\.\d$', chrom))]):
            chroms.append(chrom)
            mapped_reads.append(mapped)
        # output everything in idxstats file
        elif assembly_format == 'none':
            chroms.append(chrom)
            mapped_reads.append(mapped)

    rows = [mapped_reads[i] for i in index_natsorted(chroms)]
    column_names = [chroms[i] for i in index_natsorted(chroms)]

    return [column_names, rows]

if __name__ == '__main__':
    main()