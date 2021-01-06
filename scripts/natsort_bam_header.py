import argparse
import pysam
from collections import defaultdict
from natsort import index_natsorted
import re
import os


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description =
            'Natural sort BAM header')
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output_dir', action = 'store', metavar = 'FILE',
                        help = 'Output directory')
    parser.add_argument('--keep_alt', action = 'store_true',
                        help = ('Keep alt contigs'))
        
    args = parser.parse_args(args)

    return args



def main():

    args = parse_arguments()

    # args = parse_arguments('-i /mnt/data/nextseq190419/yi293_AGGACG_5min_UV_with_USER_HAP1/yi293_AGGACG.GTCAGTAGCGGAGAC.bam'.split())
    # natsort_header(args)

    write_bam(args)


def natsort_header(args):
    '''Natural sort BAM header
    '''
    # get bam header
    with pysam.AlignmentFile(args.input, 'rb') as input_file:
        bam_header = input_file.header.to_dict()

    if not args.keep_alt:
        SQ = []
        for chrom in bam_header['SQ']:
            # Exclude alt chromosomes (anything that doesn't end in .1 .2 _ etc.)
            if not bool(re.search('_|\.\d$', chrom['SN'])):
                SQ.append(chrom)
        
        SQ = [SQ[i] for i in index_natsorted([chrom['SN'] for chrom in SQ])]
    else:
        SQ = [bam_header['SQ'][i] for i in index_natsorted([chrom['SN'] for chrom in bam_header['SQ']])]
        
    # Natural sort to order chromosomes (ensembl fastq not natural sorted)
    bam_header['SQ'] = SQ

    return bam_header


def write_bam(args):
    '''
    Write out BAM with new header with alt chromosomes optionally filtered out

    Note:
        Need to change tid to fit new header
        tid - The target id. The target id is 0 or a positive integer mapping to 
        entries within the sequence dictionary in the header section of a TAM file or BAM file.
    '''
    # bam header #1 contains chr UCSC style #2 ensembl style
    sorted_bh = natsort_header(args)

    sbh_tid_pos = defaultdict()
    for i, item in enumerate(sorted_bh['SQ']):
        sbh_tid_pos[item['SN']] = i

    b = '.'.join([os.path.splitext(os.path.basename(args.input))[0], 'hns.bam'])
    b = os.path.join(args.output_dir, b)
    
    
    written_1 = 0
    not_written = 0
    orphan_reads = 0

    with pysam.AlignmentFile(args.input, "rb") as input_file, \
        pysam.AlignmentFile(b, "wb", header = sorted_bh) as b_out:
        
        for read in input_file.fetch():
            try:
                if not args.keep_alt:
                    if all([not bool(re.search('_|\.\d$', read.reference_name)),
                            not bool(re.search('_|\.\d$', read.next_reference_name))]):
                        # update tid
                        read.reference_id = sbh_tid_pos[read.reference_name]
                        # update mate tid
                        read.next_reference_id = sbh_tid_pos[read.next_reference_name]
                        
                        b_out.write(read)
                        written_1 += 1
                    else:
                        not_written += 1
                else:
                    # update tid
                    read.reference_id = sbh_tid_pos[read.reference_name]
                    # update mate tid
                    read.next_reference_id = sbh_tid_pos[read.next_reference_name]
                        
                    b_out.write(read)
                    written_1 += 1
            except AttributeError:
                orphan_reads += 1

    print('Reads written out:', written_1)
    print('Reads not written out:', not_written)
    if orphan_reads > 0:
        print('Orphan reads skipped:', orphan_reads)
