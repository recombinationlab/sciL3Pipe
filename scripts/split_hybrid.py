import argparse
import pysam
from collections import OrderedDict, defaultdict
from natsort import index_natsorted
import re
import os

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description =
            'Split alignments to a hybrid assembly consisting of ucsc and ensembl' +
            'style chromosome annotation.')
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output_dir', action = 'store', metavar = 'FILE',
                        help = 'Output directory')
    parser.add_argument('-1', '--chr_prefix', action = 'store', metavar = 'STR',
                        help = 'chr assembly prefix')
    parser.add_argument('-2', '--prefix', action = 'store', metavar = 'STR',
                        help = 'non-chr assembly prefix')
    parser.add_argument('-f', '--filter', action = 'store', metavar = 'STR',
                        help = 'Additional chromosomes to exclude, if more than one add' +
                               'in a comma separated list')
        
    args = parser.parse_args(args)

    return args


def main():

    args = parse_arguments()
    # /mnt/data/Projects/sciStrand-seq/yi293_AGGACG.PE.bwa.hg38.markdup.bam
    # args = parse_arguments('-i /mnt/data/nextseq190419/yi293_AGGACG_5min_UV_with_USER_HAP1/yi293_AGGACG.GTCAGTAGCGGAGAC.bam ' \
    #                        '-1 mouse -2 human -f NC_007605,hs37d5,Y'.split())
    os.makedirs(args.output_dir, exist_ok=True)
    write_split_bams(args)


def get_split_bam_headers(args):
    '''Make new BAM headers for split BAM files

        Note:
            SN - chromosome name
            LN - chromosome size
           
        Return:
            tuple with header dicts for pysam
    '''
    #get bam header
    with pysam.AlignmentFile(args.input, 'rb') as input_file:
        bam_header = input_file.header.to_dict()

    # Additional filter
    if args.filter:
        search_str = '^' + '$|^'.join(args.filter.split(',')) + '$'

    #split bam header
    
    SQ_1 = []
    SQ_2 = []
    for chrom in bam_header['SQ']:
        if args.filter:
            ad_filter = bool(re.search(search_str, chrom['SN']))
        else:
            ad_filter = False 
        # Exclude alt chromosomes (anything that doesn't end in .1 .2 _ etc.)
        if all([chrom['SN'].startswith('chr'),
                not bool(re.search('_|\.\d$', chrom['SN'])),
                not ad_filter]):
            SQ_1.append(chrom)
        elif not bool(re.search('_|\.\d$', chrom['SN'])) and not ad_filter:
            SQ_2.append(chrom)

    # Natural sort to order chromosomes (ensembl fastq not natural sorted)
    SQ_1 = [SQ_1[i] for i in index_natsorted([i['SN'] for i in SQ_1])]
    SQ_2 = [SQ_2[i] for i in index_natsorted([i['SN'] for i in SQ_2])]

    bam_1 = OrderedDict()
    bam_2 = OrderedDict()

    # Key order = HD SQ RG PG
    for k in bam_header.keys():
        # Add all fields before SQ
        if k == 'SQ':
            bam_1[k] = SQ_1
            bam_2[k] = SQ_2
        else:
            bam_1[k] = bam_header[k]
            bam_2[k] = bam_header[k]

    return (bam_1, bam_2)

def write_split_bams(args):
    '''
    Write out split BAM files based on presense of chr in reference
    Cross reference reads (mouse human) will be filtered out along with
    alt chromosomes

    Note:
        Need to change tid to fit new header
        tid - The target id. The target id is 0 or a positive integer mapping to 
        entries within the sequence dictionary in the header section of a TAM file or BAM file.
    '''
    # bam header #1 contains chr UCSC style #2 ensembl style
    bh_1, bh_2 = get_split_bam_headers(args)

    bh_1_tid_pos = defaultdict()
    for i, item in enumerate(bh_1['SQ']):
        bh_1_tid_pos[item['SN']] = i

    bh_2_tid_pos = defaultdict()
    for i, item in enumerate(bh_2['SQ']):
        bh_2_tid_pos[item['SN']] = i

    b1 = '.'.join([os.path.splitext(os.path.basename(args.input))[0], args.chr_prefix, 'bam'])
    b2 = '.'.join([os.path.splitext(os.path.basename(args.input))[0], args.prefix, 'bam'])
    b1 = os.path.join(args.output_dir, b1)
    b2 = os.path.join(args.output_dir, b2)

    # Additional filter
    if args.filter:
        search_str = '^' + '$|^'.join(args.filter.split(',')) + '$'

    written_1 = 0
    written_2 = 0
    not_written = 0
    orphan_reads = 0
    additional_filt = 0
    with pysam.AlignmentFile(args.input, "rb") as input_file, \
        pysam.AlignmentFile(b1, "wb", header = bh_1) as b1_out, \
        pysam.AlignmentFile(b2, "wb", header = bh_2) as b2_out:
        
        for read in input_file.fetch():
            try:
                if args.filter:
                    ad_filter = any([bool(re.search(search_str, read.reference_name)),
                                     bool(re.search(search_str, read.next_reference_name))])
                    if ad_filter:
                        additional_filt += 1
                else:
                    ad_filter = False 
                # cross reference reads (mouse human) will be filtered out
                if all([read.reference_name.startswith('chr'),
                        read.next_reference_name.startswith('chr'),
                        not bool(re.search('_|\.\d$', read.reference_name)),
                        not bool(re.search('_|\.\d$', read.next_reference_name)),
                        not ad_filter]):
                    # update tid
                    read.reference_id = bh_1_tid_pos[read.reference_name]
                    # update mate tid
                    read.next_reference_id = bh_1_tid_pos[read.next_reference_name]
                    
                    b1_out.write(read)
                    written_1 += 1
                elif all([not read.reference_name.startswith('chr'),
                        not read.next_reference_name.startswith('chr'),
                        not bool(re.search('_|\.\d$', read.reference_name)),
                        not bool(re.search('_|\.\d$', read.next_reference_name)),
                        not ad_filter]):
                    # update tid
                    read.reference_id = bh_2_tid_pos[read.reference_name]
                    # update mate tid
                    read.next_reference_id = bh_2_tid_pos[read.next_reference_name]

                    b2_out.write(read)
                    written_2 += 1
                else:
                    not_written += 1
            except AttributeError:
                orphan_reads += 1

    print('BAM 1 reads written out:', written_1)
    print('BAM 2 reads written out:', written_2)
    print('Reads not written out:', not_written)
    if orphan_reads > 0:
        print('Orphan reads skipped:', orphan_reads)
    if args.filter:
        print('Additional filter reads:', additional_filt)

if __name__ == '__main__':
    main()



# count = 0
# with pysam.AlignmentFile(args.input, "rb") as input_file:
#     for read in input_file.fetch(until_eof = True):
#         if count < 2:
#             print(read.reference_name)
#             print(read.tid)
#             print(read.reference_id)
#             print(read)
#             count += 1
#         else:
#             break

# def new_read(read, tid_dict):

#     a = pysam.AlignedSegment()
#     a.query_name = read.query_name
#     a.query_sequence=read.query_sequence
#     a.flag = read.flag
#     a.reference_id = tid_dict[read.reference_name]
#     a.reference_start = read.reference_start
#     a.mapping_quality = read.mapping_quality
#     a.cigar = read.cigar
#     a.next_reference_id = tid_dict[read.next_reference_name]
#     a.next_reference_start=read.next_reference_start
#     a.template_length=read.template_length
#     a.query_qualities = read.query_qualities
#     a.tags = read.tags

#     return(a)