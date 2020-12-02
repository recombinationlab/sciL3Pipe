import argparse
import pysam
from collections import OrderedDict, defaultdict
import re
import os

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description =
            'Split alignments to a hybrid assembly consisting of ucsc and ensembl' +
            'style chromosome annotation.')
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-1', '--chr_prefix', action = 'store', metavar = 'STR',
                        help = 'chr assembly prefix')
    parser.add_argument('-2', '--prefix', action = 'store', metavar = 'STR',
                        help = 'non-chr assembly prefix')
        
    args = parser.parse_args(args)

    return args


def main():

    args = parse_arguments()

    # args = parse_arguments('-i /mnt/data/Projects/sciStrand-seq/yi293_AGGACG.PE.bwa.hg38.markdup.bam ' \
    #                        '-1 mouse -2 human'.split())

    write_split_bams(args)


def get_split_bam_headers(args):
    '''Make new BAM headers for split BAM files

        Note:
            SN - chromosome name
            LN - chromosome size
            Assumes main chromosomes come first followed by alt, that will be removed 

        Return:
            touple with header dicts for pysam
    '''
    #get bam header
    with pysam.AlignmentFile(args.input, 'rb') as input_file:
        bam_header = input_file.header.to_dict()
    #split bam header

    SQ_1 = []
    SQ_2 = []
    for chrom in bam_header['SQ']:
        # Exclude alt chromosomes (anything that doesn't end in .1 .2 etc.)
        if chrom['SN'].startswith('chr') and not bool(re.search('\.\d$', chrom['SN'])):
            SQ_1.append(chrom)
        elif not bool(re.search('\.\d$', chrom['SN'])):
            SQ_2.append(chrom)

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
    # #1 contains chr UCSC style #2 ensembl style
    bh_1, bh_2 = get_split_bam_headers(args)

    bh_1_tid_pos = defaultdict()
    for i, item in enumerate(bh_1['SQ']):
        bh_1_tid_pos[item['SN']] = i

    bh_2_tid_pos = defaultdict()
    for i, item in enumerate(bh_2['SQ']):
        bh_2_tid_pos[item['SN']] = i

    b1 = '.'.join([os.path.splitext(args.input)[0], args.chr_prefix, 'bam'])
    b2 = '.'.join([os.path.splitext(args.input)[0], args.prefix, 'bam'])

    written_1 = 0
    written_2 = 0
    not_written = 0
    with pysam.AlignmentFile(args.input, "rb") as input_file, \
        pysam.AlignmentFile(b1, "wb", header = bh_1) as b1_out, \
        pysam.AlignmentFile(b2, "wb", header = bh_2) as b2_out:
        
        for read in input_file.fetch():
            # cross reference reads (mouse human) will be filtered out
            if all([read.reference_name.startswith('chr'),
                    read.next_reference_name.startswith('chr'),
                    not bool(re.search('\.\d$', read.reference_name)),
                    not bool(re.search('\.\d$', read.next_reference_name))]):
                # update tid
                read.reference_id = bh_1_tid_pos[read.reference_name]
                # update mate tid
                read.next_reference_id = bh_1_tid_pos[read.next_reference_name]
                
                b1_out.write(read)
                written_1 += 1
            elif all([not read.reference_name.startswith('chr'),
                      not read.next_reference_name.startswith('chr'),
                      not bool(re.search('\.\d$', read.reference_name)),
                      not bool(re.search('\.\d$', read.next_reference_name))]):
                # update tid
                read.reference_id = bh_2_tid_pos[read.reference_name]
                # update mate tid
                read.next_reference_id = bh_2_tid_pos[read.next_reference_name]

                b2_out.write(read)
                written_2 += 1
            else:
                not_written += 1


    print('BAM 1 reads written out:', written_1)
    print('BAM 2 reads written out:', written_2)
    print('Reads not written out:', not_written)

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