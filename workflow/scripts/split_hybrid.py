'''
Split alignments to a hybrid assembly consisting of ucsc and ensembl style chromosome annotation.
For co-assay, split DNA and RNA reads based RT UMI, which is GGGGGG for DNA.
'''



import argparse
import pysam
from collections import OrderedDict, defaultdict
from natsort import index_natsorted
import re
import os

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description = __doc__ )
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
    parser.add_argument('--dominant', action = 'store_true',
                        help = ('Output only the BAM with the most number of reads ' +
                        'and at least 10000 reads.'))
    parser.add_argument('--no_hybrid', action = 'store_false',
                        help = ('Skip spliting out hybrid reference based on chr prefix in chromsome names.'))
    parser.add_argument('--rna', action = 'store_true',
                        help = ('Split out DNA and RNA reads from co-assay.'))
    parser.add_argument('--natsort_off', action = 'store_true',
                        help = ('Turn off natural sorting of chromosomes. (Useful for yeast with Roman numeral chromosome names)'))
        
    args = parser.parse_args(args)

    return args




def main():

    args = parse_arguments()

    # /mnt/data/Projects/sciStrand-seq/yi293_AGGACG.PE.bwa.hg38.markdup.bam
    # /mnt/data/nextseq190419/yi293_AGGACG_5min_UV_with_USER_HAP1/yi293_AGGACG.GTCAGTAGCGGAGAC.bam
    # args = parse_arguments('-i /mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1_hg38_bowtie2/yi293_GAACCG/yi293_GAACCG.TTTGACTCTCACAGC.bam ' \
    #                        '-1 mouse -2 human -f NC_007605,hs37d5,Y --dominant -o /mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1_hg38_bowtie2/split/'.split())
    

    # args = parse_arguments('-i /mnt/d/SynologyDriveCloud/Projects/scil3ext/coassay/yi261_AACTAG.ACCATTTAAGAGACT.bam ' \
    #                        '-1 mouse -2 human -f NC_007605,hs37d5,Y -o /mnt/d/SynologyDriveCloud/Projects/scil3ext/coassay/ --rna --no_hybrid'.split())
    
    # args = parse_arguments('-i /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/yy401_CGCTTG.PE.bwa.hg19.markdup.bam ' \
    #                        '-1 mouse -2 human -o /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/ --filter MT'.split())
    
    if args.no_hybrid == True and args.rna == True: # split out hybrid and DNA RNA
        write_split_bams_hybrid_rna_dna(args)
    elif args.no_hybrid == False and args.rna == True: # split out DNA RNA only
        write_split_bams_rna_dna(args)
    elif args.no_hybrid == True and args.rna == False: # split out hybrid only
        write_split_bams_hybrid(args)
    else:
        print("Nothing to split!")


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

    if args.natsort_off == False:
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



def write_split_bams_rna_dna(args):
    '''
    Write out split BAM files based on presense of GGGGGG barcode or UMI
    Alt chromosome reads will also be filtered out

    Note:
        bc1 (tn5 or rt), bc2 (ligation), UMI (rt, GGGGGG for DNA), bc3 (SSS), UMI (IVT)
        bc1, bc2, bc3 together define a single cell
        UMI (rt) defines a mRNA molecule
        UMI (IVT) is not used
    '''
    # output file names
    b1_dna = '.'.join([os.path.splitext(os.path.basename(args.input))[0], 'DNA', 'bam'])
    b1_rna = '.'.join([os.path.splitext(os.path.basename(args.input))[0], 'RNA', 'bam'])
    
    b1_dna = os.path.join(args.output_dir, 'DNA', b1_dna)
    b1_rna = os.path.join(args.output_dir, 'RNA', b1_rna)

    # create dir if they don't exist
    os.makedirs(os.path.dirname(b1_dna), exist_ok=True)
    os.makedirs(os.path.dirname(b1_rna), exist_ok=True)


    # Additional filter
    if args.filter:
        search_str = '^' + '$|^'.join(args.filter.split(',')) + '$'

    # check for BAM index
    with pysam.AlignmentFile(args.input, 'rb') as input_file:
        try:
            input_file.check_index()
        except ValueError:
            print('Creating BAM index.')
            pysam.index(args.input)
        

    written_dna = 0
    written_rna = 0
    not_written = 0
    orphan_reads = 0
    additional_filt = 0
    with pysam.AlignmentFile(args.input, 'rb') as input_file, \
        pysam.AlignmentFile(b1_dna, 'wb', template = input_file) as b1_dna_out, \
        pysam.AlignmentFile(b1_rna, 'wb', template = input_file) as b1_rna_out:

        for read in input_file.fetch():
            try:
                if args.filter:
                    ad_filter = any([bool(re.search(search_str, read.reference_name)),
                                    bool(re.search(search_str, read.next_reference_name))])
                    if ad_filter:
                        additional_filt += 1
                else:
                    ad_filter = False

                # determine if read is RNA or DNA
                read_name_barcodes = read.query_name.split(',')
                if len(read_name_barcodes) != 6:
                    raise ValueError('Number of barcodes in read name does not match expectation.')
                else:
                    umi_rt = read_name_barcodes[2]
                    if umi_rt == 'GGGGGG':
                        mol_type = 'DNA'
                    else:
                        mol_type = 'RNA'

                # cross reference reads (mouse human) will be filtered out
                if all([not bool(re.search('_|\.\d$', read.reference_name)),
                        not bool(re.search('_|\.\d$', read.next_reference_name)),
                        not ad_filter]):
                    
                    if mol_type == 'DNA':
                        b1_dna_out.write(read)
                        written_dna += 1
                    else:
                        b1_rna_out.write(read)
                        written_rna += 1

                else:
                    not_written += 1
            except TypeError:
                orphan_reads += 1

    print('DNA reads:', written_dna)
    print('RNA reads:', written_rna)
    print('Reads filtered:', not_written)
    if orphan_reads > 0:
        print('Orphan reads skipped:', orphan_reads)
    if args.filter:
        print('Set filter removed reads:', additional_filt)




def write_split_bams_hybrid_rna_dna(args):
    '''
    Write out split BAM files based on presense of chr in reference and 
    based on presense of GGGGGG barcode or UMI
    Cross reference reads (mouse human) will be filtered out along with
    alt chromosomes

    Note:
        Need to change tid to fit new header
        tid - The target id. The target id is 0 or a positive integer mapping to 
        entries within the sequence dictionary in the header section of a TAM file or BAM file.

        bc1 (tn5 or rt), bc2 (ligation), UMI (rt, GGGGGG for DNA), bc3 (SSS), UMI (IVT)
        bc1, bc2, bc3 together define a single cell
        UMI (rt) defines a mRNA molecule
        UMI (IVT) is not used
    '''
    # bam header #1 contains chr UCSC style #2 ensembl style
    bh_1, bh_2 = get_split_bam_headers(args)

    bh_1_tid_pos = defaultdict()
    for i, item in enumerate(bh_1['SQ']):
        bh_1_tid_pos[item['SN']] = i

    bh_2_tid_pos = defaultdict()
    for i, item in enumerate(bh_2['SQ']):
        bh_2_tid_pos[item['SN']] = i

    # output file names
    b1_dna = '.'.join([os.path.splitext(os.path.basename(args.input))[0], args.chr_prefix, 'DNA', 'bam'])
    b2_dna = '.'.join([os.path.splitext(os.path.basename(args.input))[0], args.prefix, 'DNA', 'bam'])
    b1_rna = '.'.join([os.path.splitext(os.path.basename(args.input))[0], args.chr_prefix, 'RNA', 'bam'])
    b2_rna = '.'.join([os.path.splitext(os.path.basename(args.input))[0], args.prefix, 'RNA', 'bam'])
    
    b1_dna = os.path.join(args.output_dir, args.chr_prefix, 'DNA', b1_dna)
    b2_dna = os.path.join(args.output_dir, args.prefix, 'DNA', b2_dna)
    b1_rna = os.path.join(args.output_dir, args.chr_prefix, 'RNA', b1_rna)
    b2_rna = os.path.join(args.output_dir, args.prefix, 'RNA', b2_rna)

    # create dir if they don't exist
    os.makedirs(os.path.dirname(b1_dna), exist_ok=True)
    os.makedirs(os.path.dirname(b2_dna), exist_ok=True)
    os.makedirs(os.path.dirname(b1_rna), exist_ok=True)
    os.makedirs(os.path.dirname(b2_rna), exist_ok=True)


    # Additional filter
    if args.filter:
        search_str = '^' + '$|^'.join(args.filter.split(',')) + '$'

    # check for BAM index
    with pysam.AlignmentFile(args.input, 'rb') as input_file:
        try:
            input_file.check_index()
        except ValueError:
            print('Creating BAM index.')
            pysam.index(args.input)
        

    written_dna_1 = 0
    written_dna_2 = 0
    written_rna_1 = 0
    written_rna_2 = 0
    not_written = 0
    orphan_reads = 0
    additional_filt = 0
    with pysam.AlignmentFile(args.input, 'rb') as input_file, \
        pysam.AlignmentFile(b1_dna, 'wb', header = bh_1) as b1_dna_out, \
        pysam.AlignmentFile(b2_dna, 'wb', header = bh_2) as b2_dna_out, \
        pysam.AlignmentFile(b1_rna, 'wb', header = bh_1) as b1_rna_out, \
        pysam.AlignmentFile(b2_rna, 'wb', header = bh_2) as b2_rna_out:

        for read in input_file.fetch():
            try:
                if args.filter:
                    ad_filter = any([bool(re.search(search_str, read.reference_name)),
                                    bool(re.search(search_str, read.next_reference_name))])
                    if ad_filter:
                        additional_filt += 1
                else:
                    ad_filter = False

                # determine if read is RNA or DNA
                read_name_barcodes = read.query_name.split(',')
                if len(read_name_barcodes) != 6:
                    raise ValueError('Number of barcodes in read name does not match expectation.')
                else:
                    umi_rt = read_name_barcodes[2]
                    if umi_rt == 'GGGGGG':
                        mol_type = 'DNA'
                    else:
                        mol_type = 'RNA'

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
                    
                    if mol_type == 'DNA':
                        b1_dna_out.write(read)
                        written_dna_1 += 1
                    else:
                        b1_rna_out.write(read)
                        written_rna_1 += 1

                elif all([not read.reference_name.startswith('chr'),
                        not read.next_reference_name.startswith('chr'),
                        not bool(re.search('_|\.\d$', read.reference_name)),
                        not bool(re.search('_|\.\d$', read.next_reference_name)),
                        not ad_filter]):
                    # update tid
                    read.reference_id = bh_2_tid_pos[read.reference_name]
                    # update mate tid
                    read.next_reference_id = bh_2_tid_pos[read.next_reference_name]

                    if mol_type == 'DNA':
                        b2_dna_out.write(read)
                        written_dna_2 += 1
                    else:
                        b2_rna_out.write(read)
                        written_rna_2 += 1
                else:
                    not_written += 1
            except TypeError:
                orphan_reads += 1

    print('BAM 1 DNA reads:', written_dna_1)
    print('BAM 2 DNA reads:', written_dna_2)
    print('BAM 1 RNA reads:', written_rna_1)
    print('BAM 2 RNA reads:', written_rna_2)
    print('Reads filtered:', not_written)
    if orphan_reads > 0:
        print('Orphan reads skipped:', orphan_reads)
    if args.filter:
        print('Set filter removed reads:', additional_filt)

    if args.dominant:
        total_out = written_dna_1 + written_dna_2
        if written_dna_1/total_out <= 0.3 or written_dna_1 < 10000:
            print('DNA BAM 1 has fewer than 30% of total reads or less than 10000 reads and will be deleted.')
            os.remove(b1_dna)
        if written_dna_2/total_out <= 0.3 or written_dna_2 < 10000:
            print('DNA BAM 2 has fewer than 30% of total reads or less than 10000 reads and will be deleted.')
            os.remove(b2_dna)
        total_out = written_rna_1 + written_rna_2
        if written_rna_1/total_out <= 0.3 or written_rna_1 < 10000:
            print('RNA BAM 1 has fewer than 30% of total reads or less than 10000 reads and will be deleted.')
            os.remove(b1_rna)
        if written_rna_2/total_out <= 0.3 or written_rna_2 < 10000:
            print('RNA BAM 2 has fewer than 30% of total reads or less than 10000 reads and will be deleted.')
            os.remove(b2_rna)





def write_split_bams_hybrid(args):
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

    # output file names
    b1 = '.'.join([os.path.splitext(os.path.basename(args.input))[0], args.chr_prefix, 'bam'])
    b2 = '.'.join([os.path.splitext(os.path.basename(args.input))[0], args.prefix, 'bam'])
    b1 = os.path.join(args.output_dir, args.chr_prefix, b1)
    b2 = os.path.join(args.output_dir, args.prefix, b2)

    # create dir if they don't exist
    os.makedirs(os.path.dirname(b1), exist_ok=True)
    os.makedirs(os.path.dirname(b2), exist_ok=True)


    # Additional filter
    if args.filter:
        search_str = '^' + '$|^'.join(args.filter.split(',')) + '$'

    # check for BAM index
    with pysam.AlignmentFile(args.input, 'rb') as input_file:
        try:
            input_file.check_index()
        except ValueError:
            print('Creating BAM index.')
            pysam.index(args.input)
        

    written_1 = 0
    written_2 = 0
    not_written = 0
    orphan_reads = 0
    additional_filt = 0
    with pysam.AlignmentFile(args.input, 'rb') as input_file, \
        pysam.AlignmentFile(b1, 'wb', header = bh_1) as b1_out, \
        pysam.AlignmentFile(b2, 'wb', header = bh_2) as b2_out:

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
            except TypeError:
                orphan_reads += 1

    print('BAM 1 reads:', written_1)
    print('BAM 2 reads:', written_2)
    print('Reads filtered:', not_written)
    if orphan_reads > 0:
        print('Orphan reads skipped:', orphan_reads)
    if args.filter:
        print('Set filter removed reads:', additional_filt)

    if args.dominant:
        total_out = written_1 + written_2
        if written_1/total_out <= 0.3 or written_1 < 10000:
            print('BAM 1 has fewer than 30% of total reads or less than 10000 reads and will be deleted.')
            os.remove(b1)
        if written_2/total_out <= 0.3 or written_2 < 10000:
            print('BAM 2 has fewer than 30% of total reads or less than 10000 reads and will be deleted.')
            os.remove(b2)


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