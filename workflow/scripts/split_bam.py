'''
This script accept a bam file, an output folder, two barcode files and a cut_off value, then it will call the
bam_barcode_count function and get the total read count per barcode, then it use the cut_off value to filter the barcode,
and generate the output samfile for single cells, generate the reads distribution in the split_bam/read_distribution_barcode.

For co-assay, the script can output BAMs with both DNA and RNA reads based on Tn5 and RT barcodes. Note, Tn5 and RT barcodes are assumed
match 1:1, i.e. Tn5 barcode 1 == RT barcode 1. For co-assays, DNA and RNA reads are counted together for the cutoff parameter. 
'''


import pysam
import matplotlib as mpl
mpl.use('Agg') # non-interactive backend required for cluster nodes
import matplotlib.pyplot as plt
import argparse
import os
from collections import OrderedDict, defaultdict
import numpy as np


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(description = print(__doc__ ))
    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output_dir', action = 'store', metavar = 'DIR',
                        help = 'Output folder')
    parser.add_argument('-tb', '--tn5_barcodes', action = 'store', metavar = 'FILE',
                        help = 'File with Tn5 barcodes')
    parser.add_argument('-rb', '--rt_barcodes', action = 'store', metavar = 'FILE',
                        help = 'File with RT barcodes')
    parser.add_argument('-lb', '--ligation_barcodes', action = 'store', metavar = 'FILE',
                        help = 'File with ligation barcodes')
    parser.add_argument('-sb', '--sss_barcodes', action = 'store', metavar = 'FILE',
                        help = 'File with second strand synthesis barcodes')
    parser.add_argument('--cutoff', metavar = 'INT', type = int, default = 10000,
                        help = ('Number of unique reads cutoff for splitting single cells' +
                                'Default: 10000'))
    parser.add_argument('--rna', action='store_true', 
                        help = 'RNA reads present, distinguished by 3rd position barcode.')
    parser.add_argument('--skip_split', action='store_true',
                        help = 'Only perform barcode filtering based on cutoff and write everything into a single BAM.')

    args = parser.parse_args(args)

    return args

def main():

    args = parse_arguments()


    # args = parse_arguments('-i /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/yi333_AAGAAA.srt.bam ' \
    #                        '-o /mnt/d/split_test/ ' \
    #                        '-tb /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/barcodes/barcode_tn5.txt ' \
    #                        '-rb /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/barcodes/barcode_rt.txt ' \
    #                        '-lb /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/barcodes/barcode_ligation.txt ' \
    #                        '--cutoff 10000 ' \
    #                        '--rna'.split())

    # args = parse_arguments('-i /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/yi293_TTCGAG.PE.bwa.hg38.markdup.bam ' \
    #                        '-o /mnt/d/split_test/ ' \
    #                        '-tb /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/barcodes/barcode_tn5.txt ' \
    #                        '-lb /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/barcodes/barcode_ligation.txt ' \
    #                        '--cutoff 10000 ' \
    #                        '--skip_split'.split())
    
    
    # args = parse_arguments('-i /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/yi333_AAGAAA.srt.bam ' \
    #                        '-o /mnt/d/split_test/ ' \
    #                        '-tb /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/barcodes/barcode_tn5.txt ' \
    #                        '-rb /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/barcodes/barcode_rt.txt ' \
    #                        '-lb /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/barcodes/barcode_ligation.txt ' \
    #                        '--cutoff 10000 ' \
    #                        '--rna ' \
    #                        '--skip_split'.split())

    # create directory tree if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

    if args.rna:
        if args.rt_barcodes == None:
           print('With RNA option, RT barcodes are required')
           exit()

        if args.skip_split:
            print('Filtering BAM with RNA and DNA')
            fitler_bam_with_rna(args)
        else:
            print('Spliting BAM with RNA and DNA')
            split_bam_with_rna(args)
    elif args.sss_barcodes != None:
        if args.skip_split:
            print('Filtering BAM with SSS')
            filter_bam_with_sss(args)
        else:
            print('Splitting BAM with SSS')
            split_bam_with_sss(args)
    elif args.skip_split:
        print('Filtering BAM')
        filter_bam(args)
    else:
        print('Spliting BAM')
        split_bam(args.input, args.output_dir, args.tn5_barcodes, args.ligation_barcodes, args.cutoff)


def bam_barcode_count(bam_file, barcode_tn5_file, barcode_ligation_file):
    
    #generate the barcode list and barcode dictionary
    barcode_tn5 = []
    with open(barcode_tn5_file) as barcodes:
        for barcode in barcodes:
            barcode_tn5.append(barcode.strip())
    
    barcode_ligation = []
    with open(barcode_ligation_file) as barcodes:
        for barcode in barcodes:
            barcode_ligation.append(barcode.strip())
    
    barcode_combination = []
    barcode_dic = {}
    for tn5 in barcode_tn5:
        for ligation in barcode_ligation:
            barcode = tn5 + ligation
            barcode_combination.append(barcode)
            barcode_dic[barcode] = 0
    
    # check for BAM index
    with pysam.AlignmentFile(bam_file, "rb") as input_sam:
        try:
            input_sam.check_index()
        except ValueError:
            pysam.index(bam_file)

    #read the sam file, and count the number per barcode
    with pysam.AlignmentFile(bam_file, "rb") as input_sam:

        for line in input_sam.fetch():
            barcode = line.query_name.split(',')[0]+line.query_name.split(',')[1]
            barcode_dic[barcode] += 1

    return barcode_dic

def split_bam(bam_file, output_folder, barcode_tn5_file, barcode_ligation_file, cut_off):
    '''
    this script accept a bam file, a sample name, an output folder, two barcode files and a cut_off value,
    then it will call the bam_barcode_count function and get the total read count per barcode,
    then it use the cut_off value to filter the barcode,
    and generate the output samfile for single cells,
    generate the reads distribution in the split_bam/read_distribution_barcode;
    '''
    
    # generate the count per barcode
    barcode_count = bam_barcode_count(bam_file, barcode_tn5_file, barcode_ligation_file)
    
    # plot the barcode reads distribution and save the result to the ouput folder
    plot_name = (bam_file.split('/')[-1]).split('.')[0]
    fig = plt.figure()
    plt.hist(barcode_count.values(), bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Number of unique reads')
    plt.axvline(x=cut_off, color='r', linestyle='dashed', linewidth=1)
    plt.yscale('log')
    fig_output = output_folder + '/' + plot_name + '.png'
    
    fig.savefig(fig_output)

    #also output the barcode number and distribution to the output folder
    with open(output_folder + '/' + plot_name + '.txt', 'w') as read_dist:
        for barcode in barcode_count:
            line = barcode + ', %d\n' %(barcode_count[barcode])
            read_dist.write(line)
    
    #Generate the read distribution in the split_bam/read_distribution_barcode
    
    #filter the barcode based on the cut_off value
    barcode_filtered = []
    for barcode in barcode_count:
        if (barcode_count[barcode] >= int(cut_off)):
            barcode_filtered.append(barcode)

    #generate the output sam file and sample_list file
    with open(output_folder + '/' + plot_name + '.' + 'sample_list.txt', 'w') as sample_list_file:
        output_files = {}
        paired_reads = {}
        for barcode in barcode_filtered:
            output_file = output_folder + '/' + plot_name + '.' + barcode + '.bam'
            output_files[barcode] = open(output_file, 'w')
            sample_list_file.write(plot_name + '.' + barcode + '\n')
    
    # output each read to the output sam file
    with pysam.AlignmentFile(bam_file, "rb") as input_sam:
        for barcode in barcode_filtered:
            outname=output_files[barcode]
            paired_reads[outname] = pysam.AlignmentFile(outname, "wb", template=input_sam)
        
        for line in input_sam.fetch(): # until_eof=True to include unaligned reads
            barcode = line.query_name.split(',')[0]+line.query_name.split(',')[1]
            if barcode in barcode_filtered:
                outname=output_files[barcode]
                paired_reads[outname].write(line)

    #close the files:
    for barcode in barcode_filtered:
        output_files[barcode].close()



def filter_bam(args):
    '''
    Only filter barcodes based on cutoff and write to a single BAM file.
    '''
    
    # generate the count per barcode
    barcode_count = bam_barcode_count(args.input, args.tn5_barcodes, args.ligation_barcodes)
    
    # plot the barcode reads distribution and save the result to the ouput folder
    plot_name = (args.input.split('/')[-1]).split('.')[0]
    fig = plt.figure()
    plt.hist(barcode_count.values(), bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Number of unique reads')
    plt.axvline(x=args.cutoff, color='r', linestyle='dashed', linewidth=1)
    plt.yscale('log')
    fig_output = args.output_dir + '/' + plot_name + '.png'
    
    fig.savefig(fig_output)

    #also output the barcode number and distribution to the output folder
    with open(args.output_dir + '/' + plot_name + '.txt', 'w') as read_dist:
        for barcode in barcode_count:
            line = barcode + ', %d\n' %(barcode_count[barcode])
            read_dist.write(line)
    
    #Generate the read distribution in the split_bam/read_distribution_barcode

    #filter the barcode based on the cut_off value
    barcode_filtered = []
    for barcode in barcode_count:
        if (barcode_count[barcode] >= int(args.cutoff)):
            barcode_filtered.append(barcode)

    out_path = os.path.join(args.output_dir, os.path.splitext(os.path.basename(args.input))[0] + '.bardcode_filt.bam')

    reads_written = 0
    reads_skipped = 0
    with pysam.AlignmentFile(args.input, "rb") as input_bam, \
         pysam.AlignmentFile(out_path, "wb", template=input_bam) as out_bam:
        for line in input_bam.fetch():
            barcode = line.query_name.split(',')[0]+line.query_name.split(',')[1]
            if barcode in barcode_filtered:
                out_bam.write(line)
                reads_written += 1
            else:
                reads_skipped += 1

    print(f'Reads written out: {reads_written}')
    print(f'Reads with low barcode counts filtered: {reads_skipped}')
    


def bam_barcode_count_with_rna(bam_file, barcode_tn5_file, barcode_ligation_file, barcode_rt_file):

    # bam_file = args.input
    # barcode_tn5_file = args.tn5_barcodes
    # barcode_ligation_file = args.ligation_barcodes
    # barcode_rt_file = args.rt_barcodes
    
    #generate the barcode list and barcode dictionary
    barcode_tn5 = []
    with open(barcode_tn5_file) as barcodes:
        for barcode in barcodes:
            barcode_tn5.append(barcode.strip())
    
    barcode_ligation = []
    with open(barcode_ligation_file) as barcodes:
        for barcode in barcodes:
            barcode_ligation.append(barcode.strip())

    barcode_rt = []
    with open(barcode_rt_file) as barcodes:
        for barcode in barcodes:
            barcode_rt.append(barcode.strip())

    assert len(barcode_tn5) == len(barcode_rt), 'Number of Tn5 and RT barcodes do not agree.'
    
    Warning('Assuming Tn5 and RT barcodes match by order.')

    tn5_barcode_dic = OrderedDict()
    rt_barcode_dic = OrderedDict()
    rt_to_tn5_conversion = OrderedDict()

    for i in range(0, len(barcode_tn5)):
        for ligation in barcode_ligation:
            tn5_bc = barcode_tn5[i] + ligation
            rt_bc = barcode_rt[i] + ligation
            tn5_barcode_dic[tn5_bc] = 0
            rt_barcode_dic[rt_bc] = 0
            rt_to_tn5_conversion[rt_bc] = tn5_bc

    # check for BAM index
    with pysam.AlignmentFile(bam_file, "rb") as input_sam:
        try:
            input_sam.check_index()
        except ValueError:
            pysam.index(bam_file)

    #read the sam file, and count the number per barcode
    with pysam.AlignmentFile(bam_file, "rb") as input_sam:

        for read in input_sam.fetch():

            # determine if read is RNA or DNA
            read_name_parts = read.query_name.split(',')
            if len(read_name_parts) != 6:
                raise ValueError('Number of barcodes in read name does not match expectation.')
            else:
                bc1 = read_name_parts[0]
                bc2 = read_name_parts[1]
                umi_rt = read_name_parts[2]
                # bc3 = read_name_parts[3]
                umi_rt = read_name_parts[2]
                if umi_rt == 'GGGGGG':
                    # this is DNA molecule
                    tn5_barcode_dic[bc1+bc2] += 1
                else:
                    # this is RNA molecule
                    rt_barcode_dic[bc1+bc2] += 1

            # barcode = line.query_name.split(',')[0]+line.query_name.split(',')[1]
            # barcode_dic[barcode] += 1

    return [tn5_barcode_dic, rt_barcode_dic, rt_to_tn5_conversion]


def split_bam_with_rna(args):
    '''
    this script accept a bam file, a sample name, an output folder, two barcode files and a cut_off value,
    then it will call the bam_barcode_count function and get the total read count per barcode,
    then it use the cut_off value to filter the barcode,
    and generate the output samfile for single cells,
    generate the reads distribution in the split_bam/read_distribution_barcode;
    '''
    
    # generate the count per barcode
    tn5_barcode_count, rt_barcode_count, rt_to_tn5_conversion = bam_barcode_count_with_rna(args.input, args.tn5_barcodes, args.ligation_barcodes, args.rt_barcodes)
    # tn5_barcode_dic, rt_barcode_dic, rt_to_tn5_conversion
    tn5_barcodes = np.array(list(tn5_barcode_count.keys()))
    rt_barcodes = np.array(list(rt_barcode_count.keys()))
    
    tn5_count = np.array(list(tn5_barcode_count.values()))
    rt_count = np.array(list(rt_barcode_count.values()))


    # plot the barcode reads distribution and save the result to the ouput folder
    plot_name = (args.input.split('/')[-1]).split('.')[0]
    fig = plt.figure()
    plt.hist(tn5_count+rt_count, bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Number of unique reads')
    plt.axvline(x=args.cutoff, color='r', linestyle='dashed', linewidth=1)
    plt.yscale('log')
    fig_output = args.output_dir + '/' + plot_name + '.png'
    
    fig.savefig(fig_output)

    #also output the barcode number and distribution to the output folder
    with open(args.output_dir + '/' + plot_name + '_tn5.txt', 'w') as read_dist:
        for barcode in tn5_barcode_count:
            line = barcode + ', %d\n' %(tn5_barcode_count[barcode])
            read_dist.write(line)

    with open(args.output_dir + '/' + plot_name + '_rt.txt', 'w') as read_dist:
        for barcode in rt_barcode_count:
            line = barcode + ', %d\n' %(rt_barcode_count[barcode])
            read_dist.write(line)
    
    #Generate the read distribution in the split_bam/read_distribution_barcode
    
    #filter the barcode based on the cut_off value
    count_pass = tn5_count+rt_count >= int(args.cutoff)
    barcode_filtered = np.append(tn5_barcodes[count_pass],rt_barcodes[count_pass])
    tn5_barcode_filtered = tn5_barcodes[count_pass]
    # for barcode in barcode_count:
    #     if (barcode_count[barcode] >= int(args.cutoff)):
    #         barcode_filtered.append(barcode)


    #generate the output sam file and sample_list file based on DNA barcodes
    with open(args.output_dir + '/' + plot_name + '.' + 'sample_list.txt', 'w') as sample_list_file:
        output_files = {}
        paired_reads = {}
        for barcode in tn5_barcode_filtered:
            output_file = args.output_dir + '/' + plot_name + '.' + barcode + '.bam'
            output_files[barcode] = open(output_file, 'w')
            sample_list_file.write(plot_name + '.' + barcode + '\n')
    
    # output each read to the output sam file
    with pysam.AlignmentFile(args.input, "rb") as input_sam:
        for barcode in tn5_barcode_filtered:
            outname=output_files[barcode]
            paired_reads[outname] = pysam.AlignmentFile(outname, "wb", template=input_sam)
        
        for read in input_sam.fetch(): # until_eof=True to include unaligned reads
            # determine if read is RNA or DNA
            read_name_parts = read.query_name.split(',')
            if len(read_name_parts) != 6:
                raise ValueError('Number of barcodes in read name does not match expectation.')
            else:
                bc1 = read_name_parts[0]
                bc2 = read_name_parts[1]
                umi_rt = read_name_parts[2]
                # bc3 = read_name_parts[3]
                umi_rt = read_name_parts[2]
                if umi_rt == 'GGGGGG':
                    # this is DNA molecule
                    barcode = bc1+bc2
                    if barcode in barcode_filtered:
                        outname = output_files[barcode]
                        paired_reads[outname].write(read)
                else:
                    # this is RNA molecule
                    barcode = bc1+bc2
                    if barcode in barcode_filtered:
                        outname = output_files[rt_to_tn5_conversion[barcode]]
                        paired_reads[outname].write(read)

            # barcode = line.query_name.split(',')[0]+line.query_name.split(',')[1]
            # if barcode in barcode_filtered:
            #     outname=output_files[barcode]
            #     paired_reads[outname].write(line)

    #close the files:
    for barcode in tn5_barcode_filtered:
        output_files[barcode].close()
        


def fitler_bam_with_rna(args):
    '''
    Only filter barcodes based on cutoff and write to a single BAM file.
    '''
    
    # generate the count per barcode
    tn5_barcode_count, rt_barcode_count, rt_to_tn5_conversion = bam_barcode_count_with_rna(args.input, args.tn5_barcodes, args.ligation_barcodes, args.rt_barcodes)
    # tn5_barcode_dic, rt_barcode_dic, rt_to_tn5_conversion
    tn5_barcodes = np.array(list(tn5_barcode_count.keys()))
    rt_barcodes = np.array(list(rt_barcode_count.keys()))
    
    tn5_count = np.array(list(tn5_barcode_count.values()))
    rt_count = np.array(list(rt_barcode_count.values()))


    # plot the barcode reads distribution and save the result to the ouput folder
    plot_name = (args.input.split('/')[-1]).split('.')[0]
    fig = plt.figure()
    plt.hist(tn5_count+rt_count, bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Number of unique reads')
    plt.axvline(x=args.cutoff, color='r', linestyle='dashed', linewidth=1)
    plt.yscale('log')
    fig_output = args.output_dir + '/' + plot_name + '.png'
    
    fig.savefig(fig_output)

    #also output the barcode number and distribution to the output folder
    with open(args.output_dir + '/' + plot_name + '_tn5.txt', 'w') as read_dist:
        for barcode in tn5_barcode_count:
            line = barcode + ', %d\n' %(tn5_barcode_count[barcode])
            read_dist.write(line)

    with open(args.output_dir + '/' + plot_name + '_rt.txt', 'w') as read_dist:
        for barcode in rt_barcode_count:
            line = barcode + ', %d\n' %(rt_barcode_count[barcode])
            read_dist.write(line)
    
    #Generate the read distribution in the split_bam/read_distribution_barcode
    
    #filter the barcode based on the cut_off value
    count_pass = tn5_count+rt_count >= int(args.cutoff)
    barcode_filtered = np.append(tn5_barcodes[count_pass],rt_barcodes[count_pass])
    tn5_barcode_filtered = tn5_barcodes[count_pass]
    # for barcode in barcode_count:
    #     if (barcode_count[barcode] >= int(args.cutoff)):
    #         barcode_filtered.append(barcode)


    #generate the output sam file and sample_list file based on DNA barcodes
    with open(args.output_dir + '/' + plot_name + '.' + 'sample_list.txt', 'w') as sample_list_file:
        output_files = {}
        for barcode in tn5_barcode_filtered:
            output_file = args.output_dir + '/' + plot_name + '.' + barcode + '.bam'
            output_files[barcode] = open(output_file, 'w')
            sample_list_file.write(plot_name + '.' + barcode + '\n')
    

    out_path = os.path.join(args.output_dir, os.path.splitext(os.path.basename(args.input))[0] + '.bardcode_filt.bam')

    reads_written = 0
    reads_skipped = 0
    # output each read to the output sam file
    with pysam.AlignmentFile(args.input, "rb") as input_sam, \
         pysam.AlignmentFile(out_path, "wb", template=input_sam) as out_bam:
        
        for read in input_sam.fetch(): # until_eof=True to include unaligned reads
            # determine if read is RNA or DNA
            read_name_parts = read.query_name.split(',')
            if len(read_name_parts) != 6:
                raise ValueError('Number of barcodes in read name does not match expectation.')
            else:
                bc1 = read_name_parts[0]
                bc2 = read_name_parts[1]

                barcode = bc1+bc2
                if barcode in barcode_filtered:
                    out_bam.write(read)
                    reads_written += 1
                else:
                    reads_skipped += 1

    print(f'Reads written out: {reads_written}')
    print(f'Reads with low barcode counts filtered: {reads_skipped}')
 






def bam_barcode_count_with_sss(bam_file, barcode_tn5_file, barcode_ligation_file, barcode_sss_file=None):
    
    #generate the barcode list and barcode dictionary
    barcode_tn5 = []
    with open(barcode_tn5_file) as barcodes:
        for barcode in barcodes:
            barcode_tn5.append(barcode.strip())
    
    barcode_ligation = []
    with open(barcode_ligation_file) as barcodes:
        for barcode in barcodes:
            barcode_ligation.append(barcode.strip())
    
    barcode_combination = []
    barcode_dic = {}
    for tn5 in barcode_tn5:
        for ligation in barcode_ligation:
            barcode = tn5 + ligation
            barcode_combination.append(barcode)
            barcode_dic[barcode] = 0

    if barcode_sss_file != None:
        barcode_sss_dic = defaultdict()
        with open(barcode_sss_file) as barcodes:
            for barcode in barcodes:
                barcode_sss_dic[barcode.strip()] = barcode_dic.copy()

    # check for BAM index
    with pysam.AlignmentFile(bam_file, "rb") as input_sam:
        try:
            input_sam.check_index()
        except ValueError:
            pysam.index(bam_file)

    skipped = []
    #read the sam file, and count the number per barcode
    with pysam.AlignmentFile(bam_file, "rb") as input_sam:
        if barcode_sss_file != None:
            for line in input_sam.fetch():
                sss_barcode = line.query_name.split(',')[2]
                barcode = line.query_name.split(',')[0]+line.query_name.split(',')[1]
                try:
                    barcode_sss_dic[sss_barcode][barcode] += 1
                except:
                    skipped.append(sss_barcode)    
                     
    unq_skipped = set(skipped)
    print(f'SSS barcodes skipped: {unq_skipped}')

    if barcode_sss_file != None:
        return barcode_sss_dic
    

def split_bam_with_sss(args):
    '''
    If BAM is not split by SSS, include SSS in barcode for splitting
    '''
    
    # generate the count per barcode
    barcode_count = bam_barcode_count_with_sss(args.input, args.tn5_barcodes, args.ligation_barcodes, args.sss_barcodes)
    
    # plot the barcode reads distribution and save the result to the ouput folder
    plot_name = (args.input.split('/')[-1]).split('.')[0]
    fig = plt.figure()
    if args.sss_barcodes != None:
        tmp = [list(barcode_count[k].values()) for k in barcode_count.keys()]
        barcode_count_all = [sum(x) for x in zip(*tmp)]
        plt.hist(barcode_count_all, bins=100)
    else:
        plt.hist(barcode_count.values(), bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Number of unique reads')
    plt.axvline(x=args.cutoff, color='r', linestyle='dashed', linewidth=1)
    plt.yscale('log')
    fig_output = args.output_dir + '/' + plot_name + '.png'
    
    fig.savefig(fig_output)

    #Generate the read distribution in the split_bam/read_distribution_barcode

    #also output the barcode number and distribution to the output folder
    with open(args.output_dir + '/' + plot_name + '.txt', 'w') as read_dist:
        for sss in barcode_count:
            for barcode in barcode_count[sss]:
                line = sss + '.' + barcode + ', %d\n' %(barcode_count[sss][barcode])
                read_dist.write(line)

    #filter the barcode based on the cut_off value
    barcode_filtered = []
    for sss in barcode_count:
        for barcode in barcode_count[sss]:
            if (barcode_count[sss][barcode] >= int(args.cutoff)):
                barcode_filtered.append(sss + '.' + barcode)

    #generate the output sam file and sample_list file
    with open(args.output_dir + '/' + plot_name + '.' + 'sample_list.txt', 'w') as sample_list_file:
        output_files = {}
        paired_reads = {}
        for barcode in barcode_filtered:
            output_file = args.output_dir + '/' + barcode + '.bam'
            output_files[barcode] = open(output_file, 'w')
            sample_list_file.write(barcode + '\n')
    

    # output each read to the output sam file
    with pysam.AlignmentFile(args.input, "rb") as input_sam:
        for barcode in barcode_filtered:
            outname=output_files[barcode]
            paired_reads[outname] = pysam.AlignmentFile(outname, "wb", template=input_sam)
        
        for line in input_sam.fetch(): # until_eof=True to include unaligned reads
            barcode = line.query_name.split(',')[2] + '.' + line.query_name.split(',')[0] + line.query_name.split(',')[1]
            if barcode in barcode_filtered:
                outname=output_files[barcode]
                paired_reads[outname].write(line)

    #close the files:
    for barcode in barcode_filtered:
        output_files[barcode].close()

  

def filter_bam_with_sss(args):
    '''
    If BAM is not split by SSS, include SSS in barcode for filtering barcodes based on cutoff and write to a single BAM file.
    '''
    # generate the count per barcode
    barcode_count = bam_barcode_count_with_sss(args.input, args.tn5_barcodes, args.ligation_barcodes, args.sss_barcodes)
    
    # plot the barcode reads distribution and save the result to the ouput folder
    plot_name = (args.input.split('/')[-1]).split('.')[0]
    fig = plt.figure()

    if args.sss_barcodes != None:
        tmp = [list(barcode_count[k].values()) for k in barcode_count.keys()]
        barcode_count_all = [sum(x) for x in zip(*tmp)]
        plt.hist(barcode_count_all, bins=100)
    else:
        plt.hist(barcode_count.values(), bins=100)
    plt.ylabel('Frequency')
    plt.xlabel('Number of unique reads')
    plt.axvline(x=args.cutoff, color='r', linestyle='dashed', linewidth=1)
    plt.yscale('log')
    fig_output = args.output_dir + '/' + plot_name + '.png'
    
    fig.savefig(fig_output)

    #also output the barcode number and distribution to the output folder
    with open(args.output_dir + '/' + plot_name + '.txt', 'w') as read_dist:
        for sss in barcode_count:
            for barcode in barcode_count[sss]:
                line = sss + '.' + barcode + ', %d\n' %(barcode_count[sss][barcode])
                read_dist.write(line)
    
    #Generate the read distribution in the split_bam/read_distribution_barcode
    #filter the barcode based on the cut_off value
    barcode_filtered = []
    for sss in barcode_count:
        for barcode in barcode_count[sss]:
            if (barcode_count[sss][barcode] >= int(args.cutoff)):
                barcode_filtered.append(sss + '.' + barcode)

    out_path = os.path.join(args.output_dir, os.path.splitext(os.path.basename(args.input))[0] + '.bardcode_filt.bam')

    reads_written = 0
    reads_skipped = 0
    with pysam.AlignmentFile(args.input, "rb") as input_bam, \
         pysam.AlignmentFile(out_path, "wb", template=input_bam) as out_bam:
        for line in input_bam.fetch():
            barcode = line.query_name.split(',')[2] + '.' + line.query_name.split(',')[0] + line.query_name.split(',')[1]
            if barcode in barcode_filtered:
                out_bam.write(line)
                reads_written += 1
            else:
                reads_skipped += 1

    print(f'Reads written out: {reads_written}')
    print(f'Reads with low barcode counts filtered: {reads_skipped}')
    


if __name__ == "__main__":
    main()
    