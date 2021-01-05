import subprocess
import sys
import pysam
import gzip
import matplotlib as mpl
mpl.use('Agg') # non-interactive backend required for cluster nodes
import matplotlib.pyplot as plt
import argparse
import os


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(description =
            "This script accept a bam file, an output folder," +
            "two barcode files and a cut_off value, then it will call the" +
            "bam_barcode_count function and get the total read count per barcode," +
            "then it use the cut_off value to filter the barcode," +
            "and generate the output samfile for single cells," +
            "generate the reads distribution in the split_bam/read_distribution_barcode")

    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output_dir', action = 'store', metavar = 'DIR',
                        help = 'Output folder')
    parser.add_argument('-tb', '--tn5_barcodes', action = 'store', metavar = 'FILE',
                        help = 'File with Tn5 barcodes')
    parser.add_argument('-lb', '--ligation_barcodes', action = 'store', metavar = 'FILE',
                        help = 'File with ligation barcodes')
    parser.add_argument('--cutoff', metavar = 'INT', type = int, default = 10000,
                        help = ('Number of unique reads cutoff for splitting single cells' +
                                'Default: 10000'))

    args = parser.parse_args(args)

    return args

def main():

    args = parse_arguments()

    # create directory tree if it doesn't exist
    os.makedirs(args.output_dir, exist_ok=True)

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
        outname=output_files[barcode]
        paired_reads[outname].close()



if __name__ == "__main__":
    main()
    