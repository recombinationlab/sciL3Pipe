import pyfastx
import argparse
from Levenshtein import distance
from collections import defaultdict
import os

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(description ='This script accepts and opens' +
    'a R1 file, a R2 file, an RT primer sequence, and a SSS barcode list;' +
    'then reorders the read pairs, and extracts and matches the nearest sss' +
    'barcode and attach the sss barcode to R1 and R2 read names.' +
    'the output files include ordered R1 and R2 files,' +
    'noRTprimer.R1.fq.gz and noRTprimer.R2.fq.gz')

    parser.add_argument('-r1', '--read_1', action = 'store', metavar = 'FILE',
                        help = 'Read 1 input fastq file, output name will be '+
                        'derived from first part of _ split input name ')
    parser.add_argument('-r2', '--read_2', action = 'store', metavar = 'FILE',
                        help = 'Read 2 input fastq file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'DIR',
                        help = 'Output folder')
    parser.add_argument('-b', '--barcode_sss', action = 'store', metavar = 'FILE',
                        help = 'File containing SSS barcodes')
    parser.add_argument('-sss', '--mismatch_sss', metavar = 'INT', type = int, default = 0,
                        help = ('Number of mismatches allowed identifying SSS barcode' +
                                'Default: 0'))
    parser.add_argument('-rt', '--mismatch_rt', metavar = 'INT', type = int, default = 3,
                        help = ('Number of mismatches allowed identifying RT barcode' +
                                'Default: 3'))
    parser.add_argument('-p', '--primer', action = 'store', metavar = 'STR',
                        help = 'RT primer sequence')

    args = parser.parse_args(args)

    return args

def main():

    args = parse_arguments()

    # args = parse_arguments('-r1 /mnt/data/sci-l3/raw_fastq/yi292_S6_R1_001.fastq.gz ' \
    #                        '-r2 /mnt/data/sci-l3/raw_fastq/yi292_S6_R2_001.fastq.gz ' \
    #                        '-o /mnt/data/sci-l3/test3 ' \
    #                        '-b /mnt/data/Projects/sciStrand-seq/barcodes/barcode_sss_lib292.txt ' \
    #                        '-sss 0 ' \
    #                        '-rt 3 ' \
    #                        '-p GGGATGCAGCTCGCTCCTG'.split())


    # create directory tree if it doesn't exist
    os.makedirs(args.output, exist_ok=True)

    barcode_dict = parse_barcodes(args)
    
    SSS_barcode_extract(args.read_1, args.read_2, args.output, barcode_dict,
                        args.primer, args.mismatch_sss, args.mismatch_rt)


def SSS_barcode_extract(r1, r2, output_folder, barcode_dict, RT_primer, 
                        mismatch_sss, mismatch_RT):
    '''
    '''

    sample = os.path.basename(r1).split('.')[0]
    
    # output files
    # R1_out_file = output_folder + "/" + sample + ".R1.ordered.fastq"
    # R2_out_file = output_folder + "/" + sample + ".R2.ordered.fastq"
    noRT_R1_out_file = output_folder + "/" + sample + ".noRTprimer.R1.fastq"
    noRT_R2_out_file = output_folder + "/" + sample + ".noRTprimer.R2.fastq"
  
    total_reads = 0
    sss_not_found = 0
    junk_count = 0
    written_out = 0
    swap_written_out = 0

    r1_out_files = defaultdict()
    r2_out_files = defaultdict()

    out_file_counts = defaultdict(int)

    # Gzip after, much faster
    # open(R1_out_file, 'w') as R1_out, \
    # open(R2_out_file, 'w') as R2_out, \
    with open(noRT_R1_out_file, 'w') as noRT_R1, \
         open(noRT_R2_out_file, 'w') as noRT_R2:

        # parse input fastq files
        for s1, s2 in zip(pyfastx.Fastq(r1, build_index=False), 
                        pyfastx.Fastq(r2, build_index=False)):

            name1, seq1, qual1 = s1  
            name2, seq2, qual2 = s2

            assert name1.split(' ')[0] == name2.split(' ')[0], 'Read names do not match!'
            assert len(seq1) > 29 and len(seq2) > 29, 'Sequence is shorter than 29 bps'

            target_rt_r1 = seq1[9:29]
            target_rt_r2 = seq2[9:29]
            target_sss_r1 = seq1[4:10]
            target_sss_r2 = seq2[4:10]
            total_reads += 1
            
            dist_rt_r1 = distance(RT_primer, target_rt_r1)
            dist_rt_r2 = distance(RT_primer, target_rt_r2)

            # if RT primer is in read 1, write out
            if dist_rt_r1 <= mismatch_RT:
                # find barcode
                if mismatch_sss == 0:
                    # if barcode present, will return True
                    barcode = target_sss_r1 if barcode_dict.get(target_sss_r1, False) else False
                else: # if mismatch > 0, calculate Levenshtein distance
                    barcode = (b if distance(b, target_sss_r1) <= mismatch_sss else False for b in barcode_dict.keys())
              
                if barcode:
                    UMI = seq1[:4]
                    name1_new = '@' + barcode + ',' + UMI + ',' + name1[1:]
                    name2_new = '@' + barcode + ',' + UMI + ',' + name2[1:]
                    # split by SSS
                    try:
                        r1_out_files[barcode].write(name1_new + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                    except KeyError:
                        r1_out_files[barcode] = open(output_folder + '/' + sample + '_' + barcode + '.R1.ordered.fastq', 'w')
                        r1_out_files[barcode].write(name1_new + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                    try:
                        r2_out_files[barcode].write(name2_new + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
                    except KeyError:
                        r2_out_files[barcode] = open(output_folder + '/' + sample + '_' + barcode + '.R2.ordered.fastq', 'w')
                        r2_out_files[barcode].write(name2_new + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
                    # R1_out.write(name1_new + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                    # R2_out.write(name2_new + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
                    out_file_counts[barcode] += 1
                    written_out += 1
                else:
                    sss_not_found += 1

            # if RT primer found in R2, swap R1 R2 but keep read name
            elif dist_rt_r2 <= mismatch_RT:
                # find barcode
                if mismatch_sss == 0:
                    # if barcode present, will return True
                    barcode = target_sss_r2 if barcode_dict.get(target_sss_r2, False) else False
                else: # if mismatch > 0, calculate Levenshtein distance
                    barcode = (b if distance(b, target_sss_r2) <= mismatch_sss else False for b in barcode_dict.keys())
              
                if barcode:
                    UMI = seq1[:4]
                    name1_new = '@' + barcode + ',' + UMI + ',swap,' + name1[1:]
                    name2_new = '@' + barcode + ',' + UMI + ',swap,' + name2[1:]
                    # swap
                    # split by SSS
                    try:
                        r1_out_files[barcode].write(name2_new + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
                    except KeyError:
                        r1_out_files[barcode] = open(output_folder + '/' + sample + '_' + barcode + '.R1.ordered.fastq', 'w')
                        r1_out_files[barcode].write(name2_new + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
                    try:
                        r2_out_files[barcode].write(name1_new + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                    except KeyError:
                        r2_out_files[barcode] = open(output_folder + '/' + sample + '_' + barcode + '.R2.ordered.fastq', 'w')
                        r2_out_files[barcode].write(name1_new + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                    # R1_out.write(name2_new + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
                    # R2_out.write(name1_new + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                    out_file_counts[barcode] += 1
                    swap_written_out += 1
                else:
                    sss_not_found += 1

            # Write out reads missing the RT primer
            else:
                noRT_R1.write(name1 + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                noRT_R2.write(name2 + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
                junk_count += 1

    # close all out files
    for f in r1_out_files.values():
        f.close()
    for f in r2_out_files.values():
        f.close()

    print('Total read pairs processed:', total_reads)
    print('Read pairs without SSS barcode:', sss_not_found)
    print('Read pairs without RT primer:', junk_count)
    print('R1-R2 orientation read pairs written out:', written_out)
    print('R2-R1 orientation read pairs written out:', swap_written_out)

    for k, v in out_file_counts.items():
        print(f'SSS with barcode {k} read pair count: {v}')

def parse_barcodes(args):
    '''
    Generate the barcode dict
    If mismatched for SSS barcode = 0, then use dict to find correct barcode
    '''
    barcode_dict = defaultdict()
    with open(args.barcode_sss, 'r') as barcodes:
        for barcode in barcodes:
            barcode_dict[barcode.strip()] = True
    
    return barcode_dict


if __name__ == '__main__':
    main()