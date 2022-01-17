import pyfastx
import argparse
from Levenshtein import distance
import os
import gzip

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(description ='This script accepts and opens' +
    'trimmed R1 and ordered R2 files split according to sss barcodes, an reformatted' +
    'info file from cutadapt, a bridge sequence, a tn5 barcode list and a ligation' +
    'barcode list; then extracts and matches the nearest tn5 and ligation barcodes' +
    'and attach barcodes to R1 and R2 read names. The output files include new' +
    'R1 and R2 files, noME.R1.fq.gz')

    parser.add_argument('-r1', '--read_1', action = 'store', metavar = 'FILE',
                        help = 'Read 1 input fastq file, output name will be '+
                        'derived from first part of _ split input name ')
    parser.add_argument('-r2', '--read_2', action = 'store', metavar = 'FILE',
                        help = 'Read 2 input fastq file')
    parser.add_argument('-bc', '--bc1_bc2', action = 'store', metavar = 'FILE',
                        help = 'Barcode 1 and 2 extracted from cutadapt info file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'DIR',
                        help = 'Output folder')
    parser.add_argument('-bt', '--barcode_tn5', action = 'store', metavar = 'FILE',
                        help = 'File containing Tn5 barcodes')
    parser.add_argument('-bl', '--barcode_lig', action = 'store', metavar = 'FILE',
                    help = 'File containing ligation barcodes')
    parser.add_argument('-tn5', '--mismatch_tn5', metavar = 'INT', type = int, default = 0,
                        help = ('Number of mismatches allowed identifying Tn5 barcode' +
                                'Default: 1'))
    parser.add_argument('-l', '--mismatch_lig', metavar = 'INT', type = int, default = 3,
                        help = ('Number of mismatches allowed identifying ligation barcode' +
                                'Default: 1'))

    args = parser.parse_args(args)

    return args


def main():

    args = parse_arguments()

    # args = parse_arguments('-r1 /mnt/data/sci-l3/test/trimmed/yi292/yi292_CACGTG.R1.trimmed.fq.gz ' \
    #                     '-r2 /mnt/data/sci-l3/test/trimmed/yi292/yi292_CACGTG.R2.trimmed.fq.gz ' \
    #                     '-o /mnt/data/sci-l3/test3 ' \
    #                     '-bc /mnt/data/sci-l3/test/trimmed/yi292/yi292_CACGTG.R1.bc1.bc2.txt.gz ' \
    #                     '-bt /mnt/data/Projects/sciStrand-seq/barcodes/barcode_tn5.txt ' \
    #                     '-bl /mnt/data/Projects/sciStrand-seq/barcodes/barcode_ligation_104.txt ' \
    #                     '-tn5 1 ' \
    #                     '-l 1'.split())

    tn5_bar = tn5_barcodes(args)
    lig_bar = ligation_barcodes(args)

    tn5_ligation_barcode_attach(args.read_1, args.read_2, args.bc1_bc2, args.output, 
                                tn5_bar, lig_bar, args.mismatch_tn5, 
                                args.mismatch_lig)


def tn5_ligation_barcode_attach(r1, r2, bc1_bc2, output_folder, tn5_barcodes, 
                                ligation_barcodes, mismatch_tn5, mismatch_lig):
    '''
    Match tn5 and ligation barcodes and attach them to fastq name
    '''
    sample = os.path.basename(r1).split('.')[0]

    # Gzip after, much faster
    fq1_out = output_folder + '/' + sample + '.R1.trimmed.attached.fastq'
    fq2_out = output_folder + '/' + sample + '.R2.trimmed.attached.fastq'
    noME_out = output_folder + '/' + sample + '.noME.R1-2.fastq'

    no_tn5_lig_barcode = 0
    written_out = 0
    total_reads = 0

    # Gzip after output
    with open(fq1_out, 'wt') as R1_out, \
         open(fq2_out, 'wt') as R2_out, \
         open(noME_out, 'wt') as nome_out, \
         gzip.open(bc1_bc2, 'rt') as bc12:

        # parse input fastq files
        for s1, s2, bc in zip(pyfastx.Fastq(r1, build_index=False), 
                        pyfastx.Fastq(r2, build_index=False),
                        bc12.readlines()):

            name1, seq1, qual1 = s1  
            name2, seq2, qual2 = s2
            
            assert name1.split(' ')[0] == name2.split(' ')[0], 'Read names do not match!'
           
            total_reads += 1
        
            if len(bc.split()) < 2:
                target_tn5, target_ligation = "NA", "NA"
            else:
                target_ligation, target_tn5 = bc.split()
            
            found = False
            if target_tn5 != "NA":
                for tn5 in tn5_barcodes:
                    mm_tn5 = distance(tn5, target_tn5)
                    if (mm_tn5 <= mismatch_tn5):
                        for ligation in ligation_barcodes:
                            mm_ligation = distance(ligation, target_ligation)
                            if(mm_ligation <= mismatch_lig):
                                found = True
                                # write out read valid barcodes
                                name1_new = '@' + tn5 + ',' + ligation + ',' + name1[1:]
                                name2_new = '@' + tn5 + ',' + ligation + ',' + name2[1:]

                                R1_out.write(name1_new + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                                R2_out.write(name2_new + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')
                                # stop seaching
                                written_out += 1
                                break
                    if found:
                        break
                if not found:
                    no_tn5_lig_barcode += 1
            else:
                nome_out.write('@' + name1 + '\n' + seq1 + '\n' + '+' + '\n' + qual1 + '\n')
                nome_out.write('@' + name2 + '\n' + seq2 + '\n' + '+' + '\n' + qual2 + '\n')

    print(f'Total read pairs: {total_reads}')
    print(f'Total read pairs written: {written_out}')
    print(f'Read pairs without Tn5 or ligation barcode: {no_tn5_lig_barcode}')


def tn5_barcodes(args):
    '''Read tn5 barcodes
    '''
    barcode_tn5 = []
    with open(args.barcode_tn5) as barcodes:
        for barcode in barcodes:
            barcode_tn5.append(barcode.strip())
    return barcode_tn5

def ligation_barcodes(args):
    '''Read ligation barcodes
    '''
    barcode_ligation = []
    with open(args.barcode_lig) as barcodes:
        for barcode in barcodes:
            barcode_ligation.append(barcode.strip())
    return barcode_ligation

if __name__ == '__main__':
    main()