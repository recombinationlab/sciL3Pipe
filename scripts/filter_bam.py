
import pysam
import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline
from contextlib import contextmanager


def parse_arguments(args=None):

    parser = argparse.ArgumentParser(description =
            "This program removes reads from a BAM file according to the " +
            "filtering criteria. \n" +
            "1. Removes unmapped reads \n" +
            "2. Removes reads with invalid edit distance \n" +
            "3. Removes reads with invalid MAPQ scores \n" +
            "4. Removes reads with invalid insert size (Paired-end) \n" +
            "5. Removes reads on different chromosomes (Paired-end) \n" +
            "6. Removes reads in orientation other than FR (Paired-end) \n" +
            "The script requires that the reads are name sorted, not coordinate sorted")

    parser.add_argument('-i', '--input', action = 'store', metavar = 'FILE',
                        help = 'Input BAM file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Output BAM file containing reads that pass')
    parser.add_argument('--edit_max', action = 'store', metavar = 'Y',
                        type = int, default = 4,
                        help = ('If the edit distance between the read ' +
                                'sequence and the reference sequence falls ' +
                                'within [X, Y], keep the read. Otherwise ' +
                                'remove it. To ignore filter set both min and max ' +
                                'to 0. Default = 4'))
    parser.add_argument('--edit_min', action = 'store', metavar = 'X',
                        type = int, default = 0,
                        help = ('If the edit distance between the read ' +
                                'sequence and the reference sequence falls ' +
                                'within [X, Y], keep the read. Otherwise ' +
                                'remove it. Default = 0'))
    parser.add_argument('--insert_max', action = 'store', metavar = 'Y',
                        type = int, default = 2000,
                        help = ('If the insert size ' +
                                'falls within [X, Y], keep the read. Otherwise ' +
                                'remove it. To ignore filter set both min and max ' +
                                'to 0. Default = 2000'))
    parser.add_argument('--insert_min', action = 'store', metavar = 'X',
                        type = int, default = 0,
                        help = ('If the inset size' +
                                'falls within [X, Y], keep the read. Otherwise ' +
                                'remove it. Default = 0'))
    parser.add_argument('--mapq_min', action = 'store', metavar = 'M',
                        type = int, default = 0,
                        help = ('If the MAPQ score of this read falls within ' +
                                '[M, N], keep the read. Otherwise remove ' +
                                'it. Default = 0'))
    parser.add_argument('--mapq_max', action = 'store', metavar = 'N',
                        type = int, default = 255,
                        help = ('If the MAPQ score of this read falls within ' +
                                '[M, N], keep the read. Otherwise remove ' +
                                'it. Default = 255'))
    parser.add_argument('--paired', action = 'store_true',
                        help = ('Set if the BAM file contains a paired-end ' +
                                'alignment. Defaults to single-read alignment.'))
    parser.add_argument('--ignore_orientation', action = 'store_true',
                        help = ('Keep paired end reads not in FR or RF orientation'))
    parser.add_argument('--keep_trans', action = 'store_true',
                        help = ('Keep paired end reads not on same chromosome'))
    parser.add_argument('--hybrid_reference', action = 'store_true',
                        help = ('Specifies if alignment was performed with a hybrid' +
                                'reference (Ref 1 = 1, 2, 3, etc.; Ref 2 = chr1, chr2, chr3, etc.'))

    args = parser.parse_args(args)

    if args.edit_max < args.edit_min:
        exit("Invalid edit-distance parameters. Max: " +
             str(args.edit_max) + ", min: " + str(args.edit_min))
    if args.edit_max == args.edit_min:
        print("No edit distance filter")

    if args.mapq_max < args.mapq_min:
        exit("Invalid MAPQ parameters. Max: " +
             srt(args.mapq_max) + ", min: " + str(args.mapq_min))

    if args.insert_max < args.insert_min:
        exit("Invalid insert size parameters. Max: " +
             str(args.insert_max) + ", min: " + str(args.insert_min))
    elif args.insert_max == args.insert_min:
        print('No insert size filter')

    if any([args.ignore_orientation, args.keep_trans, args.hybrid_reference]) and not args.paired:
        print("Paired-end options ignored for single-end alignment")

    return args


def main():
    args = parse_arguments()

    # # test args
    # args = parse_arguments('--edit_max 0 --edit_min 0 --hybrid_reference ' \
    #                          '-i /mnt/data/Projects/sci_lianti/yi293_TTCGAG.PE.bwa.hg38.collate.bam ' \
    #                          '--insert_max 1000 --insert_min 0 --mapq_max 255 --mapq_min 0 ' \
    #                          '-o /mnt/data/Projects/sci_lianti/yi293_TTCGAG.PE.bwa.hg38.filt.bam ' \
    #                          '--paired'.split())

    # mapq = get_mapq_scores(args)
    # ed = get_edit_distance(args)
    # insert_size = get_insert_size(args)

    if args.paired:
        filter_paired_reads(args)
    else:
        filter_single_reads(args)


def plot_insert_size(insert_size, is_min, is_max, save_plt):
    '''
    Plot stacked bar plot with insert size with other filters overlaid

    Args:
        insert_size(np.array): numpy array of values to plot
        save_plt(str): output path name to save plot
    '''

    if type(insert_size) == list:
        out_png = save_plt + '.insert_size_pass.png'
        plt.hist(insert_size, bins=100, density=False, histtype='bar')
    else:
        out_png = save_plt + '.insert_size.png'
        plt.hist(insert_size, bins=100, density=False, histtype='bar', stacked=True, 
                    label=['Pass', 'Fail'])
        plt.legend(prop={'size': 10})
    plt.ylabel('Read count')
    plt.xlabel('Insert size')
    # if both min and max == 0, insert filter ignored
    if is_min > 0:
        plt.axvline(x=is_min, color='r', linestyle='dashed', linewidth=1)
    if is_max > 0:
        plt.axvline(x=is_max, color='r', linestyle='dashed', linewidth=1)
    plt.yscale('log')
    plt.savefig(out_png, bbox_inches='tight')
    plt.close()

def plot_mapq(mapq_scores, min_mapq, max_mapq, save_plt):
    '''
    Plot counts histogram of BAM MapQ scores
    with reads that would be filtered out based on other criteria 

    Args:
        mapq_scores(list): list of MapQ scores to plot
        min_mapq(int): Cutoff for minimum allowed MapQ score
        max_mapq(int): Cutoff for Maximum allowed MapQ score
        save_plt(str): output path name to save plot
    '''
    out_png = save_plt + '.MAPQ.png'

    plt.hist(mapq_scores, bins=max(max(mapq_scores[0]), max(mapq_scores[1])), 
    density=False, histtype='bar', 
    stacked=True, label=['Pass', 'Fail'])
    plt.legend(prop={'size': 10})
    plt.ylabel('Read count')
    plt.xlabel('MapQ score')
    plt.axvline(x=min_mapq, color='r', linestyle='dashed', linewidth=1)
    if max_mapq < 255:
        plt.axvline(x=max_mapq, color='r', linestyle='dashed', linewidth=1)
    plt.savefig(out_png, bbox_inches='tight')
    plt.close()

def plot_mismatches(edit_distance, min_nm, max_nm, save_plt):
    '''
    Plot NM flag, number of mismatches in a read

    Args:
        edit_distance(list): list of mismatches to plot
        min_nm(int): minimum allowed edit distance
        max_nm(int): maximum allowed edit distance
        save_plt(str): output path name to save plot
    '''

    out_png = save_plt + '.NM.png'
    plt.hist(edit_distance, bins=20, range=(0,20), density=False, histtype='bar', 
    stacked=True, label=['Pass', 'Fail'])
    plt.legend(prop={'size': 10})
    plt.ylabel('Read count')
    plt.xlabel('Edit distance')
    if min_nm > 0:
        plt.axvline(x=min_nm, color='r', linestyle='dashed', linewidth=1)
    if max_nm > 0:
        plt.axvline(x=max_nm, color='r', linestyle='dashed', linewidth=1)
    plt.xlim(0, 20)
    plt.savefig(out_png, bbox_inches='tight')
    plt.close()


def filter_single_reads(args):

    valid_reads = 0
    invalid_reads = 0

    # pass fail based on other filters
    mapq_score_fail = []
    mapq_score_pass = []
    edit_distance_fail = []
    edit_distance_pass = []

    none_context = contextmanager(lambda: iter([None]))()

    with pysam.AlignmentFile(args.input, "rb") as input_bam, \
         (pysam.AlignmentFile(args.output, "wb", template = input_bam) if args.output else none_context) as output_bam:

        for read in input_bam.fetch(until_eof = True):
            valid_edit_distance = has_valid_edit_distance(read, args.edit_min, args.edit_max),
            valid_mapq_score = has_valid_mapq_score(read, args.mapq_min, args.mapq_max)

            # For MAPQ plotting
            if all([not read.is_unmapped,
                    valid_edit_distance]):
                mapq_score_pass.append(read.mapping_quality)
            elif not read.is_unmapped:
                mapq_score_fail.append(read.mapping_quality)

            # Edit distance plotting
            if all([not read.is_unmapped,
                    valid_mapq_score]):
                edit_distance_pass.append(int(read.get_tag('NM')))
            elif not read.is_unmapped:
                edit_distance_fail.append(int(read.get_tag('NM')))

            if all([not read.is_unmapped,
                    valid_edit_distance,
                    valid_mapq_score]):

                valid_reads += 1
                if args.output:
                    output_bam.write(read)
            else:
                invalid_reads += 1

    print("Valid reads written out:", valid_reads)
    print("Invalid reads filtered:", invalid_reads)

    plot_mapq(np.array([mapq_score_pass, mapq_score_fail], dtype=object), 
              args.mapq_min, args.mapq_max, args.input)

    plot_mismatches(np.array([edit_distance_pass, edit_distance_fail], dtype=object),
                    args.edit_min, args.edit_max, args.input)



def filter_paired_reads(args):
    '''
    Main filtering of paired end reads
    '''
   
    EOF = False
    valid_reads = 0
    invalid_reads = 0
    singletons = 0
    valid_count = 0
    trans_count = 0
    unmapped_reads = 0
    insert_size_filt = 0
    edit_distance_filt = 0
    mapq_score_filt = 0

    # pass fail based on other filters
    mapq_score_fail = []
    mapq_score_pass = []
    edit_distance_fail = []
    edit_distance_pass = []
    insert_size_fail = []
    insert_size_pass = []

    none_context = contextmanager(lambda: iter([None]))()

    with pysam.AlignmentFile(args.input, 'rb') as input_bam, \
        (pysam.AlignmentFile(args.output, 'wb', template = input_bam) if args.output else none_context) as output_bam:

        reads = input_bam.fetch(until_eof = True)
        current_read = next(reads)
        for read in reads:

            if read.query_name == current_read.query_name:

                # only check mapped reads
                if all([not read.is_unmapped,
                        not current_read.is_unmapped]):

                    if args.ignore_orientation:
                        valid_orientation = True
                    else:
                        valid_orientation = is_FR_RF(current_read, read)
                        if valid_orientation:
                            valid_count += 1
                        
                    no_trans_reads = has_same_chromosome(current_read, read, 
                                                            args.hybrid_reference,
                                                            args.keep_trans)
                    if not no_trans_reads:
                        trans_count += 1

                    valid_insert_size = has_valid_insert_size(read, args.insert_min, 
                                                            args.insert_max)
                    
                    if valid_insert_size:
                        insert_size_filt += 1

                    valid_edit_distance = all([has_valid_edit_distance(read, 
                                                                    args.edit_min, 
                                                                    args.edit_max),
                                            has_valid_edit_distance(current_read, 
                                                                    args.edit_min, 
                                                                    args.edit_max)])
                
                    if valid_edit_distance:
                        edit_distance_filt += 1

                    valid_mapq_score = all([has_valid_mapq_score(read, 
                                                                args.mapq_min, 
                                                                args.mapq_max),
                                            has_valid_mapq_score(current_read, 
                                                                args.mapq_min, 
                                                                args.mapq_max)])
                    if valid_mapq_score:
                        mapq_score_filt += 1
                # else:
                #     valid_edit_distance = False
                #     valid_insert_size = False
                #     valid_orientation = False
                #     no_trans_reads = False
                #     valid_mapq_score = False

                # For MAPQ plotting
                if all([not read.is_unmapped,
                        not current_read.is_unmapped,
                        valid_edit_distance,
                        valid_orientation,
                        no_trans_reads,
                        valid_insert_size]):
                    mapq_score_pass.append(min(read.mapping_quality, 
                                            current_read.mapping_quality))
                elif all([not read.is_unmapped,
                          not current_read.is_unmapped]):
                    mapq_score_fail.append(min(read.mapping_quality, 
                                            current_read.mapping_quality))

                # Edit distance plotting
                if all([not read.is_unmapped,
                        not current_read.is_unmapped,
                        valid_mapq_score,
                        valid_orientation,
                        no_trans_reads,
                        valid_insert_size]):
                    edit_distance_pass.append(max(int(read.get_tag('NM')), 
                                               int(current_read.get_tag('NM'))))
                elif all([not read.is_unmapped,
                          not current_read.is_unmapped]):
                    edit_distance_fail.append(max(int(read.get_tag('NM')), 
                                               int(current_read.get_tag('NM'))))

                # Insert size plotting
                if all([not read.is_unmapped,
                        not current_read.is_unmapped,
                        valid_edit_distance,
                        valid_mapq_score,
                        valid_orientation,
                        no_trans_reads]):
                    insert_size_pass.append(abs(read.template_length))
                elif all([not read.is_unmapped,
                          not current_read.is_unmapped]):
                    insert_size_fail.append(abs(read.template_length))


                # Filtering and BAM write out
                if all([not read.is_unmapped,
                        not current_read.is_unmapped,
                        valid_edit_distance,
                        valid_mapq_score,
                        valid_orientation,
                        no_trans_reads,
                        valid_insert_size]):
                        
                    valid_reads += 1 
                    if args.output:
                        output_bam.write(current_read)
                        output_bam.write(read)
                elif all([not read.is_unmapped,
                          not current_read.is_unmapped]):
                    unmapped_reads += 1
                else:
                    invalid_reads += 1

                # reset counter
                try:
                    current_read = next(reads)
                except:
                    StopIteration
                    EOF = True
            else:
                singletons += 1
                current_read = read
    if not EOF:
        singletons += 1

    print("Valid read pairs written out:", valid_reads)
    print("Invalid read pairs filtered:", invalid_reads)
    print("Singleton reads filtered:", singletons)
    print("Unmapped reads:", unmapped_reads)
    if args.keep_trans and args.hybrid_reference:
        print("Hybrid read pairs count:", trans_count)
    elif not args.keep_trans and args.hybrid_reference:
        print("Trans or hybrid read pairs count:", trans_count)
    elif not args.keep_trans and not args.hybrid_reference:
        print("Trans read pairs count:", trans_count)

    if not args.ignore_orientation:
        print("Proper orientation read pairs:", valid_count)

    print("Valid insert size read pairs:", insert_size_filt)
    print("Valid edit distance read pairs:", edit_distance_filt)
    print("Valid mapq score read pairs:", mapq_score_filt)

    plot_insert_size(np.array([insert_size_pass,insert_size_fail], dtype=object), 
                     args.insert_min, args.insert_max, args.input)

    plot_insert_size(insert_size_pass, 
                     args.insert_min, args.insert_max, args.input)

    plot_mapq(np.array([mapq_score_pass, mapq_score_fail], dtype=object), 
              args.mapq_min, args.mapq_max, args.input)

    plot_mismatches(np.array([edit_distance_pass, edit_distance_fail], dtype=object),
                    args.edit_min, args.edit_max, args.input)


def has_valid_edit_distance(read, min, max):
    '''
    If min and max are equal (e.g. 0,0) skip this filter
    '''
    score = int(read.get_tag('NM'))
    if min == max:
        return True
    else:
        return min <= score <= max

def has_valid_insert_size(read, min, max):
    '''
    Can be specified while running bowtie2
    If min and max are equal (e.g. 0,0) skip this filter
    '''
    if min == max:
        return True
    else:
        return min <= abs(read.template_length) <= max

def has_same_chromosome(read, mate, hybrid_reference, keep_trans):
    '''
    Exclude trans reads, if hybrid_reference exclude cross reference reads

    Args:
        read(pysam): read
        mate(pysam): read mate
        hybrid_reference(boolean): if hybrid, remove reads that map between refs
        keep_trans(boolean): keep trans reads within reference
    '''
    if hybrid_reference:
        if keep_trans:
            if read.reference_name.startswith('chr') and mate.reference_name.startswith('chr'):
                return True
            elif not read.reference_name.startswith('chr') and not mate.reference_name.startswith('chr'):
                return True
            else:
                return False
    else:
        if keep_trans:
            return True
        else:
            return read.tid == mate.tid

def is_FR_RF(read, mate):
    '''
    read.is_proper_pair - test if reads are in FR or RF orientation
    '''
    # # read forward, mate reverse
    # if not read.is_reverse and mate.is_reverse:
    #     return read.reference_start <= mate.reference_start
    # # read reverse, mate forward
    # elif read.is_reverse and not mate.is_reverse:
    #     return mate.reference_start <= read.reference_start
    assert read.is_proper_pair == mate.is_proper_pair, 'Pairs do not match, make sure BAM was name sorted'

    return read.is_proper_pair

def has_valid_mapq_score(read, lo, hi):
    return lo <= read.mapping_quality <= hi


################################################################################
# Test functions
################################################################################


def get_mapq_scores(args):
    '''
    [TEST] Get BAM MapQ scores to test plotting function

    Args:
        args(Namespace): argparse parsed arguments
    '''
    mapq_scores = []

    with pysam.AlignmentFile(args.input, 'rb') as input_bam:
        for read in input_bam.fetch(until_eof = True):
            mapq_scores.append(read.mapping_quality)

    return mapq_scores


def get_edit_distance(args):
    '''
    [TEST] Get BAM NM flag (edit distance) to test plotting function

    Args:
        args(Namespace): argparse parsed arguments
    '''
    edit_distance = []
    unmapped = 0
    mapped = 0

    with pysam.AlignmentFile(args.input, 'rb') as input_bam:
        for read in input_bam.fetch(until_eof = True):
            try:
                edit_distance.append(int(read.get_tag('NM')))
                mapped += 1
            except KeyError:
                unmapped += 1
    
    print('Unmapped reads:', unmapped)
    print('Mapped reads:', mapped)

    return edit_distance

def get_insert_size(args):
    '''
    [TEST] Get BAM insert size to test plotting function

    Args:
        args(Namespace): argparse parsed arguments
    '''
    insert_size = []
    insert_size_mapq = []
    trans = 0
    proper_pair = 0
    total_pairs = 0
    unmapped = 0
    proper_pair_2 = 0
    singletons = 0
    long_pair = 0
    with pysam.AlignmentFile(args.input, 'rb') as input_bam:
        reads = input_bam.fetch(until_eof = True)
        current_read = next(reads)

        for read in reads:
            # if has_same_chromosome(read, current_read, True, True):
            #     print(read.reference_name, current_read.reference_name)
            if read.query_name == current_read.query_name:
                total_pairs += 1
                if not read.is_unmapped:
                # if read.is_proper_pair and not read.is_unmapped:
                    if read.is_proper_pair:
                        proper_pair += 1
                    if is_FR_RF(read, current_read):
                        proper_pair_2 += 1
                    if read.tid == current_read.tid:
                        if read.mapping_quality >= 20:
                            insert_size_mapq.append(abs(read.template_length))
                            if abs(read.template_length) >= 1000:
                                long_pair += 1
                        else:
                            insert_size.append(abs(read.template_length))
                        # if read.template_length > 200000000:
                        #     print(str(read))
                        #     print(str(current_read))
                    else:
                        trans += 1
                else:
                    unmapped += 1
                # reset counter
                try:
                    current_read = next(reads)
                except:
                    StopIteration
                    EOF = True
            else:
                singletons += 1
                current_read = read

    print('Trans:', trans)
    print('Total pairs:', total_pairs)
    print('Proper pair:', proper_pair)
    print('Unmapped:', unmapped)
    print('Proper pair 2:', proper_pair_2)
    print('Singletons:', singletons)
    print('Long pair:', long_pair)

    return np.array([insert_size_mapq, insert_size], dtype=object)




if __name__ == '__main__':
    main()