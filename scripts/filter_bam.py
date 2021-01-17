
import pysam
import argparse
import sys
import matplotlib as mpl
mpl.use('Agg') # non-interactive backend required for cluster nodes
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline
from contextlib import contextmanager
import re 
import tempfile
import os
import shutil

# TODO: Test single-end filtering 

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
                        required=True, help = 'Input BAM file')
    parser.add_argument('-o', '--output', action = 'store', metavar = 'FILE',
                        help = 'Output BAM file containing reads that pass. ' +
                        'If not specified, no output produced (useful for QC plots only)')
    parser.add_argument('--output_failed', action = 'store', metavar = 'FILE',
                        help = 'Output BAM file containing reads that failed')
    parser.add_argument('--edit_max', action = 'store', metavar = 'Y',
                        type = int, default = 0,
                        help = ('If the edit distance between the read ' +
                                'sequence and the reference sequence falls ' +
                                'within [X, Y], keep the read. Otherwise ' +
                                'remove it. To ignore filter set both min and max ' +
                                'to 0. Default = 0q'))
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
    parser.add_argument('--max_ratio', action = 'store', metavar = 'X',
                        type = float, default = 0,
                        help = ('If the ratio of soft clipped bases over matched is ' +
                                'X, filter out the read. e.g. 50 soft clipped bases to ' +
                                '150 total matched bases would produce a ratio of 0.5, ' +
                                'which if X=0.5, the read(s) would be filtered out.' +
                                'Default = 0'))
    parser.add_argument('--threads', action = 'store', metavar = 'X',
                        type = int, default = 1,
                        help = ('Number of threads to use if supplied file is not name ' +
                                'sorted. Default = 1'))
    parser.add_argument('--paired', action = 'store_true',
                        help = ('Set if the BAM file contains a paired-end ' +
                                'alignment. Defaults to single-read alignment.'))
    parser.add_argument('--sort', action = 'store_true',
                        help = ('Coordinate sort the final BAM output. If input was ' +
                        'coordinate sorted, final output will automatically be coordinate ' +
                        'sorted without setting this flag'))
    parser.add_argument('--ignore_orientation', action = 'store_true',
                        help = ('Keep paired end reads not in FR or RF orientation'))
    parser.add_argument('--keep_trans', action = 'store_true',
                        help = ('Keep paired end reads not on same chromosome'))
    parser.add_argument('--no_plot', action = 'store_true',
                        help = ('Do not produce any plots'))
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
    # args = parse_arguments('--edit_max 100 --edit_min 0 --hybrid_reference ' \
    #                          '-i /mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1/split/subset/yi293_GAACCG.CTTTCTCTCGACTTG.human.bam ' \
    #                          '--threads 4 ' \
    #                          '-o /mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1/split/subset/filtered_0.5/test.bam ' \
    #                          '--max_ratio 0.5 ' \
    #                          '--paired'.split())

    #                          -insert_max 2000 --insert_min 0 --mapq_max 255 --mapq_min 0

    # /mnt/data/sci-l3/test/alignments/yi292/yi292_GTGCAG.PE.bwa.hg38.collate.bam
    # /mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1/split/subset/filtered_0.5/yi293_GAACCG.CATTGTGTCGGAGAC.human.sc.bam

    # mapq = get_mapq_scores(args)
    # ed = get_edit_distance(args)
    # insert_size = get_insert_size(args)

    if args.output:
        os.makedirs(os.path.dirname(args.output), exist_ok=True)

    if args.paired:
        filter_paired_reads(args)
    else:
        filter_single_reads(args)



def valid_bam_sort(bam_input):
    '''
    Check BAM file is collate or name sorted, required for filtering

    Args:
        bam_input(str): BAM input path

    return(boolean): True if name sorted, False if unsorted or coordinate sorted
    '''
    with pysam.AlignmentFile(bam_input, 'rb') as input_file:
        sort_status = input_file.header['HD']
        if sort_status['SO'] == 'coordinate':
            print('BAM coordinate sorted, will collate before filtering.')
            return False
        if sort_status['SO'] == 'unsorted':
            try:
                # collate
                if sort_status['GO'] == 'query':
                    return True
            except KeyError:
                print('BAM not sorted, will collate before filtering.')
                return False


def collate_bam(bam_input, bam_output, threads):
    '''
    Group reads pairs together

    Args:
        bam_input(str): BAM input path
        bam_output(str): collated BAM output path
        threads(int): number of threads to use
    '''
    pysam.collate('-@', str(threads), '-o', bam_output, bam_input)

def coordinate_sort(bam_input, bam_output, threads):
    '''
    Coordinate sort BAM file
    
    Args:
        bam_input(str): BAM input path
        bam_output(str): collated BAM output path
        threads(int): number of threads to use
    '''
    pysam.sort('-@', str(threads), '-o', bam_output, bam_input)

def plot_insert_size(insert_size, is_min, is_max, save_plt, max_pass=2000):
    '''
    Plot stacked bar plot with insert size with other filters overlaid

    Args:
        insert_size(np.array): numpy array of values to plot
        is_min(int): minumin insert size
        is_max(int): max insert size
        save_plt(str): output path name to save plot
        max_pass(int): for pass plots, maximum x-axis value
    '''

    if type(insert_size) == list:
        out_png = save_plt + '.insert_size_pass.png'
        plt.hist(insert_size, bins=100, density=False, histtype='bar', range=(0,max_pass))
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

    plt.hist(mapq_scores, 
             bins=max(max(mapq_scores[0], default=0), max(mapq_scores[1], default=0)), 
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


def plot_soft_clipped(soft_clipped, save_plt, ratio_filt=False):
    '''
    Plot histogram of soft clipped bases

    Args:
        soft_clipped(np.array): list of soft clipped bases to plot
        max_sc(int): maximum allowed soft clipping
        save_plt(str): output path name to save plot
    '''
    if ratio_filt:
        out_png = save_plt + '.soft_clipping_pass.png'
        plt.hist(soft_clipped, bins=200, density=False, histtype='bar',
        stacked=True, label=['Pass', 'Fail'])
        plt.legend(prop={'size': 10})
    else:
        out_png = save_plt + '.soft_clipping.png'
        plt.hist(soft_clipped, bins=200, density=False, histtype='bar',
        stacked=True, label=['Pass', 'Fail'])
        plt.legend(prop={'size': 10})
    plt.yscale('log')
    plt.ylabel('Read count')
    plt.xlabel('Soft clipped bases')
    plt.savefig(out_png, bbox_inches='tight') 
    plt.close()


def soft_clipped_vs_matched(sc, m, slope, save_plt):
    '''
    Plot heatmap of soft clipped vs matched bases

    Args:
        sc(np.array): soft clipped 
        m(np.array): matched
        slope(double): slope of line on plot
        save_plt(str): output path name to save plot
    '''
    out_png = save_plt + '.soft_clipping_vs_matched.png'
    heatmap, xedges, yedges = np.histogram2d(m, sc, bins=100)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    im = plt.imshow(np.log2(heatmap.T+1), extent=extent, origin='lower')
    plt.ylabel('Soft clipped')
    plt.xlabel('Matched')
    xl=plt.xlim()
    plt.axline((0, 0), slope=slope)
    plt.xlim(xl[0])
    plt.colorbar(im)
    plt.savefig(out_png, bbox_inches='tight') 
    plt.close()

def matched_vs_matched(m_r1, m_r2, save_plt):
    '''
    Plot heatmap of matched read 1 bases vs matched read 2 bases

    Args:
        m_r1(np.array): matched read 1
        m_r2(np.array): matched read 2
        save_plt(str): output path name to save plot
    '''
    out_png = save_plt + '.matched_r1_vs_matched_r2.png'
    heatmap, xedges, yedges = np.histogram2d(m_r2, m_r1, bins=50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    im = plt.imshow(np.log2(heatmap.T+1), extent=extent, origin='lower')
    plt.ylabel('Matched read 1')
    plt.xlabel('Matched read 2')
    plt.colorbar(im)
    plt.savefig(out_png, bbox_inches='tight') 
    plt.close()

def none_context(a=None):
    '''
    contextmanager can only be used once, so need to define it for multiple uses
    '''
    return contextmanager(lambda: (x for x in [a]))()

def filter_single_reads(args):

    if args.sort:
        tmp_dir = tempfile.mkdtemp()
        bam_input = args.input
        bam_filt = os.path.join(tmp_dir, 'filtered.bam')

    valid_reads = 0
    invalid_reads = 0
    unmapped = 0

    # pass fail based on other filters
    mapq_score_fail = []
    mapq_score_pass = []
    edit_distance_fail = []
    edit_distance_pass = []
    soft_clip_fail = []
    soft_clip_ratio_fail = []
    soft_clip_pass = []
    soft_clip_ratio_pass = []
    sc = []
    m = []


    with pysam.AlignmentFile(args.input, "rb") as input_bam, \
         (pysam.AlignmentFile(args.output, "wb", template = input_bam) if args.output else none_context()) as output_bam, \
         (pysam.AlignmentFile(args.output_failed, "wb", template = input_bam) if args.output_failed else none_context()) as output_failed_bam:

        for read in input_bam.fetch(until_eof = True):
            valid_edit_distance = has_valid_edit_distance(read, args.edit_min, args.edit_max),
            valid_mapq_score = has_valid_mapq_score(read, args.mapq_min, args.mapq_max)
            valid_ratio = has_valid_soft_clip_matched(args.max_ratio, read)

            # For MAPQ plotting
            if all([not read.is_unmapped,
                    valid_edit_distance,
                    valid_ratio]):
                mapq_score_pass.append(read.mapping_quality)
            elif not read.is_unmapped:
                mapq_score_fail.append(read.mapping_quality)

            # Edit distance plotting
            if all([not read.is_unmapped,
                    valid_mapq_score,
                    valid_ratio]):
                edit_distance_pass.append(int(read.get_tag('NM')))
            elif not read.is_unmapped:
                edit_distance_fail.append(int(read.get_tag('NM')))

            # Soft clipping matched plotting
            sc.append(soft_clipped(read.cigarstring))
            m.append(matched(read.cigarstring))

            if all([not read.is_unmapped,
                    valid_edit_distance,
                    valid_mapq_score]):
                soft_clip_pass.append(soft_clipped(read.cigarstring))
            elif not read.is_unmapped:
                soft_clip_fail.append(soft_clipped(read.cigarstring))
                
            if all([not read.is_unmapped,
                    valid_ratio]):
                soft_clip_ratio_pass.append(soft_clipped(read.cigarstring))   
            elif not read.is_unmapped:
                soft_clip_ratio_fail.append(soft_clipped(read.cigarstring))

            
                    
            if all([not read.is_unmapped,
                    valid_edit_distance,
                    valid_mapq_score,
                    valid_ratio]):

                valid_reads += 1
                if args.output:
                    output_bam.write(read)
            elif read.is_unmapped:
                unmapped += 1
            else:
                invalid_reads += 1
                if args.output_failed:
                    output_failed_bam.write(read)

    print("Valid reads written out:", valid_reads)
    print("Invalid reads filtered:", invalid_reads)
    print("Unmapped reads:", unmapped)

    # sort and clean up tmp files 
    if args.sort:
        # sort final output BAM
        coordinate_sort(bam_filt, args.output, args.threads)
        # delete temp folder
        shutil.rmtree(tmp_dir)

    if not args.no_plot:

        if args.output is None:
            png_out = args.input
        else:
            png_out = args.output

        plot_mapq(np.array([mapq_score_pass, mapq_score_fail], dtype=object), 
                args.mapq_min, args.mapq_max, png_out)

        plot_mismatches(np.array([edit_distance_pass, edit_distance_fail], dtype=object),
                        args.edit_min, args.edit_max, png_out)

        plot_soft_clipped(np.array([soft_clip_pass, soft_clip_fail], dtype=object),
                          png_out)

        plot_soft_clipped(np.array([soft_clip_pass, soft_clip_ratio_fail], dtype=object),
                          png_out, ratio_filt=True)

        soft_clipped_vs_matched(sc, m, args.max_ratio, png_out)


def filter_paired_reads(args):
    '''
    Main filtering of paired end reads

    Note:
        BAM needs to be name sorted for speed
        mate: This method is too slow for high-throughput processing. 
        If a read needs to be processed with its mate, work from a 
        read name sorted file or, better, cache reads.
    '''

    # check if BAM is collate or name sorted
    bam_sorted = valid_bam_sort(args.input)
    if not bam_sorted:
        # create temp dir for collate file
        tmp_dir = tempfile.mkdtemp()
        bam_input = os.path.join(tmp_dir, 'collated.bam')
        bam_filt = os.path.join(tmp_dir, 'filtered.bam')
        collate_bam(args.input, bam_input, args.threads)
    elif args.sort:
        tmp_dir = tempfile.mkdtemp()
        bam_input = args.input
        bam_filt = os.path.join(tmp_dir, 'filtered.bam')
    else:
        bam_input = args.input
        bam_filt = args.output


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
    soft_clip_filt = 0

    # pass fail based on other filters
    mapq_score_fail = []
    mapq_score_pass = []
    edit_distance_fail = []
    edit_distance_pass = []
    insert_size_fail = []
    insert_size_pass = []
    soft_clip_fail = []
    soft_clip_ratio_fail = []
    soft_clip_pass = []
    soft_clip_ratio_pass = []
    sc = []
    m_r1 = []
    m_r2 = []


    with pysam.AlignmentFile(bam_input, 'rb') as input_bam, \
        (pysam.AlignmentFile(bam_filt, 'wb', template = input_bam) if args.output else none_context()) as output_bam, \
        (pysam.AlignmentFile(args.output_failed, 'wb', template = input_bam) if args.output_failed else none_context()) as output_failed_bam:

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

                    valid_ratio = has_valid_soft_clip_matched(args.max_ratio, 
                                                                current_read,
                                                                read)
                    if valid_ratio:
                        soft_clip_filt += 1

                    # For MAPQ plotting
                    if all([not read.is_unmapped,
                            not current_read.is_unmapped,
                            valid_edit_distance,
                            valid_orientation,
                            no_trans_reads,
                            valid_insert_size,
                            valid_ratio]):
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
                            valid_insert_size,
                            valid_ratio]):
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
                            no_trans_reads,
                            valid_ratio]):
                        insert_size_pass.append(abs(read.template_length))
                    elif all([not read.is_unmapped,
                                not current_read.is_unmapped]):
                        insert_size_fail.append(abs(read.template_length))

                    # Soft clipping vs matched plotting
                    sc.append(soft_clipped(current_read.cigarstring) + soft_clipped(read.cigarstring))
                    m_r1.append(matched(current_read.cigarstring))
                    m_r2.append(matched(read.cigarstring))

                    if all([not read.is_unmapped,
                            not current_read.is_unmapped,
                            valid_edit_distance,
                            valid_mapq_score,
                            valid_orientation,
                            no_trans_reads,
                            valid_insert_size]):
                        soft_clip_pass.append(soft_clipped(read.cigarstring)+
                                                soft_clipped(current_read.cigarstring))
                    elif all([not read.is_unmapped,
                                not current_read.is_unmapped]):
                        soft_clip_fail.append(soft_clipped(read.cigarstring)+
                                                soft_clipped(current_read.cigarstring))

                    if all([not read.is_unmapped,
                            not current_read.is_unmapped,
                            valid_ratio]):
                        soft_clip_ratio_pass.append(soft_clipped(read.cigarstring)+
                                                soft_clipped(current_read.cigarstring))
                    elif all([not read.is_unmapped,
                                not current_read.is_unmapped]):
                            soft_clip_ratio_fail.append(soft_clipped(read.cigarstring)+
                                                        soft_clipped(current_read.cigarstring))

                    # Filtering and BAM write out
                    if all([not read.is_unmapped,
                            not current_read.is_unmapped,
                            valid_edit_distance,
                            valid_mapq_score,
                            valid_orientation,
                            no_trans_reads,
                            valid_insert_size,
                            valid_ratio]):
                        
                        valid_reads += 1 
                        if args.output:
                            output_bam.write(current_read)
                            output_bam.write(read)
                    elif any([read.is_unmapped,
                                current_read.is_unmapped]):
                        unmapped_reads += 1
                    else:
                        invalid_reads += 1
                        if args.output_failed:
                            output_failed_bam.write(current_read)
                            output_failed_bam.write(read)

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
    print("Valid soft clipping/matched ratio read pairs:", soft_clip_filt)

    # sort and clean up tmp files 
    if not bam_sorted or args.sort:
        # sort final output BAM
        coordinate_sort(bam_filt, args.output, args.threads)
        # delete temp folder
        shutil.rmtree(tmp_dir)


    if not args.no_plot:
        if args.output is None:
            png_out = args.input
        else:
            png_out = args.output

        plot_insert_size(np.array([insert_size_pass,insert_size_fail], dtype=object), 
                        args.insert_min, args.insert_max, png_out)

        plot_insert_size(insert_size_pass, 
                        args.insert_min, args.insert_max, png_out, max_pass=args.insert_max)

        plot_mapq(np.array([mapq_score_pass, mapq_score_fail], dtype=object), 
                args.mapq_min, args.mapq_max, png_out)

        plot_mismatches(np.array([edit_distance_pass, edit_distance_fail], dtype=object),
                        args.edit_min, args.edit_max, png_out)

        plot_soft_clipped(np.array([soft_clip_pass, soft_clip_fail], dtype=object),
                        png_out)

        plot_soft_clipped(np.array([soft_clip_ratio_pass, soft_clip_ratio_fail], dtype=object), 
                        png_out, ratio_filt=True)

        soft_clipped_vs_matched(sc, np.array(m_r1) + np.array(m_r2), args.max_ratio, png_out)

        matched_vs_matched(m_r1, m_r2, png_out)


def has_valid_edit_distance(read, min, max):
    '''
    If min and max are equal (e.g. 0,0) skip this filter
    '''
    score = int(read.get_tag('NM'))
    if min == max:
        return True
    else:
        return min <= score <= max

def has_valid_insert_size(read, in_min, in_max):
    '''
    Can be specified while running bowtie2
    If min and max are equal (e.g. 0,0) skip this filter
    '''
    if in_min == in_max:
        return True
    else:
        return in_min <= abs(read.template_length) <= in_max

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
            return read.tid == mate.tid
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

def soft_clipped(cigar):
    '''
    Extract soft clipped bases from cigar
    '''
    sc = re.findall("\d+S", cigar)
    if len(sc) > 0:
        total_sc = 0
        for i in sc:
            total_sc += int(i[:-1])
    else:
        total_sc = 0
    return total_sc

def matched(cigar):
    '''
    Extract number of aligned bases from cigar
    '''
    m = re.findall("\d+M", cigar)
    if len(m) > 0:
        total_m = 0
        for i in m:
            total_m += int(i[:-1])
    else:
        total_m = 0
    return total_m


def has_valid_soft_clipping(max_sc, read, mate=None):
    '''
    Filter reads with more than max allowed soft clipping
    '''
    scr = int(soft_clipped(read.cigarstring))
    if mate == None:
        return scr <= max_sc
    else:
        scm = int(soft_clipped(mate.cigarstring))

        return scr+scm <= max_sc

def has_valid_soft_clip_matched(ratio, read, mate=None):
    '''
    Filter reads with a higher ratio of soft-clipped/matched
    '''
    if mate == None:
        sc = soft_clipped(read.cigarstring)
        m = matched(read.cigarstring)
    else:
        sc = soft_clipped(read.cigarstring) + soft_clipped(mate.cigarstring)
        m = matched(read.cigarstring) + matched(mate.cigarstring)

    r = float(sc)/float(m)

    return r < ratio



################################################################################
# Test functions
################################################################################



def get_soft_clipped(args):
    '''
    [TEST] Get BAM length of soft-clipped bases to test plotting function

    Args:
        args(Namespace): argparse parsed arguments
    '''

    total_pairs = 0
    unmapped_pairs = 0
    singletons = 0
    sc = []
    m = []
    m1 = []
    m2 = []

    with pysam.AlignmentFile(args.input, 'rb') as input_bam:
        reads = input_bam.fetch(until_eof = True)
        current_read = next(reads)

        for read in reads:
    
            total_pairs += 1
            # print(read.cigarstring, read.is_unmapped, current_read.cigarstring, current_read.is_unmapped)
            
            if read.query_name == current_read.query_name:
                if all([not read.is_unmapped,
                        not current_read.is_unmapped]):
                
                    # sum soft clipped bases of both pairs
                    sc.append(soft_clipped(current_read.cigarstring) +  soft_clipped(read.cigarstring))
                    m1_value = matched(current_read.cigarstring)
                    m2_value = matched(read.cigarstring)
                    m.append(m1_value + m2_value)
                    m1.append(m1_value)
                    m2.append(m2_value)
                else:
                    unmapped_pairs += 1

                try:
                    current_read = next(reads)
                except:
                    StopIteration
                    EOF = True
            else:
                singletons += 1
                current_read = read
 

    print(f'Total pairs: {total_pairs}')
    print(f'Unmapped pairs: {unmapped_pairs}')
    print(f'Singletons: {singletons}')

    return [sc,m,m1,m2]
    



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