'''
This program removes reads from a BAM file according to the filtering criteria.
1. Removes unmapped reads
2. Removes reads with invalid edit distance
3. Removes reads with invalid MAPQ scores
4. Removes reads with invalid insert size (Paired-end)
5. Removes reads on different chromosomes (Paired-end)
6. Removes reads in orientation other than FR (Paired-end)
The script requires that the reads are name sorted, not coordinate sorted
'''


import pysam
import argparse
import sys
import matplotlib as mpl
mpl.use('Agg') # non-interactive backend required for cluster nodes
import matplotlib.pyplot as plt
import numpy as np
# %matplotlib inline
from contextlib import contextmanager
from collections import OrderedDict, defaultdict
import re 
import tempfile
import os
import shutil
import io
import base64
from natsort import index_natsorted


# TODO: Test single-end filtering

def parse_arguments(args=None):

    parser = argparse.ArgumentParser(description = print(__doc__))

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
                                'to 0. Default = 0'))
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
                        help = ('Coordinate sort the final BAM output.'))
    parser.add_argument('--ignore_orientation', action = 'store_true',
                        help = ('Keep paired end reads not in FR or RF orientation'))
    parser.add_argument('--keep_trans', action = 'store_true',
                        help = ('Keep paired end reads not on same chromosome'))
    parser.add_argument('--no_plot', action = 'store_true',
                        help = ('Do not produce any plots'))
    parser.add_argument('--hybrid_reference', action = 'store_true',
                        help = ('Specifies if alignment was performed with a hybrid ' +
                                'reference (Ref 1 = 1, 2, 3, etc.; Ref 2 = chr1, chr2, chr3, etc.'))
    parser.add_argument('--mark_duplicate', action = 'store_true',
                        help = ('Read pairs with identical R1 chrom and coordinates will ' +
                                'be marked as duplicates [Paired only]'))
    parser.add_argument('--QC_limit', action = 'store', metavar = 'X',
                    type = int, default = 0,
                    help = ('Number of read (pairs) to limit QC. Useful to limit running ' +
                            'time on very large BAM files. Default = 0'))
    parser.add_argument('-f', '--filter', action = 'store', metavar = 'STR',
                        help = 'Additional chromosomes to exclude, if more than one add ' +
                               'in a comma separated list')
    parser.add_argument('--filter_alt', action = 'store_true',
                        help = ('Filter out alt chromosomes'))
    parser.add_argument('--filter_regions', action = 'store', metavar= 'FILE',
                        help = ('A bed file with regions from which to exclude any overlapping reads.'))
    parser.add_argument('--natsort_off', action = 'store_true',
                        help = ('Turn off natural sorting of chromosomes. (Useful for yeast with Roman numeral chromosome names)'))


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

    if not args.paired and args.mark_duplicate:
        print('Mark duplicate option ignored for single-end alignment')

    return args


def main():
    args = parse_arguments()

    # test args
    # args = parse_arguments('--edit_max 100 --edit_min 0 --hybrid_reference ' \
    #                          '-i /mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1/split/subset/yi293_GAACCG.CTTTCTCTCGACTTG.human.bam ' \
    #                          '-o /mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1/split/subset/yi293_GAACCG.CTTTCTCTCGACTTG.human.dedup.bam ' \
    #                          '--threads 4 ' \
    #                          '--max_ratio 0.5 ' \
    #                          '--paired ' \
    #                          '--QC_limit 10000 ' \
    #                          '--mark_duplicate'.split())
    
    # args = parse_arguments('--edit_max 4 ' \
    #                        '--edit_min 0 ' \
    #                        '-i /mnt/data/Projects_sync/sciStrand-seq/yy401_CGCTTG.PE.bwa.hg19.markdup.bam ' \
    #                        '-o /mnt/data/Projects_sync/sciStrand-seq/yy401_CGCTTG.PE.bwa.hg19.markdup.filt.bam ' \
    #                        '--threads 4 ' \
    #                        '--insert_max 2000 ' \
    #                        '--insert_min 0 ' \
    #                        '--mapq_max 255 ' \
    #                        '--mapq_min 20 ' \
    #                        '--sort ' \
    #                        '--filter_alt ' \
    #                        '--paired '.split())
    
    # '-i /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/yy401_CGCTTG.PE.bwa.hg19.markdup.bam ' \
    # '-o /mnt/d/SynologyDriveCloud/Projects/sciStrand-seq/yy401_CGCTTG.PE.bwa.hg19.markdup.filt.bam ' \
    # '--filter_regions /mnt/d/SynologyDriveCloud/Projects/scistrandfig/inst/extdata/blacklists/hg19-blacklist.v2.ensembl.bed ' \
    #    '--filter NC_007605 ' \

    # '-o /mnt/data/nextseq190419/yi293_GAACCG_1min_UV_with_USER_HAP1/split/subset/filtered_0.5/test.bam ' \
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


def region_filter(bam_input, bam_output, threads, filter_regions):
    '''
    Filter reads overlapping regions provided as a bed file

    Args:
        bam_intput(str): BAM input path
        bam_output(str): filtered BAM output path
        threads(int): number of threads to use
        filter_regions(str): BED file path with regions to exlude
    '''
    # check for BAM index
    with pysam.AlignmentFile(bam_input, 'rb') as input_file:
        try:
            input_file.check_index()
        except ValueError:
            print('Creating BAM index.')
            pysam.index(bam_input)

    tmp = pysam.view('-@', str(threads), '-L', filter_regions,'-U', bam_output, '-b', bam_input)




def fig_to_base64(fig):
    '''
    Convert matplotlib plot to base64 for html output
    '''
    img = io.BytesIO()
    fig.savefig(img, format='svg',
                bbox_inches='tight')
    img.seek(0)

    return base64.b64encode(img.getvalue())


def write_html_report(out_path, image_list):
    '''
    Args:
        out_path(str): output html path
        image_list(list): List of base64 encoded images to save in report
    '''
    html = ''
    for img in image_list:
        img_dec = img.decode('utf-8')
        html += f'<img src="data:image/svg+xml;base64, {img_dec}">\n'
    with open(out_path, 'w') as html_out:
        html_out.write(html)


def plot_insert_size(insert_size, is_min, is_max, max_pass=2000):
    '''
    Plot stacked bar plot with insert size with other filters overlaid

    Args:
        insert_size(np.array): numpy array of values to plot
        is_min(int): minimum insert size
        is_max(int): max insert size
        max_pass(int): for pass plots, maximum x-axis value
    '''

    if type(insert_size) == list:
        # out_png = save_plt + '.insert_size_pass.png'
        plt.hist(insert_size, bins=100, density=False, histtype='stepfilled', range=(0,max_pass))
        plt.title('Insert size filter pass')
    else:
        # out_png = save_plt + '.insert_size.png'
        plt.hist(insert_size, bins=100, density=False, histtype='stepfilled', stacked=True, 
                 label=['Pass', 'Fail other', 'Fail'], color=['#1f77b4', '#ff7f0e', '#d62728'])
        plt.legend(prop={'size': 10})
        plt.title('Insert size')
    plt.ylabel('Read count')
    plt.xlabel('Insert size')
    # if both min and max == 0, insert filter ignored
    if is_min > 0:
        plt.axvline(x=is_min, color='r', linestyle='dashed', linewidth=1)
    if is_max > 0:
        plt.axvline(x=is_max, color='r', linestyle='dashed', linewidth=1)

    if max_pass > 0:
        plt.yscale('log')

    out = fig_to_base64(plt)
    plt.close()

    return out


def plot_mapq(mapq_scores, min_mapq, max_mapq):
    '''
    Plot counts histogram of BAM MapQ scores
    with reads that would be filtered out based on other criteria 

    Args:
        mapq_scores(list): list of MapQ scores to plot
        min_mapq(int): Cutoff for minimum allowed MapQ score
        max_mapq(int): Cutoff for Maximum allowed MapQ score
    '''


    plt.hist(mapq_scores,
             bins=max(10, max(mapq_scores[0], default=0), max(mapq_scores[1], default=0)), 
             density=False, histtype='stepfilled', 
             stacked=True, label=['Pass', 'Fail other', 'Fail'], color=['#1f77b4', '#ff7f0e', '#d62728'])
    plt.legend(prop={'size': 10})
    plt.ylabel('Read count')
    plt.xlabel('MapQ score')
    plt.axvline(x=min_mapq, color='r', linestyle='dashed', linewidth=1)
    if max_mapq < 255:
        plt.axvline(x=max_mapq, color='r', linestyle='dashed', linewidth=1)

    out = fig_to_base64(plt)
    plt.close()

    return out


def plot_mismatches(edit_distance, min_nm, max_nm):
    '''
    Plot NM flag, number of mismatches in a read

    Args:
        edit_distance(list): list of mismatches to plot
        min_nm(int): minimum allowed edit distance
        max_nm(int): maximum allowed edit distance
    '''

    fig = plt.figure()
    plt.hist(edit_distance, bins=20, range=(0,20), density=False, histtype='stepfilled', 
    stacked=True, label=['Pass', 'Fail other', 'Fail'], color=['#1f77b4', '#ff7f0e', '#d62728'])
    plt.legend(prop={'size': 10})
    plt.ylabel('Read count')
    plt.xlabel('Edit distance')
    plt.title('NM, mismatched within an aligned read')
    if min_nm > 0:
        plt.axvline(x=min_nm, color='r', linestyle='dashed', linewidth=1)
    if max_nm > 0:
        plt.axvline(x=max_nm, color='r', linestyle='dashed', linewidth=1)
    plt.xlim(0, 20)
    
    out = fig_to_base64(plt)
    plt.close()

    return out


def plot_soft_clipped(soft_clipped, read='', ratio_filt=False):
    '''
    Plot histogram of soft clipped bases

    Args:
        soft_clipped(np.array): list of soft clipped bases to plot
        read(str): read 1 or read 2
    '''
    plt.hist(soft_clipped, bins=max(10,max([max(i, default=0) for i in soft_clipped])), 
             density=False, histtype='stepfilled',
             stacked=True, label=['Pass', 'Fail other', 'Fail'], color=['#1f77b4', '#ff7f0e', '#d62728'])
    if ratio_filt:
        plt.title('Soft clipping filter pass ' + read)
    else:
        plt.title('Soft clipping ' + read)
    plt.legend(prop={'size': 10})
    plt.yscale('log')
    plt.ylabel('Read count')
    plt.xlabel('Soft clipped bases')

    out = fig_to_base64(plt)
    plt.close()

    return out





def soft_clipped_vs_matched(sc, m, slope):
    '''
    Plot heatmap of soft clipped vs matched bases

    Args:
        sc(np.array): soft clipped 
        m(np.array): matched
        slope(double): slope of line on plot
    '''
    # out_png = save_plt + '.soft_clipping_vs_matched.png'
    heatmap, xedges, yedges = np.histogram2d(m, sc, bins=100)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    im = plt.imshow(np.log2(heatmap.T+1), extent=extent, origin='lower')
    plt.ylabel('Soft clipped')
    plt.xlabel('Matched')
    plt.title('Soft clipping vs matched')
    xl=plt.xlim()
    plt.axline((0, 0), slope=slope)
    plt.xlim(xl[0])
    plt.colorbar(im)

    out = fig_to_base64(plt)
    plt.close()

    return out


def matched_vs_matched(m_r1, m_r2):
    '''
    Plot heatmap of matched read 1 bases vs matched read 2 bases

    Args:
        m_r1(np.array): matched read 1
        m_r2(np.array): matched read 2
    '''
    # out_png = save_plt + '.matched_r1_vs_matched_r2.png'
    heatmap, xedges, yedges = np.histogram2d(m_r2, m_r1, bins=50)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

    im = plt.imshow(np.log2(heatmap.T+1), extent=extent, origin='lower')
    plt.ylabel('Matched read 1')
    plt.xlabel('Matched read 2')
    plt.title('Matched r1 vs matched r2')
    plt.colorbar(im)

    out = fig_to_base64(plt)
    plt.close()

    return out


def none_context(a=None):
    '''
    contextmanager can only be used once, so need to define it for multiple uses
    '''
    return contextmanager(lambda: (x for x in [a]))()



def get_bam_headers(args, bam_input=None):
    '''Make new BAM header for BAM files

        Note:
            SN - chromosome name
            LN - chromosome size
           
        Return:
            ordered header dicts for pysam
    '''
    if bam_input is None:
        bam_input = args.input
    #get bam header
    with pysam.AlignmentFile(bam_input, 'rb') as input_file:
        bam_header = input_file.header.to_dict()

    # Additional filter
    if args.filter:
        search_str = '^' + '$|^'.join(args.filter.split(',')) + '$'

    #bam header
    
    SQ_1 = []
    for chrom in bam_header['SQ']:
        if args.filter:
            ad_filter = bool(re.search(search_str, chrom['SN']))
        else:
            ad_filter = False 
        # Exclude alt chromosomes (anything that doesn't end in .1 .2 _ etc.)
        if args.filter_alt:
            if not bool(re.search('_|\.\d$', chrom['SN'])) and not ad_filter:
                SQ_1.append(chrom)
        elif not ad_filter:
            SQ_1.append(chrom)


    if args.natsort_off == False:
        # Natural sort to order chromosomes (ensembl fastq not natural sorted)
        SQ_1 = [SQ_1[i] for i in index_natsorted([i['SN'] for i in SQ_1])]

    bam_1 = OrderedDict()

    # Key order = HD SQ RG PG
    for k in bam_header.keys():
        # Add all fields before SQ
        if k == 'SQ':
            bam_1[k] = SQ_1
        else:
            bam_1[k] = bam_header[k]

    return (bam_1)


def filter_single_reads(args):

    if args.filter_regions:
        if os.path.exists(args.filter_regions):
            tmp_dir = tempfile.mkdtemp()
            bam_input = os.path.join(tmp_dir, 'region_filt.bam')
            region_filter(args.input, bam_input, args.threads, args.filter_regions)
    else:
        bam_input = args.input

    if args.sort:
        tmp_dir = tempfile.mkdtemp()
        bam_filt = os.path.join(tmp_dir, 'filtered.bam')
    else:
        bam_filt = args.output

    valid_reads = 0
    invalid_reads = 0
    unmapped = 0
    total = 0
    additional_filt = 0
    no_nm_tag = 0 

    # pass fail based on other filters
    mapq_score_fail = []
    mapq_score_fail_other = []
    mapq_score_pass = []
    edit_distance_fail = []
    edit_distance_fail_other = []
    edit_distance_pass = []
    soft_clip_fail = []
    soft_clip_fail_other = []
    soft_clip_pass = []

    sc = []
    m = []

    bh_1 = get_bam_headers(args, bam_input)

    bh_1_tid_pos = defaultdict()
    for i, item in enumerate(bh_1['SQ']):
        bh_1_tid_pos[item['SN']] = i

    # Additional filter
    if args.filter:
        search_str = '^' + '$|^'.join(args.filter.split(',')) + '$'

    with pysam.AlignmentFile(bam_input, "rb") as input_bam, \
         (pysam.AlignmentFile(bam_filt, "wb", header = bh_1) if args.output else none_context()) as output_bam, \
         (pysam.AlignmentFile(args.output_failed, "wb", template = input_bam) if args.output_failed else none_context()) as output_failed_bam:

        for read in input_bam.fetch(until_eof = True):
            
            if args.QC_limit > 0 and total > args.QC_limit:
                break
            
            total += 1

            valid_edit_distance = has_valid_edit_distance(read, args.edit_min, args.edit_max),
            valid_mapq_score = has_valid_mapq_score(read, args.mapq_min, args.mapq_max)
            valid_ratio = has_valid_soft_clip_matched(args.max_ratio, read)

            # For MAPQ plotting
            if all([not read.is_unmapped,
                    valid_edit_distance,
                    valid_mapq_score,
                    valid_ratio]):
                mapq_score_pass.append(read.mapping_quality)
            elif all([not read.is_unmapped,
                      valid_edit_distance,
                      valid_ratio]):
                mapq_score_fail_other.append(read.mapping_quality)
            else:
                mapq_score_fail.append(read.mapping_quality)

            # Edit distance plotting
            try:
                if all([not read.is_unmapped,
                        valid_edit_distance,
                        valid_mapq_score,
                        valid_ratio]):
                    edit_distance_pass.append(int(read.get_tag('NM')))
                elif all([not read.is_unmapped,
                        valid_mapq_score,
                        valid_ratio]):
                    edit_distance_fail_other.append(int(read.get_tag('NM')))
                else:
                    edit_distance_fail.append(int(read.get_tag('NM')))
            except KeyError:
                # tag is not present in BAM file
                no_nm_tag += 1 
                

            # Soft clipping matched plotting
            sc.append(soft_clipped(read.cigarstring))
            m.append(matched(read.cigarstring))

            if all([not read.is_unmapped,
                    valid_edit_distance,
                    valid_mapq_score,
                    valid_ratio]):
                soft_clip_pass.append(soft_clipped(read.cigarstring))
            elif all([not read.is_unmapped,
                    valid_edit_distance,
                    valid_mapq_score]):
                soft_clip_fail_other.append(soft_clipped(read.cigarstring))
            else:
                soft_clip_fail.append(soft_clipped(read.cigarstring))
                
            
            # filter by chromosome names
            if args.filter:
                ad_filter = any(bool(re.search(search_str, read.reference_name)))
            else:
                ad_filter = False 

            if args.filter_alt:
                if all([not bool(re.search('_|\.\d$', read.reference_name)),
                        not ad_filter]):
                    keep_chrom = True
                else:
                    keep_chrom = False
                    additional_filt += 1
            else:
                if not ad_filter:
                    keep_chrom = True 
                else:
                    keep_chrom = False
                    additional_filt += 1
                    
            if all([not read.is_unmapped,
                    valid_edit_distance,
                    valid_mapq_score,
                    valid_ratio,
                    keep_chrom]):
                
                # update tid
                read.reference_id = bh_1_tid_pos[read.reference_name]

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
    print("Additional filtered reads: ", additional_filt)

    # sort and clean up tmp files 
    if args.sort:
        # sort final output BAM
        coordinate_sort(bam_filt, args.output, args.threads)
        # delete temp folder
        shutil.rmtree(tmp_dir)

    if not args.no_plot:

        if args.output is None:
            html_out = args.input + '.html'
        else:
            html_out = args.output + '.html'

        p1 =plot_mapq(np.array([mapq_score_pass, mapq_score_fail_other, mapq_score_fail], dtype=object), 
                     args.mapq_min, args.mapq_max)

        if no_nm_tag == 0:
            p2 = plot_mismatches(np.array([edit_distance_pass, edit_distance_fail_other, edit_distance_fail], dtype=object),
                            args.edit_min, args.edit_max)
        else:
            p2 = fig_to_base64(plt.figure())

        p3 = plot_soft_clipped(np.array([soft_clip_pass, soft_clip_fail_other, soft_clip_fail], dtype=object))


        p4 = soft_clipped_vs_matched(sc, m, args.max_ratio)

        write_html_report(html_out, [p1, p2, p3, p4])



def read_pair_generator(bam, region_string=None):
    '''
    Generate read pairs in a BAM file or within a region string.
    Reads are added to read_dict until a pair is found.

    Adapted from: https://www.biostars.org/p/306041/
    '''
    secondary_supplementary = 0
    read_dict = defaultdict(lambda: [None, None])
    for read in bam.fetch(until_eof = True, region=region_string):
        if read.is_secondary or read.is_supplementary: #not read.is_proper_pair or
            secondary_supplementary += 1
            continue
        qname = read.query_name
        if qname not in read_dict:
            if read.is_read1:
                read_dict[qname][0] = read
            else:
                read_dict[qname][1] = read
        else:
            if read.is_read1:
                yield read, read_dict[qname][1]
            else:
                yield read_dict[qname][0], read
            del read_dict[qname]
    
    singletons = len(read_dict)

    print('Reads without a pair: ', singletons)
    print('Secondary or supplementary alignments: ', secondary_supplementary)


def filter_paired_reads(args):
    '''
    Main filtering of paired end reads

    Note:
        BAM needs to be name sorted for speed
        mate: This method is too slow for high-throughput processing. 
        If a read needs to be processed with its mate, work from a 
        read name sorted file or, better, cache reads.
    '''
    
    if args.filter_regions:
        if os.path.exists(args.filter_regions):
            tmp_dir = tempfile.mkdtemp()
            bam_input = os.path.join(tmp_dir, 'region_filt.bam')
            region_filter(args.input, bam_input, args.threads, args.filter_regions)
        else:
            print('Region filter path not valid.')
    else:
        bam_input = args.input

    if args.sort:
        tmp_dir = tempfile.mkdtemp()
        bam_filt = os.path.join(tmp_dir, 'filtered.bam')
    else:
        bam_filt = args.output

    total = 0
    written_reads = 0
    valid_reads = 0
    invalid_reads = 0
    valid_count = 0
    trans_count = 0
    unmapped_reads = 0
    insert_size_filt = 0
    edit_distance_filt = 0
    mapq_score_filt = 0
    soft_clip_filt = 0
    additional_filt = 0
    no_nm_tag = 0 

    # pass fail based on other filters
    mapq_score_fail = []
    mapq_score_fail_other = []
    mapq_score_pass = []
    edit_distance_fail = []
    edit_distance_fail_other = []
    edit_distance_pass = []
    insert_size_fail = []
    insert_size_fail_other = []
    insert_size_pass = []
    
    soft_clip_fail_r1 = []
    soft_clip_fail_r2 = []
    soft_clip_fail_other_r1 = []
    soft_clip_fail_other_r2 = []
    # soft_clip_ratio_fail_r1 = []
    # soft_clip_ratio_fail_r2 = []

    soft_clip_pass_r1 = []
    soft_clip_pass_r2 = []
    # soft_clip_ratio_pass_r1 = []
    # soft_clip_ratio_pass_r2 = []
    
    sc_r1 = []
    sc_r2 = []
    m_r1 = []
    m_r2 = []

    # for deduplication
    unique_r1 = set()
    duplicates = 0

    bh_1 = get_bam_headers(args, bam_input)

    bh_1_tid_pos = defaultdict()
    for i, item in enumerate(bh_1['SQ']):
        bh_1_tid_pos[item['SN']] = i

    # Additional filter
    if args.filter:
        search_str = '^' + '$|^'.join(args.filter.split(',')) + '$'

    with pysam.AlignmentFile(bam_input, 'rb') as input_bam, \
        (pysam.AlignmentFile(bam_filt, 'wb', header = bh_1) if args.output else none_context()) as output_bam, \
        (pysam.AlignmentFile(args.output_failed, 'wb', template = input_bam) if args.output_failed else none_context()) as output_failed_bam:

        for read, mate in read_pair_generator(input_bam):

            if args.QC_limit > 0 and total > args.QC_limit:
                break
            
            total += 1
            
            assert read.query_name == mate.query_name, f'{read.query_name} and {mate.query_name} do not match'

            # only check mapped reads
            if all([not read.is_unmapped,
                    not mate.is_unmapped]):

                if args.ignore_orientation:
                    valid_orientation = True
                else:
                    valid_orientation = is_FR_RF(mate, read)
                    if valid_orientation:
                        valid_count += 1
                    
                no_trans_reads = has_same_chromosome(mate, read, 
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
                                        has_valid_edit_distance(mate, 
                                                                args.edit_min, 
                                                                args.edit_max)])
            
                if valid_edit_distance:
                    edit_distance_filt += 1

                valid_mapq_score = all([has_valid_mapq_score(read, 
                                                            args.mapq_min, 
                                                            args.mapq_max),
                                        has_valid_mapq_score(mate, 
                                                            args.mapq_min, 
                                                            args.mapq_max)])
                if valid_mapq_score:
                    mapq_score_filt += 1

                valid_ratio = has_valid_soft_clip_matched(args.max_ratio, 
                                                            mate,
                                                            read)
                if valid_ratio:
                    soft_clip_filt += 1


                # filter by chromosome names (if True will be filtered)
                if args.filter:
                    ad_filter = any([bool(re.search(search_str, read.reference_name)),
                                    bool(re.search(search_str, mate.reference_name))])    
                else:
                    ad_filter = False 

                if args.filter_alt:
                    if all([not bool(re.search('_|\.\d$', read.reference_name)),
                            not bool(re.search('_|\.\d$', mate.reference_name)),
                            not ad_filter]):
                        keep_chrom = True
                    else:
                        keep_chrom = False
                        additional_filt += 1
                else:
                    if not ad_filter:
                        keep_chrom = True
                    else:
                        keep_chrom = False
                        additional_filt += 1



                # For MAPQ plotting
                if all([not read.is_unmapped,
                        not mate.is_unmapped,
                        valid_mapq_score,
                        valid_edit_distance,
                        valid_orientation,
                        no_trans_reads,
                        valid_insert_size,
                        valid_ratio]):
                    mapq_score_pass.append(min(read.mapping_quality, 
                                               mate.mapping_quality))
                elif all([not read.is_unmapped,
                          not mate.is_unmapped,
                          valid_mapq_score]):
                    mapq_score_fail_other.append(min(read.mapping_quality, 
                                                     mate.mapping_quality))
                else:
                    mapq_score_fail.append(min(read.mapping_quality, 
                                            mate.mapping_quality))

                # Edit distance plotting
                try:
                    if all([not read.is_unmapped,
                            not mate.is_unmapped,
                            valid_mapq_score,
                            valid_edit_distance,
                            valid_orientation,
                            no_trans_reads,
                            valid_insert_size,
                            valid_ratio]):
                        edit_distance_pass.append(max(int(read.get_tag('NM')), 
                                                int(mate.get_tag('NM'))))
                    elif all([not read.is_unmapped,
                            not mate.is_unmapped,
                            valid_edit_distance]):
                        edit_distance_fail_other.append(max(int(read.get_tag('NM')), 
                                                int(mate.get_tag('NM'))))
                    else:
                        edit_distance_fail.append(max(int(read.get_tag('NM')), 
                                                int(mate.get_tag('NM'))))
                except KeyError:
                    no_nm_tag += 1 

                # Insert size plotting
                if all([not read.is_unmapped,
                        not mate.is_unmapped,
                        valid_mapq_score,
                        valid_edit_distance,
                        valid_orientation,
                        no_trans_reads,
                        valid_insert_size,
                        valid_ratio]):
                    insert_size_pass.append(abs(read.template_length))
                elif all([not read.is_unmapped,
                        not mate.is_unmapped,
                        valid_insert_size]):
                    insert_size_fail_other.append(abs(read.template_length))
                else:
                    insert_size_fail.append(abs(read.template_length))

                # Soft clipping vs matched plotting
                sc_r1.append(soft_clipped(read.cigarstring))
                sc_r2.append(soft_clipped(mate.cigarstring))
                m_r1.append(matched(read.cigarstring))
                m_r2.append(matched(mate.cigarstring))

                if all([not read.is_unmapped,
                        not mate.is_unmapped,
                        valid_mapq_score,
                        valid_edit_distance,
                        valid_orientation,
                        no_trans_reads,
                        valid_insert_size,
                        valid_ratio]):
                    soft_clip_pass_r1.append(soft_clipped(read.cigarstring))
                    soft_clip_pass_r2.append(soft_clipped(mate.cigarstring))
                elif all([not read.is_unmapped,
                          not mate.is_unmapped,
                          valid_ratio]):
                    soft_clip_fail_other_r1.append(soft_clipped(read.cigarstring))
                    soft_clip_fail_other_r2.append(soft_clipped(mate.cigarstring))
                else:
                    soft_clip_fail_r1.append(soft_clipped(read.cigarstring))
                    soft_clip_fail_r2.append(soft_clipped(mate.cigarstring))



                # Filtering and BAM write out
                if all([not read.is_unmapped,
                        not mate.is_unmapped,
                        valid_edit_distance,
                        valid_mapq_score,
                        valid_orientation,
                        no_trans_reads,
                        valid_insert_size,
                        valid_ratio,
                        keep_chrom]):
                    
                    # update tid                
                    read.reference_id = bh_1_tid_pos.get(read.reference_name)
                    read.next_reference_id = bh_1_tid_pos.get(read.next_reference_name)
                    # update mate tid
                    mate.reference_id = bh_1_tid_pos.get(mate.reference_name)
                    mate.next_reference_id = bh_1_tid_pos.get(mate.next_reference_name)
                    
                    if args.mark_duplicate:
                        # Deduplicated using read 1 coordinates
                        read_id = '-'.join([str(read.reference_name), 
                                            str(read.reference_start), 
                                            str(read.reference_end)])

                        if read_id in unique_r1:
                            read.is_duplicate = True
                            duplicates += 1
                        else:
                            read.is_duplicate = False
                            unique_r1.add(read_id)

                    valid_reads += 1 
                    if args.output:
                        output_bam.write(read)
                        output_bam.write(mate)
                        written_reads += 1 
                # elif any([read.is_unmapped,
                #           mate.is_unmapped]):
                #     unmapped_reads += 1
                else:
                    invalid_reads += 1
                    if args.output_failed:
                        output_failed_bam.write(read)
                        output_failed_bam.write(mate)
            else:
                unmapped_reads += 1


    print("Valid read pairs: ", valid_reads)
    print("Written out read pairs: ", written_reads)
    print("Invalid read pairs filtered: ", invalid_reads)
    print("Unmapped reads pairs: ", unmapped_reads)
    if args.keep_trans and args.hybrid_reference:
        print("Hybrid read pairs count: ", trans_count)
    elif not args.keep_trans and args.hybrid_reference:
        print("Trans or hybrid read pairs count: ", trans_count)
    elif not args.keep_trans and not args.hybrid_reference:
        print("Trans read pairs count: ", trans_count)

    if not args.ignore_orientation:
        print("Proper orientation read pairs: ", valid_count)

    print("Valid insert size read pairs: ", insert_size_filt)
    print("Valid edit distance read pairs: ", edit_distance_filt)
    print("Valid mapq score read pairs: ", mapq_score_filt)
    print("Valid soft clipping/matched ratio read pairs: ", soft_clip_filt)
    print("Duplicate read pairs: ", duplicates)
    print("Additional filtered by chromosome name: ", additional_filt)

    # sort and clean up tmp files 
    if args.sort:
        # sort final output BAM
        if args.output:
            coordinate_sort(bam_filt, args.output, args.threads)
        # delete temp folder
        shutil.rmtree(tmp_dir)


    if not args.no_plot:
        if args.output is None:
            html_out = args.input + '.html'
        else:
            html_out = args.output + '.html'

        p1 = plot_insert_size(np.array([insert_size_pass,insert_size_fail_other, insert_size_fail], dtype=object), 
                        args.insert_min, args.insert_max)

        p2 = plot_insert_size(insert_size_pass, args.insert_min, args.insert_max, 
                              max_pass=args.insert_max)

        p3 = plot_mapq(np.array([mapq_score_pass, mapq_score_fail_other, mapq_score_fail], dtype=object), 
                       args.mapq_min, args.mapq_max)

        if no_nm_tag == 0:
            p4 = plot_mismatches(np.array([edit_distance_pass, edit_distance_fail_other, edit_distance_fail], dtype=object),
                                args.edit_min, args.edit_max)
        else:
            p4 = fig_to_base64(plt.figure())

        p5 = plot_soft_clipped(np.array([soft_clip_pass_r1, soft_clip_fail_other_r1, soft_clip_fail_r1], dtype=object),
                               read='read 1')

        p6 = plot_soft_clipped(np.array([soft_clip_pass_r2, soft_clip_fail_other_r2, soft_clip_fail_r2], dtype=object),
                               read='read 2')


        p7 = plot_soft_clipped(np.array([np.array(soft_clip_pass_r1) + np.array(soft_clip_pass_r2), 
                                         np.array(soft_clip_fail_other_r1) + np.array(soft_clip_fail_other_r2),
                                         np.array(soft_clip_fail_r1) + np.array(soft_clip_fail_r2)], dtype=object),
                               read='read 1 and 2')

        p8 = soft_clipped_vs_matched(np.array(sc_r1) + np.array(sc_r2), 
                                     np.array(m_r1) + np.array(m_r2), args.max_ratio)

        p9 = matched_vs_matched(m_r1, m_r2)

        write_html_report(html_out, [p1, p2, p3, p4, p5, p6, p7, p8, p9])




def has_valid_edit_distance(read, min, max):
    '''
    If min and max are equal (e.g. 0,0) skip this filter
    '''
    try:
        score = int(read.get_tag('NM'))
        if min == max:
            return True
        else:
            return min <= score <= max
    except KeyError:
        return True

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

    if ratio == 0:
        return True
    else:
        return r < ratio



################################################################################
# Test / unused functions
################################################################################

# Singleton reads filtered:  0
# Secondary or supplementary alignments:  103
# Valid read pairs written out:  24029
# Invalid read pairs filtered:  1809
# Unmapped reads pairs:  136
# Trans or hybrid read pairs count:  79
# Proper orientation read pairs:  25728
# Valid insert size read pairs:  25826
# Valid edit distance read pairs:  25838
# Valid mapq score read pairs:  25838
# Valid soft clipping/matched ratio read pairs:  24108
# Duplicate read pairs:  3105

def filter_paired_reads_old(args):
    '''
    !! secondary alignments mess this up !!

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
    secondary_supplemental = 0
    no_nm_tag = 0

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

        while current_read.is_secondary or current_read.is_supplementary:
            secondary_supplemental += 1
            current_read = next(reads)

        for read in reads:
            while read.is_secondary or read.is_supplementary:
                secondary_supplemental += 1
                read = next(reads)

            if read.query_name == current_read.query_name:
                
                # only check primary mapped reads
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
                    try:
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
                    except KeyError:
                        no_nm_tag += 1 

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
                else:
                    unmapped_reads += 1
                    
                # reset counter
                try:
                    current_read = next(reads)
                except:
                    StopIteration
                    EOF = True
            else:
                if current_read.is_secondary or current_read.is_supplementary or \
                read.is_secondary or read.is_supplementary:
                    secondary_supplemental += 1
                else:
                    singletons += 1
                current_read = read
    if not EOF:
        if current_read.is_secondary or current_read.is_supplementary or \
           read.is_secondary or read.is_supplementary:
            secondary_supplemental += 1
        else:
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
    print("Secondary or supplemental reads: ", secondary_supplemental)

    # sort and clean up tmp files 
    if not bam_sorted or args.sort:
        # sort final output BAM
        if args.output:
            coordinate_sort(bam_filt, args.output, args.threads)
        # delete temp folder
        shutil.rmtree(tmp_dir)

    if not args.no_plot:
        if args.output is None:
            html_out = args.input + '.html'
        else:
            html_out = args.output + '.html'

        p1 = plot_insert_size(np.array([insert_size_pass,insert_size_fail], dtype=object), 
                        args.insert_min, args.insert_max)

        p2 = plot_insert_size(insert_size_pass, args.insert_min, args.insert_max, 
                              max_pass=args.insert_max)

        p3 = plot_mapq(np.array([mapq_score_pass, mapq_score_fail], dtype=object), 
                       args.mapq_min, args.mapq_max)

        
        p4 = plot_mismatches(np.array([edit_distance_pass, edit_distance_fail], dtype=object),
                            args.edit_min, args.edit_max)

        p5 = plot_soft_clipped(np.array([soft_clip_pass, soft_clip_fail], dtype=object))

        p6 = plot_soft_clipped(np.array([soft_clip_ratio_pass, soft_clip_ratio_fail], dtype=object), 
                               ratio_filt=True)

        p7 = soft_clipped_vs_matched(sc, np.array(m_r1) + np.array(m_r2), args.max_ratio)

        p8 = matched_vs_matched(m_r1, m_r2)

        write_html_report(html_out, [p1, p2, p3, p4, p5, p6, p7, p8])


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

        for read, mate in read_pair_generator(input_bam):
            if total_pairs < 10:
                print(read.flag)
                read.is_duplicate = True
                print(read.is_duplicate)
                print(read.query_name, read.flag)

                print('-'.join([str(read.reference_name), str(read.reference_start), str(read.reference_end)]))
                print()
                # print(read.get_tags())
                total_pairs += 1
                # print(read.cigarstring, read.is_unmapped, mate.cigarstring, mate.is_unmapped)
                
                if read.query_name == mate.query_name:
                    if all([not read.is_unmapped,
                            not mate.is_unmapped]):
                    
                        # sum soft clipped bases of both pairs
                        sc.append(soft_clipped(mate.cigarstring) +  soft_clipped(read.cigarstring))
                        m1_value = matched(mate.cigarstring)
                        m2_value = matched(read.cigarstring)
                        m.append(m1_value + m2_value)
                        m1.append(m1_value)
                        m2.append(m2_value)
                    else:
                        unmapped_pairs += 1
            else:
                break


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
        for read, mate in read_pair_generator(input_bam):
            # if has_same_chromosome(read, mate, True, True):
            #     print(read.reference_name, mate.reference_name)
            if read.query_name == mate.query_name:
                total_pairs += 1
                if not read.is_unmapped:
                # if read.is_proper_pair and not read.is_unmapped:
                    if read.is_proper_pair:
                        proper_pair += 1
                    if is_FR_RF(read, mate):
                        proper_pair_2 += 1
                    if read.tid == mate.tid:
                        if read.mapping_quality >= 20:
                            insert_size_mapq.append(abs(read.template_length))
                            if abs(read.template_length) >= 1000:
                                long_pair += 1
                        else:
                            insert_size.append(abs(read.template_length))
                        # if read.template_length > 200000000:
                        #     print(str(read))
                        #     print(str(mate))
                    else:
                        trans += 1
                else:
                    unmapped += 1
                # reset counter
                try:
                    mate = next(reads)
                except:
                    StopIteration
                    EOF = True
            else:
                singletons += 1
                mate = read

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