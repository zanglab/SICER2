# Authors: Dustin E Schones, Chongzhi Zang, Weiqun Peng and Keji Zhao

# Modified by: Jin Yong Yoo

import multiprocessing as mp
from functools import partial
from math import *
import sys
import numpy as np

from sicer.lib import GenomeData


def get_bed_coords(chrom_reads, chrom_length, fragment_size, chrom, verbose):
    """
    *This takes into account the identical tags
    *Tags on different strands are positioned differently
        -> tag start (int(sline[1])) + fragment_size/2
        <- tag start (int(sline[2])) - 1 - fragment_size/2, the extra -1 is because that bed format has open-half, the sline[2] is not included.
    The stored positions are not the midpoint rather than the start
    The interface is no longer the same as that for getBedCoords(file)
    input:
        file:  the file that has the raw tag data from one chromosome
        fragment_size: the fragment size after CHIP experiment.
    output:
        return: a sorted list of positions which might have redundent entries
    """

    postive_tag_counts = 0
    negative_tag_counts = 0
    shift = int(round(fragment_size / 2))
    taglist = []
    print_return = ""
    for read in chrom_reads:
        chrom = read[0]
        start = read[1]
        end = read[2]
        name = read[3]
        score = read[4]
        strand = read[5]
        if (start < 0):
            if verbose:
                print_return += ("Ilegitimate read with start less than zero is ignored \n"
                                 + chrom + "\t" + str(start) + "\t" + str(
                            end) + "\t" + name + "\t" + str(score) + "\t" + strand + "\n")
        elif (end >= chrom_length):
            if verbose:
                print_return += (
                            "Ilegitimate read with end beyond chromosome length " + str(chrom_length) + " is ignored \n"
                            + chrom + "\t" + str(start) + "\t" + str(
                        end) + "\t" + name + "\t" + str(score) + "\t" + strand + "\n")
        else:
            if (strand == '+'):
                position = start + shift
                # If the position is beyond limit then don't shift.
                if (position >= chrom_length):
                    position = chrom_length - 1
                taglist.append(position)
                postive_tag_counts += 1

            elif (strand == '-'):
                position = end - 1 - shift
                # in case the shift move the positions
                # beyond zero, use zero
                if (position < 0):
                    position = 0  # UCSC genome coordinate is 0-based
                taglist.append(position)
                negative_tag_counts += 1
    taglist.sort()

    total_tag_counts = postive_tag_counts + negative_tag_counts
    print_return += 'Total count of ' + chrom + ' tags: ' + str(total_tag_counts)
    if verbose:
        print_return += ('  ('+str(postive_tag_counts) + ' positive tags, ' + str(negative_tag_counts) + ' negative tags)')

    return (taglist, print_return)


def Generate_windows_and_count_tags(taglist, chrom, chrom_length, window_size):
    """
    taglist: sorted list of positions that includes every tag on a chromosome
    window_size: the artificial bin size for binning the tags
    bed_vals: a dictionary keyed by the start of tag_containing
        windows, with value being the tag count in the window.

    In this function, the bins are set up using an absolute coordinate
    system.  Namely [0, window_size-1),[window_size,
    2*window_size-1). If the last window goes beyond the limit of the chromosome,
    that window is ignored.

    The result writen into the file is guaranteed to be already sorted
    within a chromosome.
    """
    chrom_graph = []
    total_tag_count = 0

    if (len(taglist) > 0):
        current_window_start = (taglist[0] // window_size) * window_size
        tag_count_in_current_window = 1

        for i in range(1, len(taglist)):
            start = (taglist[i] // window_size) * window_size

            if start == current_window_start:
                tag_count_in_current_window += 1

            elif start > current_window_start:
                # All the tags in the previous window have been counted
                current_window_end = current_window_start + window_size - 1
                # if the window goes beyond the chromsome limit, it is discarded.

                if current_window_end < chrom_length:
                    # write the window to file
                    output = (chrom, current_window_start, current_window_end, tag_count_in_current_window)
                    chrom_graph.append(output)
                    total_tag_count += tag_count_in_current_window

                current_window_start = start
                tag_count_in_current_window = 1

            else:
                sys.stderr.write("Error!")
                sys.exit(1)

        current_window_end = current_window_start + window_size - 1
        # if the window goes beyond the chromsome limit, it is discarded.
        if current_window_end < chrom_length:
            output = (chrom, current_window_start, current_window_end, tag_count_in_current_window)
            chrom_graph.append(output)
            total_tag_count += tag_count_in_current_window

    return (chrom_graph, total_tag_count)


def makeGraphFile(args, filtered, chrom, chrom_length):
    file = args.treatment_file.replace('.bed', '')  # removes the .bed extension

    bed_file_name = file + '_' + chrom   # name of the ChIP-seq reads
    if filtered:
        bed_file_name = bed_file_name + '_filtered.npy'
    else:
        bed_file_name = bed_file_name + '.npy'

    chrom_reads = np.load(bed_file_name, allow_pickle=True)

    tag_list, print_return = get_bed_coords(chrom_reads, chrom_length, args.fragment_size, chrom, args.verbose)
    
    chrom_graph, tag_count = Generate_windows_and_count_tags(tag_list, chrom, chrom_length, args.window_size)


    file_save_name = file + '_' + chrom
    if filtered:
        file_save_name += '_filtered_graph.npy'
    else:
        file_save_name += '_graph.npy'
    
    #graph_dtype = np.dtype([('chrom', 'U6'), ('start', np.int32), ('end', np.int32), ('count', np.int32)])
    np_chrom_graph = np.array(chrom_graph, dtype=object)
    np.save(file_save_name, np_chrom_graph)
    return (tag_count, print_return)


def main(args, pool, filtered=False):
    chroms = GenomeData.species_chroms[args.species]
    chrom_lengths = GenomeData.species_chrom_lengths[args.species]

    list_of_args = []
    for i, chrom in enumerate(chroms):
        if chrom in chrom_lengths.keys():
            chrom_length = chrom_lengths[chrom]
        else:
            print("Can not find the length of ", chrom)
        list_of_args.append((chrom, chrom_length))

    # Use multiprocessing to partition the gneome in windows and generate the summary files in parallel processes
    #pool = mp.Pool(processes=min(args.cpu, len(chroms)))
    makeGraphFile_partial = partial(makeGraphFile, args, filtered)
    makeGraphFile_result = pool.starmap(makeGraphFile_partial, list_of_args)
    #pool.close()

    total_tag_count = 0
    for result in makeGraphFile_result:
        total_tag_count += result[0]
        print(result[1])

    return (total_tag_count)
