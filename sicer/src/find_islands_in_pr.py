# Authors: Chongzhi Zang, Weiqun Peng, Dustin E Schones and Keji Zhao

# Modified by: Jin Yong Yoo

import multiprocessing as mp
import os
from functools import partial
from math import *

import numpy as np

from sicer.lib import Background_island_probscore_statistics
from sicer.lib import GenomeData

"""
Take in coords for bed_gaph type summary files and find 'islands' of modifications.
There are a number options here that can be turned on or off depending on need
Right now:

(1) scan for all 'islands' where the single window or consecutive
windows

(2) remove all single window 'islands' -- this can be commented out
when looking for more localized signals (such as TF binding sites?)

(3) try to combine islands that are within gap distance of each other.
This gap distance is supplemented by a window_buffer just so we don't
do anything stupid with window sizes

(4) Remove all single window combined islands -- if step (2) was done,
this is redundant

(5) Lastly, filter out all the islands we've found that have a total
score < islands_minimum_tags
"""


# Factorial
def fact(m):
    value = 1.0;
    if m != 0:
        while m != 1:
            value = value * m;
            m = m - 1;
    return value;


# Return the log of a factorial, using Srinivasa Ramanujan's approximation when m>=20
def factln(m):
    if m < 20:
        value = 1.0;
        if m != 0:
            while m != 1:
                value = value * m;
                m = m - 1;
        return log(value);
    else:
        return m * log(m) - m + log(m * (1 + 4 * m * (1 + 2 * m))) / 6.0 + log(pi) / 2;


def poisson(i, average):
    if i < 20:
        return exp(-average) * average ** i / fact(i);
    else:
        exponent = -average + i * log(average) - factln(i);
        return exp(exponent);


def combine_proximal_islands(islands, gap, window_size_buffer=3):
    """
    islands: a list of tuples of following format: (chrom, start, end, score)
    Therefore, "islands[index][1]" would mean the start position of the window at the given index
    Extend the regions found in the find_continuous_region function.
    If gap is not allowed, gap = 0, if one window is allowed, gap = window_size (200)

    Return a list of combined regions.
    """

    proximal_island_dist = gap + window_size_buffer;

    if len(islands) == 0:
        return []
    final_islands = []
    current_island = islands[0];

    if len(islands) == 1:
        final_islands = islands;
    else:
        for index in range(1, len(islands)):
            dist = islands[index][1] - current_island[2];
            if dist <= proximal_island_dist:
                current_island[2] = islands[index][2];
                current_island[3] += islands[index][3];
            else:
                final_islands.append(current_island);
                current_island = islands[index];
        # The last island:
        final_islands.append(current_island);

    return final_islands;


def find_region_above_threshold(island_list, score_threshold):
    filtered_islands = [];
    for island in island_list:
        if island[3] >= (score_threshold - .0000000001):
            filtered_islands.append(island);
    return filtered_islands;


def filter_ineligible_windows(chrom_graph, min_tags_in_window, average):
    '''Filters windows that have tag count lower than the minimum threshold count and calculates score for windows that meet the minimum count.
        Score is defined as s = -log(Poisson(read_count,lambda))'''

    filtered_chrom_graph = []
    for window in chrom_graph:
        read_count = window[3]
        score = -1
        if (read_count >= min_tags_in_window):
            prob = poisson(read_count, average);
            if prob < 1e-250:
                score = 1000;  # outside of the scale, take an arbitrary number.
            else:
                score = -log(prob)
        eligible_window = (window[0], window[1], window[2], score)
        if score > 0:
            filtered_chrom_graph.append(eligible_window);

    np_filtered_chrom_graph = np.array(filtered_chrom_graph, dtype=object)
    return np_filtered_chrom_graph


def filter_and_find_islands(min_tags_in_window, gap_size, score_threshold, average, verbose, graph_file):
    '''Function for handling multiprocessing. Calls functions for filtering windows and finding islands.'''
    number_of_islands = 0
    print_return = ""
    chrom_graph = np.load(graph_file, allow_pickle=True)

    if (len(chrom_graph) > 0):
        chrom = chrom_graph[0][0]
        filtered_chrom_graph = filter_ineligible_windows(chrom_graph, min_tags_in_window, average)
        islands = combine_proximal_islands(filtered_chrom_graph, gap_size, 2);
        islands = find_region_above_threshold(islands, score_threshold);
        number_of_islands += len(islands)
        if not (len(islands) > 0):
            if verbose:
                print_return += chrom + " does not have any islands meeting the required significance"
        np.save(graph_file, islands)

    return (graph_file, number_of_islands, print_return)


def main(args, total_read_count, pool):
    print("Species: ", args.species);
    print("Window_size: ", args.window_size);
    print("Gap size: ", args.gap_size);
    print("E value is:", args.e_value);
    print("Total read count:", total_read_count)
    chroms = GenomeData.species_chroms[
        args.species];  # list of chromsomes for the given species (e.g. chr1, chr2, ... , chrx)
    genome_length = sum(GenomeData.species_chrom_lengths[args.species].values());  # list of length of each chromsomes
    effective_genome_length = int(args.effective_genome_fraction * genome_length);
    average = float(total_read_count) * args.window_size / effective_genome_length;  # average read count

    print("Genome Length: ", genome_length);
    print("Effective genome Length: ", effective_genome_length);
    print("Window average:", average);

    window_pvalue = 0.20;
    bin_size = 0.001;
    background = Background_island_probscore_statistics.Background_island_probscore_statistics(total_read_count,
                                                                                               args.window_size,
                                                                                               args.gap_size,
                                                                                               window_pvalue,
                                                                                               effective_genome_length,
                                                                                               bin_size);
    min_tags_in_window = background.min_tags_in_window

    print("Window pvalue:", window_pvalue)
    print("Minimum num of tags in a qualified window: ", min_tags_in_window)  # first threshold cutoff

    print("\nDetermining the score threshold from random background...");
    # determine threshold from random background
    score_threshold = background.find_island_threshold(args.e_value);
    print("The score threshold is:", score_threshold);

    # generate the probscore summary graph file, only care about enrichment
    # filter the summary graph to get rid of windows whose scores are less than window_score_threshold

    file = args.treatment_file.replace('.bed', '')
    list_of_graph_files = []
    for chrom in chroms:
        list_of_graph_files.append(file + '_' + chrom + '_graph.npy')

    # Use multiprocessing to filter windows with tag count below minimum requirement
    print(
        "Generating the enriched probscore summary graph and filtering the summary graph to eliminate ineligible windows... ");
    #pool = mp.Pool(processes=min(args.cpu, len(chroms)))
    filter_and_find_islands_partial = partial(filter_and_find_islands, min_tags_in_window, args.gap_size,
                                              score_threshold, average, args.verbose)
    filtered_islands_result = pool.map(filter_and_find_islands_partial, list_of_graph_files)
    #pool.close()

    file_name = args.treatment_file.replace('.bed', '')
    outfile_path = os.path.join(args.output_directory, (file_name + '-W' + str(args.window_size)
                                                        + '-G' + str(args.gap_size) + '.scoreisland'))
    total_number_islands = 0
    path_to_filtered_graph = []
    with open(outfile_path, 'w') as outfile:
        for i in range(0, len(filtered_islands_result)):
            filtered_chrom_graph = np.load(filtered_islands_result[i][0],allow_pickle=True)
            path_to_filtered_graph.append(filtered_islands_result[i][0])
            total_number_islands += filtered_islands_result[i][1]
            if (filtered_islands_result[i][2] != ""):
                print(filtered_islands_result[i][2])
            for window in filtered_chrom_graph:
                chrom = window[0]
                line = (chrom + '\t' + str(window[1]) + '\t' + str(window[2])
                        + '\t' + str(window[3]) + '\n')
                outfile.write(line)

    print("Total number of islands: ", total_number_islands);
