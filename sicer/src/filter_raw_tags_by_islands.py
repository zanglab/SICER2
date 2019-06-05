#!/usr/bin/env python3

# Authors: Chongzhi Zang, Weiqun Peng

# Modified by: Jin Yong Yoo

import bisect
import multiprocessing as mp
import os
from functools import partial
from math import *

import numpy as np

from sicer.lib import GenomeData


def tag_position(read, fragment_size):
    shift = int(round(fragment_size / 2))
    if (read[5] == '+'):
        return read[1] + shift
    elif (read[5] == '-'):
        return read[2] - 1 - shift


def filter_tags_by_islands(file_name, fragment_size, chrom):
    island_list = np.load(file_name + '_' + chrom + '_island_summary.npy', allow_pickle=True)
    read_list = np.load(file_name + '_' + chrom + '.npy', allow_pickle=True)
    filtered_reads = []
    if (len(island_list) > 0):
        island_start_list = []
        island_end_list = []
        for island in island_list:
            island_start_list.append(island[1])
            island_end_list.append(island[2])

        island_start_list.sort()
        island_end_list.sort()

        for read in read_list:
            position = tag_position(read, fragment_size)
            if bisect.bisect_right(island_start_list, position) - bisect.bisect_left(island_end_list, position) == 1:
                filtered_reads.append(read)

    np_filtered_reads = np.array(filtered_reads, dtype=object)
    np.save(file_name + '_' + chrom + '_filtered.npy', np_filtered_reads)


def main(args, pool):
    chroms = GenomeData.species_chroms[args.species];
    treatment_file = args.treatment_file.replace('.bed', '')

    # Use multiprocessing to filter raw tags by islands in parallel processes
    #pool = mp.Pool(processes=min(args.cpu, len(chroms)))
    filter_tags_by_islands_partial = partial(filter_tags_by_islands, treatment_file, args.fragment_size)
    pool.map(filter_tags_by_islands_partial, chroms)
    #pool.close()

    output_file_name = treatment_file + '-W' + str(args.window_size)
    if (args.subcommand == "SICER"):
        output_file_name += '-G' + str(args.gap_size)
    output_file_name += '-FDR' + str(args.false_discovery_rate) + '-islandfiltered.bed'
    outfile_path = os.path.join(args.output_directory, output_file_name)
    with open(outfile_path, 'w') as outfile:
        for chrom in chroms:
            filtered_bed = np.load(treatment_file + '_' + chrom + '_filtered.npy', allow_pickle=True)
            for read in filtered_bed:
                output_line = str(read[0]) + '\t' + str(read[1]) + '\t' + str(read[2]) + '\t' + str(
                    read[3]) + '\t' + str(read[4]) + '\t' + str(read[5]) + '\n'
                outfile.write(output_line)

