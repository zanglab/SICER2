#!/usr/bin/env python
#Authors: Chongzhi Zang, Weiqun Peng

# Modified by: Jin Yong Yoo

import multiprocessing as mp
import os
from functools import partial
from math import *

import numpy as np

from sicer.lib import GenomeData


# Function designed for handling multiprocessing. Executes the redundancy removal algorithm
# for each independent chromosome
def find_union_islands(no_control, temp_dir_1, temp_dir_2, chrom):
    file_name_1 = temp_dir_1 + '_' + chrom
    file_name_2 = temp_dir_2 + '_' + chrom
    if (no_control == True):
        file_name_1 += '_graph.npy'
        file_name_2 += '_graph.npy'
    else:
        file_name_1 += '_island_summary.npy'
        file_name_2 += '_island_summary.npy'

    island_list_1 = np.load(file_name_1, allow_pickle=True)
    island_list_2 = np.load(file_name_2, allow_pickle=True)
    if (len(island_list_1) == 0):
        island_list = island_list_2
    elif (len(island_list_2) == 0):
        island_list = island_list_1
    else:
        island_list = np.concatenate((island_list_1, island_list_2))

    union_island_list = []
    if (len(island_list) > 0):
        island_list = island_list[np.argsort(island_list[:, 1])]
        current = island_list[0]
        i = 1
        while i < len(island_list):
            compare = island_list[i]
            assert current[1] <= compare[2]
            if compare[1] > current[2]:
                union_island_list.append(current)
                current = compare
                i += 1
            else:
                current[2] = max(current[2], compare[2])
                i += 1
        union_island_list.append(current)
    np_union_island_list = np.array(union_island_list, dtype=object)
    np.save(chrom + '_union_output.npy', np_union_island_list)


def main(args, temp_dir_1, temp_dir_2, pool):
    chroms = GenomeData.species_chroms[args.species]

    # Partially fill out the full directory of the files we want to access
    temp_dir_1 += '/' + args.treatment_file[0].replace('.bed', '')
    temp_dir_2 += '/' + args.treatment_file[1].replace('.bed', '')

    no_control = args.control_file is None

    # Use multiprocessing module to run parallel processes for each chromosome
    #pool = mp.Pool(processes=min(args.cpu, len(chroms)))
    find_union_islands_partial = partial(find_union_islands, no_control, temp_dir_1, temp_dir_2)
    pool.map(find_union_islands_partial, chroms)
    #pool.close()

    outfile_name = (args.treatment_file[0].replace('.bed', '') + '-vs-' + args.treatment_file[1].replace('.bed', '') + '-W' + str(
        args.window_size))
    if (args.subcommand == "SICER"):
        outfile_name += '-G' + str(args.gap_size) + '-E' + str(args.e_value) + '-union.island'
    else:
        outfile_name += '-union.island'
    outfile_path = os.path.join(args.output_directory, outfile_name)

    with open(outfile_path, 'w') as outfile:
        for chrom in chroms:
            union_island_list = np.load(chrom + '_union_output.npy', allow_pickle=True)
            for island in union_island_list:
                output_line = island[0] + '\t' + str(island[1]) + '\t' + str(island[2]) + '\n'
                outfile.write(output_line)

