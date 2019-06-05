#!/usr/bin/env python
# Author: Jin Yong Yoo

import multiprocessing as mp
import os
from math import *

import numpy as np

from sicer.lib import GenomeData


def get_counts(graph_file):
    chrom_graph = np.load(graph_file, allow_pickle=True)
    count = 0
    for line in chrom_graph:
        count += line[3]
    return count


def main(args, output_file_name, pool):
    chroms = GenomeData.species_chroms[args.species];
    scaling_factor = 1000000
    total_count = 0  # total count of islands
    file = args.treatment_file.replace('.bed', '')

    filtered_mode = "islandfiltered" in output_file_name
    list_of_graph_files = []
    for chrom in chroms:
        if filtered_mode:
            list_of_graph_files.append(file + '_' + chrom + '_filtered_graph.npy')
        else:
            list_of_graph_files.append(file + '_' + chrom + '_graph.npy')

    # Use multiprocessing to count the number of islands
    #pool = mp.Pool(processes=min(args.cpu, len(chroms)))
    count_results = pool.map(get_counts, list_of_graph_files)
    #pool.close()
    for count in count_results:
        total_count += count

    scaling_factor = total_count / 1000000 * (args.window_size / 1000)
    outfile_path = os.path.join(args.output_directory, output_file_name)

    # Normalize tag count using the scaling factor and generate a file in WIG format
    with open(outfile_path, 'w') as outfile:
        if filtered_mode:
            file = file + '-islandfiltered'
        outfile.write("track type=wiggle_0 name=" + file + "\n")
        for i in range(0, len(chroms)):
            chrom_graph = np.load(list_of_graph_files[i], allow_pickle=True)
            if (len(chrom_graph) > 0):
                outfile.write("variableStep chrom=" + chroms[i] + " span=" + str(args.window_size) + "\n")
                for window in chrom_graph:
                    normalized_tag_count = round(float(window[3] / scaling_factor), 2)
                    start_coord = window[1] + 1
                    outfile.write(str(start_coord) + "\t" + str(normalized_tag_count) + "\n")
