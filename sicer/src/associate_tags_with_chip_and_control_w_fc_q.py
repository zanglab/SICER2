# Authors: Chongzhi Zang, Weiqun Peng

# Modified by: Jin Yong Yoo

import multiprocessing as mp
import os
from functools import partial
from math import *

import numpy as np
import scipy
import scipy.stats

from sicer.lib import GenomeData;
from sicer.lib import associate_tags_with_regions


def associate_tag_count_to_regions(args, scaling_factor, control_library_size, genomesize, chrom):
    island_file = args.treatment_file.replace('.bed', '') + '_' + chrom + '_graph.npy'
    treatment_file = args.treatment_file.replace('.bed', '') + '_' + chrom + '.npy'
    control_file = args.control_file.replace('.bed', '') + '_' + chrom + '.npy'

    island_list = np.load(island_file, allow_pickle=True)
    island_start_list = [island[1] for island in island_list]
    island_end_list = [island[2] for island in island_list]
    
    pvalue_array = np.empty(len(island_list), dtype=np.float64)

    total_chip_count = 0
    island_chip_readcount_list = [0] * len(island_list)
    treatment_reads = np.load(treatment_file, allow_pickle=True)
    for read in treatment_reads:
        position = associate_tags_with_regions.tag_position(read, args.fragment_size)
        index = associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position)
        if index >= 0:  # if the read is found in a valid island
            island_chip_readcount_list[index] += 1
            total_chip_count += 1

    total_control_count = 0
    island_control_readcount_list = [0] * len(island_list)
    control_reads = np.load(control_file, allow_pickle=True)
    for read in control_reads:
        position = associate_tags_with_regions.tag_position(read, args.fragment_size)
        index = associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position)
        if index >= 0:
            island_control_readcount_list[index] += 1
            total_control_count += 1

    summary_list = []
    #pvalue_list = []
    for index in range(0, len(island_list)):
        island = island_list[index]
        observation_count = island_chip_readcount_list[index]
        control_count = island_control_readcount_list[index]
        average = 0
        if (control_count > 0):
            average = control_count * scaling_factor
        else:
            length = island[2] - island[1] + 1
            average = length * control_library_size * 1.0 / genomesize
            average = min(0.25, average) * scaling_factor;
        fc = float(observation_count) / float(average)
        if (observation_count > average):
            pvalue = scipy.stats.poisson.sf(observation_count, average)
        else:
            pvalue = 1

        pvalue_array[index] = pvalue
        output_line = (island[0], island[1], island[2], observation_count, control_count, pvalue, fc, 0.0)
        summary_list.append(output_line)

    np_output_lines = np.array(summary_list, dtype=object)
    file_name = args.treatment_file.replace('.bed', '') + '_' + chrom + '_' + 'island_summary.npy'
    np.save(file_name, np_output_lines)
    pvalue_save_name = chrom + '_pvalue.npy'
    np.save(pvalue_save_name, pvalue_array)
    return pvalue_save_name


def main(args, chip_library_size, control_library_size, pool):
    chroms = GenomeData.species_chroms[args.species];
    genomesize = sum(GenomeData.species_chrom_lengths[args.species].values());
    genomesize = args.effective_genome_fraction * genomesize;

    print("ChIP library read count:", chip_library_size)
    print("Control library read count:", control_library_size)

    totalchip = 0;
    totalcontrol = 0;
    scaling_factor = chip_library_size * 1.0 / control_library_size

    # Use multiprocessing to associate each read with an island
    #pool = mp.Pool(processes=min(args.cpu, len(chroms)))
    associate_tag_count_to_regions_partial = partial(associate_tag_count_to_regions, args, scaling_factor,
                                                     control_library_size, genomesize)
    p_value_files = pool.map(associate_tag_count_to_regions_partial, chroms)
    #pool.close()

    # Get the list of p-value from each parallel processes and concatenate them into one list of all p-values
    p_value_list = np.array([])
    for p_value_file in p_value_files:
        chrom_p_value_list = np.load(p_value_file, allow_pickle=True)
        p_value_list = np.concatenate([p_value_list, chrom_p_value_list])
        os.remove(p_value_file)
    p_value_rank_array = scipy.stats.rankdata(p_value_list)
    total_num_of_pvalue = len(p_value_list)

    index = 0
    file_name = args.treatment_file.replace('.bed', '')
    output_file_name = file_name + '-W' + str(args.window_size)
    if (args.subcommand == "SICER"):
        output_file_name += '-G' + str(args.gap_size) + '-islands-summary'
    elif (args.subcommand == "RECOGNICER"):
        output_file_name += '-islands-summary'
    outfile_path = os.path.join(args.output_directory, output_file_name)
    with open(outfile_path, 'w') as outfile:
        for chrom in chroms:
            island_file_name = file_name + '_' + chrom + '_' + 'island_summary.npy'
            island = np.load(island_file_name, allow_pickle=True)
            #modified_island = []
            for i in range(len(island)):
                line = island[i]
                totalchip += int(line[3])
                totalcontrol += int(line[4])
                alpha_stat = p_value_list[index] * total_num_of_pvalue / p_value_rank_array[index];
                if alpha_stat > 1:
                    alpha_stat = 1;
                    
                island[i][7] = alpha_stat
                outputline = (line[0] + '\t' + str(line[1]) + '\t' + str(line[2]) + '\t' + str(line[3]) + '\t' + str(
                    line[4]) + '\t' +
                              str(line[5]) + '\t' + str(line[6]) + '\t' + str(alpha_stat)+ '\n')
                outfile.write(outputline)
                #modified_island.append(tuple(line))
                index += 1

            np.save(island_file_name, island)

    print("Total number of chip reads on islands is:", totalchip)
    print("Total number of control reads on islands is:", totalcontrol)
