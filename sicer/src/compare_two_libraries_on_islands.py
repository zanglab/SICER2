#!/usr/bin/env python

# Authors: Chongzhi Zang, Weiqun Peng

# Modified by: Jin Yong Yoo

import multiprocessing as mp
import os
from functools import partial
from math import *

import numpy as np
import scipy.stats

from sicer.lib import GenomeData
from sicer.lib import Utility
from sicer.lib import associate_tags_with_regions


def calc_pvalue(chip_read_count, control_read_count, scaling_factor, pseudo_count):
    """
    Currently using poisson distribution

    scaling_factor: the factor that accounts for the differences of control library and ChIP library. effective control read count
    is control_read_count * scaling factor
    pseudocount: when control_read_count is zero, replace zero with pseudocount to alleviate the impact of statistical fluctuation

    output: pvalue
    """
    if control_read_count > 0:
        average = control_read_count * scaling_factor
    else:
        average = pseudo_count * scaling_factor
    if chip_read_count > average:
        pvalue = scipy.stats.poisson.sf(chip_read_count, average)[()]
    else:
        pvalue = 1
    return pvalue


def fdr(pvalue_list):
    """
    Calculate the multiple testing corrected p-value using BH
    """
    fdr_list = []
    pvaluearray = scipy.array(pvalue_list)
    totalnumber = len(pvalue_list)
    pvaluerankarray = scipy.stats.rankdata(pvaluearray)
    for i in range(totalnumber):
        fdr_value = pvalue_list[i] * totalnumber / pvaluerankarray[i]
        if fdr_value > 1:
            fdr_value = 1
        fdr_list.append(fdr_value)
    return fdr_list


def associate_tags_count_to_regions(args, path_A, path_B, scaling_factor, chrom):
    island_list = np.load(chrom + '_union_output.npy', allow_pickle=True)
    island_start_list = []
    island_end_list = []

    for island in island_list:
        island_start_list.append(island[1])
        island_end_list.append(island[2])

    totalA = 0
    island_A_readcount_list = [0] * len(island_list)
    treatment_A_file_name = args.treatment_file[0].replace('.bed', '') + '_' + chrom + '.npy'
    treatment_A_reads = np.load(os.path.join(path_A, treatment_A_file_name), allow_pickle=True)
    for read in treatment_A_reads:
        position = associate_tags_with_regions.tag_position(read, args.fragment_size)
        index = associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position)
        if index >= 0:
            island_A_readcount_list[index] += 1
            totalA += 1

    totalB = 0
    island_B_readcount_list = [0] * len(island_list)
    treatment_B_file_name = args.treatment_file[1].replace('.bed', '') + '_' + chrom + '.npy'
    treatment_B_reads = np.load(os.path.join(path_B, treatment_B_file_name), allow_pickle=True)
    for read in treatment_B_reads:
        position = associate_tags_with_regions.tag_position(read, args.fragment_size)
        index = associate_tags_with_regions.find_readcount_on_islands(island_start_list, island_end_list, position)
        if index >= 0:
            island_B_readcount_list[index] += 1
            totalB += 1

        # Calculate the p value.
    pvalue_A_vs_B_list = []
    pvalue_B_vs_A_list = []
    for index in range(len(island_list)):
        island = island_list[index]
        Acount = island_A_readcount_list[index]
        Bcount = island_B_readcount_list[index]
        pvalue_A_vs_B = calc_pvalue(Acount, Bcount, scaling_factor, 1)
        pvalue_A_vs_B_list.append(pvalue_A_vs_B)
        pvalue_B_vs_A = calc_pvalue(Bcount, Acount, 1 / scaling_factor, 1)
        pvalue_B_vs_A_list.append(pvalue_B_vs_A)

    np_island_A_readcount_list = np.array(island_A_readcount_list)
    np_island_B_readcount_list = np.array(island_B_readcount_list)
    np.save(chrom + '_readcount_A.npy', np_island_A_readcount_list)
    np.save(chrom + '_readcount_B.npy', np_island_B_readcount_list)

    np_pvalue_A_vs_B_list = np.array(pvalue_A_vs_B_list)
    np_pvalue_B_vs_A_list = np.array(pvalue_B_vs_A_list)
    np.save(chrom + '_pvalue_AB.npy', np_pvalue_A_vs_B_list)
    np.save(chrom + '_pvalue_BA.npy', np_pvalue_B_vs_A_list)

    return (totalA, totalB)


def main(args, path_to_tempdir_1, path_to_tempdir_2, A_library_size, B_library_size, pool):
    chroms = GenomeData.species_chroms[args.species]

    print("Library size of ", args.treatment_file[0], ":  ", A_library_size)
    print("Library size of ", args.treatment_file[1], ":  ", B_library_size)

    library_scaling_factor = A_library_size * 1.0 / B_library_size  # A vs B

    #pool = mp.Pool(processes=min(args.cpu, len(chroms)))
    associate_tag_count_to_regions_partial = partial(associate_tags_count_to_regions, args, path_to_tempdir_1,
                                                     path_to_tempdir_2, library_scaling_factor)
    tag_counts = pool.map(associate_tag_count_to_regions_partial, chroms)
    #pool.close()

    total_read_count_A = 0  # Count of all the reads of library A that belong in islands
    total_read_count_B = 0
    for count in tag_counts:
        total_read_count_A += count[0]
        total_read_count_B += count[1]

    print("Total number of A reads on islands is: ", total_read_count_A)
    print("Total number of B reads on islands is: ", total_read_count_B)

    island_A_readcount = []
    island_B_readcount = []
    pvalue_A_vs_B_list = []
    pvalue_B_vs_A_list = []

    for chrom in chroms:
        # These numpy arrays are the arrays stored by the parallel processes
        # Goal is to combine them into 4 arrays
        A_readcount = np.load(chrom + '_readcount_A.npy', allow_pickle=True)
        B_readcount = np.load(chrom + '_readcount_B.npy', allow_pickle=True)
        pvalue_AB_list = np.load(chrom + '_pvalue_AB.npy', allow_pickle=True)
        pvalue_BA_list = np.load(chrom + '_pvalue_BA.npy', allow_pickle=True)

        for i in range(0, len(A_readcount)):
            island_A_readcount.append(A_readcount[i])
            island_B_readcount.append(B_readcount[i])
            pvalue_A_vs_B_list.append(pvalue_AB_list[i])
            pvalue_B_vs_A_list.append(pvalue_BA_list[i])

    # Calculate the FDR
    fdr_A_vs_B_list = fdr(pvalue_A_vs_B_list)
    fdr_B_vs_A_list = fdr(pvalue_B_vs_A_list)

    # Output the islands read counts, normalized read counts, fc, pvalue both ways
    scaling_factor = 1000000
    pseudo_count = 1
    outfile_name = (args.treatment_file[0].replace('.bed', '') + '-and-' + args.treatment_file[1].replace('.bed', '') +
                    '-W' + str(args.window_size))
    if (args.subcommand == "SICER"):
        outfile_name += ('-G' + str(args.gap_size) + '-summary')
    else:
        outfile_name += '-summary'
    outfile_path = os.path.join(args.output_directory, outfile_name)
    with open(outfile_path, 'w') as outfile:
        outline = (
                    '#chrom' + "\t" + 'start' + "\t" + 'end' + "\t" + "Readcount_A" + "\t" + 'Normalized_Readcount_A' + "\t" + 'ReadcountB' + "\t" + 'Normalized_Readcount_B'
                    + "\t" + "Fc_A_vs_B" + "\t" + "pvalue_A_vs_B" + "\t" + "FDR_A_vs_B" + "\t" + "Fc_B_vs_A" + "\t" + "pvalue_B_vs_A" + "\t" + "FDR_B_vs_A" + "\n")
        outfile.write(outline)
        j = 0;
        for chrom in chroms:
            island_list = np.load(chrom + '_union_output.npy', allow_pickle=True)
            complete_island_list = []
            for index in range(len(island_list)):
                island = island_list[index]
                Acount = island_A_readcount[j]
                Bcount = island_B_readcount[j]
                normalized_A = Acount / float(A_library_size) * scaling_factor
                normalized_B = Bcount / float(B_library_size) * scaling_factor
                fc_A_vs_B = ((Acount + pseudo_count) * 1.0 / (Bcount + pseudo_count)) / library_scaling_factor
                fc_B_vs_A = ((Bcount + pseudo_count) * 1.0 / (Acount + pseudo_count)) * library_scaling_factor
                outline = (island[0] + "\t" + str(island[1]) + "\t" + str(island[2]) + "\t" + str(Acount) + "\t" + str(
                    normalized_A) + "\t" + str(Bcount) + "\t" + str(normalized_B) +
                           "\t" + str(fc_A_vs_B) + "\t" + str(pvalue_A_vs_B_list[j]) + "\t" + str(
                            fdr_A_vs_B_list[j]) + "\t" + str(fc_B_vs_A) + "\t" + str(pvalue_B_vs_A_list[j]) +
                           "\t" + str(fdr_B_vs_A_list[j]) + "\n")
                outfile.write(outline)
                island = (island[0], island[1], island[2], Acount, normalized_A, Bcount, normalized_B, fc_A_vs_B,
                          pvalue_A_vs_B_list[j], fdr_A_vs_B_list[j], fc_B_vs_A, pvalue_B_vs_A_list[j],
                          fdr_B_vs_A_list[j])
                complete_island_list.append(island)
                j += 1
            np_island_list = np.array(complete_island_list, dtype=object)
            np.save(chrom + '_union_island_summary.npy', np_island_list)

    # Calculate the correlations using normalized read counts
    A_array = np.array(island_A_readcount, dtype=float)
    B_array = np.array(island_B_readcount, dtype=float)

    # Normalization to reads per million
    A_array = A_array / float(A_library_size * scaling_factor)
    B_array = B_array / float(B_library_size * scaling_factor)
    pearson = scipy.stats.pearsonr(A_array, B_array)
    print("Pearson's correlation is: ", pearson[0], " with p-value ", pearson[1])
    spearman = scipy.stats.spearmanr(A_array, B_array)
    print("Spearman's correlation is: ", spearman[0], " with p-value ", spearman[1])
