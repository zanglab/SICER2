# Authors: Chongzhi Zang, Weiqun Peng

# Modified by: Jin Yong Yoo
# Modified by: Mengxue Tian (2025)

import multiprocessing as mp
import os
import re
import sys
import subprocess
from functools import partial
import numpy as np

from sicer.lib import GenomeData

'''Filters redundant reads according to the cutoff value by taking a sorted list and comparing adjacent reads'''


def remove_redundant_1chrom_single_strand_sorted(reads, cutoff):
    current_start = 0
    current_end = 0
    current_count = 1
    total = 0
    retained = 0
    mask = []

    for i, read in enumerate(reads):
        total += 1
        start = read[1]
        end = read[2]
        if start != current_start:
            retained += 1
            current_start = start
            current_end = end
            current_count = 1
        elif end != current_end:
            retained += 1
            current_start = start
            current_end = end
            current_count = 1
        else:
            current_count += 1
            if current_count <= cutoff:
                retained += 1
            else:
                mask.append(i)

    return (total, retained, mask)


'''Separates reads by positive and negative strands before filtering redudant reads.
    Saves the filtered reads as a numpy binary file in temporary directory created in run_SICER.
    This is because python cannot pass extremely large objects between parallel processes'''


def strand_broken_remove(chrom, cutoff, file, chrom_reads, input_type="SE"):
    """
    Remove redundant reads on this chromosome.

    SE mode:
        - Split reads by strand (+ / -)
        - Sort and deduplicate separately
        - Print plus/minus counts

    PE mode:
        - Ignore strand
        - Mix all reads together
        - Deduplicate only by (start, end)
        - Print total/retained only
    """

    print_return = ""

    # ======== PE MODE ========
    if input_type == "PE":
        # Sort by start, then end
        sorted_reads = np.sort(chrom_reads, order=['start', 'end'])

        # Deduplicate by (start, end)
        total = 0
        retained = 0
        mask = []
        current_start = None
        current_end = None
        current_count = 0

        for i, read in enumerate(sorted_reads):
            total += 1
            start = read['start']
            end = read['end']

            if start != current_start or end != current_end:
                # New unique fragment
                retained += 1
                current_start = start
                current_end = end
                current_count = 1
            else:
                # Repeated fragment
                current_count += 1
                if current_count > cutoff:
                    mask.append(i)

        filtered_reads = np.delete(sorted_reads, obj=mask)

        np.save(f"{file}_{chrom}.npy", filtered_reads)

        print_return += (
            f"{chrom:<5s}"
            f"{total:^25d}"
            f"{retained:^25d}"
        )

        return (print_return, retained)

    # ======== SE MODE (original behavior) ========
    sorted_reads = np.sort(chrom_reads, order=['strand','start','end'])

    mark = 0
    for i in range(len(sorted_reads)):
        if sorted_reads[i][5] != '+':
            mark = i
            break

    (p_total, p_retained, p_mask) = remove_redundant_1chrom_single_strand_sorted(sorted_reads[:mark], cutoff)
    (m_total, m_retained, m_mask) = remove_redundant_1chrom_single_strand_sorted(sorted_reads[mark:], cutoff)
    m_mask = [m_mask[i]+mark for i in range(len(m_mask))]
    mask = p_mask + m_mask
    filtered_reads = np.delete(sorted_reads, obj=mask)

    print_return += ('{:<5s}{:^25d}{:^25d}{:^25d}{:^25d}'.format(
        chrom, p_total, p_retained, m_total, m_retained))

    np.save(f"{file}_{chrom}.npy", filtered_reads)
    total_retained = p_retained + m_retained

    return (print_return, total_retained)


'''Opens the given file and greps all the reads of the given chromsome and stores them in a list of tuples
    Tuple format: (chrom,start,end,name,score,strand)'''


def match_by_chrom(file, chrom):
    match = chrom + "[[:space:]]"
    matched_reads = subprocess.Popen(['grep', match, file], stdout=subprocess.PIPE) #Use Popen so that if no matches are found, it doesn't throw an exception
    chrom_reads = str(matched_reads.communicate()[0],'utf-8').splitlines()  # generates a list of each reads, which are represented by a string value
    file_name = os.path.basename(file)
    read_dtype = np.dtype([('chrom', 'U6'), ('start', np.int32), ('end', np.int32), ('name', 'U20'), ('score', np.int32), ('strand', 'U1')])
    processed_reads = np.empty(len(chrom_reads), dtype=read_dtype)

    for i, reads in enumerate(chrom_reads):
        reads = re.split('\t', reads)
        if (len(reads) < 6):
            sys.stderr.write(
                "Error: Input BED files must have the first six fields. Check " + file_name + " to see if it has the following fields: chrom, chromStart, chromEnd, name, score, and strand\n")
            sys.exit(1)
        reads[1] = int(reads[1])
        reads[2] = int(reads[2])
        processed_reads[i] = tuple(reads)

    return processed_reads


'''Function designed for handling multiprocessing. Separates all reads by chromosome
    and then filters redudant reads'''


def find_and_filter_reads(path_to_file, cutoff, input_type, chrom):
    file_name = os.path.basename(path_to_file)
    file_name = file_name.replace('.bed', '')

    chrom_reads = match_by_chrom(path_to_file, chrom)  # Separates all reads by chromosome
    return strand_broken_remove(chrom, cutoff, file_name, chrom_reads, input_type)



'''path_to_file: complete path to the .bed file that needs to processed for redudant reads'''


def main(args, path_to_file, pool):
    chroms = GenomeData.species_chroms[args.species]  # list of chromosomes of the given species
    cutoff = args.redundancy_threshold

    # Use multiprocessing module to run parallel processes for each chromosome
    # NOTE: we pass args.input_type down so SE/PE has the same removing redundant logic
    find_and_filter_reads_partial = partial(find_and_filter_reads, path_to_file, cutoff, args.input_type)
    filtered_result = pool.map(find_and_filter_reads_partial, chroms)
    # pool.close()

    total_read_count = 0

    # Print header
    if args.input_type == "PE":
        # Fragment-level (PE-like) mode: only total/retained fragments
        print('-' * 60)
        print('{:<5s}{:^25s}{:^25s}'.format(
            "chrom", "Total fragments", "Retained fragments"
        ))
        print('-' * 60)
    else:
        # SE mode: plus/minus reads separately
        print('-' * 105)
        print('{:<5s}{:^25s}{:^25s}{:^25s}{:^25s}'.format(
            "chrom",
            "Total plus reads",
            "Retained plus reads",
            "Total minus reads",
            "Retained minus reads"
        ))
        print('-' * 105)

    # Print per-chromosome summary from strand_broken_remove
    for result in filtered_result:
        print(result[0])
        total_read_count += result[1]

    return total_read_count
