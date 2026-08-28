# Authors: Dustin E Schones, Chongzhi Zang, Weiqun Peng and Keji Zhao

# Modified by: Jin Yong Yoo
# Modified by: Mengxue Tian (2025)

import multiprocessing as mp
from functools import partial
from math import *
import sys
import numpy as np

from sicer.lib import GenomeData







def get_bed_coords(chrom_reads, chrom_length, fragment_size, chrom, verbose, input_type):
    """
    Convert BED reads on a chromosome into a list of 1D tag positions.

    Two modes are supported:
      - SE (single-end): use a shift model based on fragment_size/2 and strand
         -> '+' strand: position = start + fragment_size/2
         -> '-' strand: position = end - 1 - fragment_size/2
        The extra -1 is because BED end is open (0-based, end not included).

      - PE (paired-end): assume each BED entry already represents a fragment
        (start, end). The tag position is the fragment midpoint:
         -> position = (start + end) // 2
        Strand is ignored in this mode.

    Returns:
        taglist: sorted list of integer positions (one per tag)
        print_return: summary string for logging
    """
       
    # ---------- Sanity checks & auto-guard for wrong SE usage ----------
    # We only do these heuristics in SE mode; PE mode pass.
    if input_type == "SE":
        # calculate fragment length distribution
        frag_lengths = chrom_reads['end'] - chrom_reads['start']
        if len(frag_lengths) > 0:
            median_len = np.median(frag_lengths)
            strands = chrom_reads['strand']
            unique_strands = set(strands)

            # “missing strand”：only have '.' or nothing, no '+' / '-'
            missing_strand = unique_strands.issubset({'.', ''})

            # Case 1: looks like raw SE read（short reads），but no strand → give error
            if missing_strand and median_len < 100:
                msg = (
                    f"ERROR: {chrom}: Input is in SE mode but appears to be raw read-level BED "
                    f"without strand information (median length ~{median_len:.1f} bp).\n"
                    "SE mode requires proper +/- strand to shift reads by fragment_size/2.\n"
                    "Please provide 6-column BED with strand, or convert to fragment-level BED "
                    "and run with --input_type PE.\n"
                )
                sys.stderr.write(msg)
                raise ValueError(
                    "SE read-level BED is missing strand information. "
                    "Provide strand or use fragment-level BED with --input_type PE."
                )

            # Case 2: no strand，but fragment is long → more like fragment BED（similar with PE）
            # use as "SE", and give warning
            if missing_strand and median_len > 150:
                msg = (
                    f"WARNING: {chrom}: SE mode was requested, but input BED has no strand and "
                    f"appears to be fragment-level (median length ~{median_len:.1f} bp).\n"
                    "Consider rerunning with --input_type PE if this is truly PE/fragment data.\n"
                )
                sys.stderr.write(msg)

            # Case 3: have strand，but all '+'，and length is long → highly possible to be PE fragment but input-type is SE
            if input_type == "SE" and unique_strands == {'+'}:
                p90 = np.percentile(frag_lengths, 90)
                if p90 > (fragment_size * 1.5):
                    msg = (
                        f"WARNING: {chrom}: Detected all '+' strands and large fragment lengths "
                        f"(90th percentile ~{p90:.1f} bp, fragment_size={fragment_size}).\n"
                        "It looks like PE fragment BED, but input_type=SE.\n"
                        "Please rerun with --input_type PE to avoid incorrect shifting.\n"
                    )
                    sys.stderr.write(msg)
                    raise ValueError(
                        "Input BED appears to be PE fragment data but --input_type=SE was used.\n"
                        "Please rerun with --input_type PE."
                    )


    positive_tag_counts = 0
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

        # Basic coordinate sanity checks
        if start < 0:
            if verbose:
                print_return += (
                    "Illegitimate read with start < 0 is ignored:\n"
                    f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"
                )
            continue

        if end >= chrom_length:
            if verbose:
                print_return += (
                    "Illegitimate read with end beyond chromosome length "
                    f"{chrom_length} is ignored:\n"
                    f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n"
                )
            continue

        # Compute tag position depending on input_type
        if input_type == "PE":
            # Paired-end: use fragment midpoint, ignore strand
            position = (start + end) // 2
            if position < 0:
                position = 0
            if position >= chrom_length:
                position = chrom_length - 1
            taglist.append(position)
            positive_tag_counts += 1  # count all tags as "positive" for reporting

        else:  # SE mode: original SICER behavior (shift by fragment_size/2)
            if strand == '+':
                position = start + shift
                # If the position is beyond limit then don't shift.
                if position >= chrom_length:
                    position = chrom_length - 1
                taglist.append(position)
                positive_tag_counts += 1

            elif strand == '-':
                position = end - 1 - shift
                # If the shift moves the position below 0, clamp to 0
                if position < 0:
                    position = 0  # UCSC genome coordinate is 0-based
                taglist.append(position)
                negative_tag_counts += 1

            else:
                # Unknown strand; in SE mode we can ignore or treat as '+'
                # Here we ignore silently.
                continue

    taglist.sort()

    total_tag_counts = positive_tag_counts + negative_tag_counts
    print_return += 'Total count of ' + chrom + ' tags: ' + str(total_tag_counts)
    if verbose:
        print_return += (
            '  (' + str(positive_tag_counts) + ' positive tags, '
            + str(negative_tag_counts) + ' negative tags)'
        )

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

    tag_list, print_return = get_bed_coords(
    chrom_reads,
    chrom_length,
    args.fragment_size,
    chrom,
    args.verbose,
    args.input_type  # new argument
    )
    
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

