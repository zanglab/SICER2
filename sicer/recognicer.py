# Developed at Zang Lab, University of Virginia (2019)
# Author: Jin Yong Yoo

# Python Imports
import os
from copy import deepcopy

# SICER Internal Imports
from sicer.background_stat import BackgroundStatistics
from sicer.shared.genome_data import GenomeData
from sicer.bed_reader import BEDReader
from sicer.generate_windows import generate_windows
from sicer.coarsegraining import find_islands_by_coarsegraining
from sicer.associate_tags_with_control import associate_tags_with_control
from sicer.filter_islands_by_fdr import filter_islands_by_fdr
from sicer.recover_significant_reads import recover_significant_reads
from sicer.find_union_islands import find_union_islands
from sicer.compare_two_libraries import compare_two_libraries
from sicer.file_writers import WigFileWriter, IslandFileWriter, BEDFileWriter, DiffExprIslandFileWriter

def run_recognicer(args, df_run=False): 

    genome_data = GenomeData(args.species)
    base_name = os.path.splitext(os.path.basename(args.treatment_file))[0]

    treatment_reader = BEDReader(args.treatment_file, genome_data, args.cpu, args.redundancy_threshold)
    treatment_reads = treatment_reader.read_file()

    if args.control_file is not None:
        control_reader = BEDReader(args.control_file, genome_data, args.cpu, args.redundancy_threshold)
        control_reads = control_reader.read_file()

    windows = generate_windows(treatment_reads, genome_data, args.fragment_size, args.window_size, args.cpu)
    WigFileWriter(base_name, args.output_directory, windows, args.window_size, False).write()

    genome_length = sum(genome_data.chrom_length.values())
    effective_genome_length = int(args.effective_genome_fraction * genome_length)
    avg_tag_count = windows.getTotalTagCount() * args.window_size / effective_genome_length

    min_tag_threshold = int(avg_tag_count) + 1

    print("Minimum number of tags in a qualified window:", min_tag_threshold)

    islands = find_islands_by_coarsegraining(windows, genome_data, min_tag_threshold, args.window_size, 
                            args.step_size, args.step_score, args.cpu)

    IslandFileWriter(base_name, args.output_directory, "cgisland",
                                islands, args.window_size).write()

    if args.control_file is not None:
        genome_size = args.effective_genome_fraction * genome_length
        scaling_factor = treatment_reads.getReadCount() / control_reads.getReadCount()

        islands = associate_tags_with_control(islands, treatment_reads, control_reads,
                                                genome_size, scaling_factor,
                                                args.fragment_size, args.cpu
                                            )

        IslandFileWriter(base_name, args.output_directory, "summary",
                                islands, args.window_size).write()

        filtered_islands = filter_islands_by_fdr(islands, args.false_discovery_rate, args.cpu, False)

        IslandFileWriter(base_name, args.output_directory, "fdr-filtered", filtered_islands,
                                args.window_size, None, args.false_discovery_rate).write()

    if args.significant_reads:
        sig_reads = recover_significant_reads(filtered_islands, treatment_reads, args.fragment_size, args.cpu)
        BEDFileWriter(base_name, args.output_directory, sig_reads,args.window_size,
                        args.false_discovery_rate).write()

        sig_windows = generate_windows(sig_reads, genome_data, args.fragment_size, args.window_size, args.cpu)
        WigFileWriter(base_name, args.output_directory, sig_windows, 
                        args.window_size, True, args.false_discovery_rate).write()

    if df_run:
        return treatment_reads, islands

def run_recognicer_df(args):
    genome_data = GenomeData(args.species)

    # Create deep copy of the 'args' object for each treatment
    args_1 = deepcopy(args)
    args_2 = deepcopy(args)

    # Format each args for SICER run
    args_1.treatment_file = str(args.treatment_file[0])
    args_2.treatment_file = str(args.treatment_file[1])

    if args.control_file:
        args_1.control_file = str(args.control_file[0])
        args_2.control_file = str(args.control_file[1])

    # Execute SICER for each treatment
    treatment_reads_1, islands_1 = run_recognicer(args_1, True)
    treatment_reads_2, islands_2 = run_recognicer(args_2, True)

    file_name_1 = os.path.splitext(os.path.basename(args.treatment_file[0]))[0]
    file_name_2 = os.path.splitext(os.path.basename(args.treatment_file[1]))[0]

    union_islands = find_union_islands(islands_1, islands_2, genome_data, args.cpu)

    writer = DiffExprIslandFileWriter(
        file_name_1, file_name_2, args.output_directory, "union-island",
        union_islands, args.window_size
    )

    writer.write()

    df_islands = compare_two_libraries(
        genome_data, treatment_reads_1, treatment_reads_2, 
        union_islands, args.fragment_size, args.cpu
    )

    writer.islands = df_islands
    writer.file_type = "summary"
    writer.write()

    a_vs_b_filtered = filter_islands_by_fdr(
        df_islands, args.false_discovery_rate_df, args.cpu, True, True
    )
    
    writer.islands = a_vs_b_filtered
    writer.file_type = "fdr-filtered-increased"
    writer.fdr = args.false_discovery_rate_df
    writer.write()

    b_vs_a_filtered = filter_islands_by_fdr(
        df_islands, args.false_discovery_rate_df, args.cpu, True, False
    )

    writer.islands = b_vs_a_filtered
    writer.file_type = "fdr-filtered-decreased"
    writer.write()
