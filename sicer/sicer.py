# Developed at Zang Lab, University of Virginia (2019)
# Author: Jin Yong Yoo

# Python Imports
import os
import shutil
import sys
import tempfile

# SICER Internal Imports
from sicer.background_stat import BackgroundStatistics
from sicer.shared.genome_data import GenomeData
from sicer.bed_reader import BEDReader
from sicer.generate_windows import generate_windows
from sicer.find_islands import find_islands
from sicer.associate_tags_with_control import associate_tags_with_control
from sicer.utility.file_writers import WigFileWriter

WINDOW_PVALUE = 0.20
BIN_SIZE = 0.001

def run_sicer(args, df_run=False): 

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

    print("Calculating background statistics...")
    background_stat = BackgroundStatistics(
                        windows.getTotalTagCount(), 
                        args.window_size,
                        args.gap_size,
                        WINDOW_PVALUE,
                        effective_genome_length,
                        BIN_SIZE
                    )

    min_tag_threshold = background_stat.min_tags_in_window
    score_threshold = background_stat.find_island_threshold(args.e_value);

    print("Minimum number of tags in a qualified window:", min_tag_threshold)
    print("Score threshold:", score_threshold);

    islands = find_islands(windows, genome_data, min_tag_threshold, score_threshold, 
                            args.gap_size, avg_tag_count, args.cpu)

    if args.control_file is not None:
        genome_size = args.effective_genome_fraction * genome_length
        scaling_factor = treatment_reads.getReadCount() / control_reads.getReadCount()

        islands = associate_tags_with_control(islands, treatment_reads, control_reads,
                                                genome_size, scaling_factor,
                                                args.fragment_size, args.cpu
                                            )






    