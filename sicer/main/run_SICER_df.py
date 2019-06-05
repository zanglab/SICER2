# Developed by Zang Lab at University of Virginia - 2018

#Author: Jin Yong Yoo

import copy
import os
import shutil
import sys
import tempfile
import multiprocessing as mp

curr_path = os.getcwd()

# From SICER Package
from sicer.lib import GenomeData
from sicer.main import run_SICER
from sicer.src import find_union_islands
from sicer.src import compare_two_libraries_on_islands
from sicer.src import filter_islands_by_significance


def main(args):
    # Checks if there is a control library
    control_lib_exists = True
    if (args.control_file is None):
        control_lib_exists = False

    # Create deep copy of the 'args' object for each treatment
    args_1 = copy.deepcopy(args)
    args_2 = copy.deepcopy(args)

    # Format each args for SICER run
    args_1.treatment_file = str(args.treatment_file[0])
    args_2.treatment_file = str(args.treatment_file[1])
    args_1.df = False
    args_2.df = False

    if (control_lib_exists):
        args_1.control_file = str(args.control_file[0])
        args_2.control_file = str(args.control_file[1])

        # Execute run_SICER for each treatment library
    temp_dir_1, library_size_file1 = run_SICER.main(args_1, True)
    temp_dir_2, library_size_file2 = run_SICER.main(args_2, True)


    try:
        temp_dir = tempfile.mkdtemp()
        # Change current working directory to temp_dir
        os.chdir(temp_dir)
    except:
        sys.exit(
            "Temporary directory required for SICER_df cannot be created. Check if directories can be created in %s."
            % curr_path)
    try:
        num_chroms = len(GenomeData.species_chroms[args.species])
        pool = mp.Pool(processes=min(args.cpu, num_chroms))

        # Find the union island between two treatment files. It will generate a summary file
        print("\n")
        args.treatment_file[0] = os.path.basename(args.treatment_file[0])
        args.treatment_file[1] = os.path.basename(args.treatment_file[1])
        print("Finding all the union islands of ", args.treatment_file[0], "and ", args.treatment_file[1], "...")
        find_union_islands.main(args, temp_dir_1, temp_dir_2, pool)
        print("\n")

        # Compare two treatment libraries
        print("Comparing two treatment libraries...")
        compare_two_libraries_on_islands.main(args, temp_dir_1, temp_dir_2, library_size_file1, library_size_file2, pool)
        print("\n")

        print("Identifying significantly increased islands using BH corrected p-value cutoff...")
        filter_islands_by_significance.main(args, 9, pool)
        print("\n")

        print("Identifying significantly decreased islands using BH-corrected p-value cutoff...")
        filter_islands_by_significance.main(args, 12, pool)
        print("\n")

        pool.close()
        pool.join()

    finally:
        print("Removing all temporary directories and all files in it.")
        shutil.rmtree(temp_dir)
        shutil.rmtree(temp_dir_1)
        shutil.rmtree(temp_dir_2)
