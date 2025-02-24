# Developed by Zang Lab at University of Virginia - 2018 & 2024

#Author: Jin Yong Yoo, Shiyu Jiang

import os
import shutil
import sys
import tempfile
import multiprocessing as mp
import logging
import time

curr_path = os.getcwd()

# From SICER Package
from sicer.lib import GenomeData
from sicer.src import remove_redundant_reads
from sicer.src import run_make_graph_file_by_chrom
from sicer.src import find_islands_in_pr
from sicer.src import associate_tags_with_chip_and_control_w_fc_q
from sicer.src import filter_islands_by_significance
from sicer.src import make_normalized_wig
from sicer.src import filter_raw_tags_by_islands
from sicer.src import clipper_fdr_control

''' args: ArgumentParser object formed form command line parameters
    df_run: If df_run is true, then this instance of SICER is called by SICER-DF module.
            Default value is False.'''


def main(args, df_run=False):
    logger = logging.getLogger("SICER 2")
    s_logger = logging.getLogger("s_logger")

    # Checks if there is a control library
    control_lib_exists = True
    if args.control_file is None:
        control_lib_exists = False
        logger.info("No control library is provided. SICER 2 will run without a control library.\n")
    else:
        logger.info("Control library is provided. SICER 2 will run with a control library.\n")

    # Creates temporary directory to contain all intermediate files. # TODO: define none
    if args.temp_directory is None:
        try:
            temp_dir = tempfile.mkdtemp()
            # Change current working directory to temp_dir
            os.chdir(temp_dir)
            logger.info(f"Temporary directory is set to {temp_dir}.")
        except:
            logger.error("Temporary directory required for SICER 2 cannot be created. "
                         "Check if directories can be created in %s. Error: %s" % (curr_path, str(e)))
            sys.exit(-1)
    else:
        if df_run:
            temp_dir = os.path.join(args.temp_directory, str(time.time()))
        else:
            temp_dir = args.temp_directory
        try:
            # check if the temp directory is empty
            if os.path.exists(temp_dir):
                if len(os.listdir(temp_dir)) > 0:
                    logger.error("Temporary directory required for SICER 2 is not empty. "
                                 "Remove this directory for temp storage.")
                    sys.exit(-1)
                    # shutil.rmtree(temp_dir)
            os.makedirs(temp_dir)
            os.chdir(temp_dir)
            logger.info(f'Temporary directory is set to {temp_dir}.')
        except Exception as e:
            logger.error("Temporary directory required for SICER 2 cannot be created. "
                         "Check if directories can be created in %s. Error: %s" % (curr_path, str(e)))
            sys.exit(-1)

    try:
        # Step 0: create Pool object for parallel-Processing
        num_chroms = len(GenomeData.species_chroms[args.species])
        pool = mp.Pool(processes=min(args.cpu, num_chroms))

        # Step 1: Remove redundancy reads in input file according to input threshold
        # Output is the total number of reads retained. Represents size of library.
        treatment_file_name = os.path.basename(args.treatment_file)
        logger.run(f"Preprocess the {treatment_file_name} file "
                   f"to remove redundancy with threshold of {args.redundancy_threshold}.")
        total_treatment_read_count = remove_redundant_reads.main(args, args.treatment_file, pool)
        args.treatment_file = treatment_file_name

        # Step 2: Remove redundancy reads in control library according to input threshold
        if (control_lib_exists):
            control_file_name = os.path.basename(args.control_file)
            logger.run(f"Preprocess the {control_file_name} file to remove redundancy with threshold of "
                       f"{args.redundancy_threshold}.")
            total_control_read_count = remove_redundant_reads.main(args, args.control_file, pool)
            args.control_file = control_file_name

        # Step 3: Partition the genome in windows and generate graph files for each chromsome
        s_logger.info(' ')
        logger.run("Partition the genome in windows and generate summary files.")
        total_tag_in_windows = run_make_graph_file_by_chrom.main(args, pool)

        # Step4+5: Normalize and generate WIG file
        s_logger.info(' ')
        logger.run('Normalize graphs by total island filtered reads per million and generating summary WIG file.')
        output_WIG_name = treatment_file_name.replace('.bed', '') + "-W" + str(args.window_size) + "-normalized.wig"
        make_normalized_wig.main(args, output_WIG_name, pool)

        # Step 6: Find candidate islands exhibiting clustering
        s_logger.info(' ')
        logger.run("Find candidate islands exhibiting clustering.")
        total_number_islands = find_islands_in_pr.main(args, total_tag_in_windows, pool)

        # Running SICER with a control library
        if control_lib_exists:
            if args.fdr_method.lower() == "bh":
                fdr_method = "Benjamini-Hochberg"
            elif args.fdr_method.lower() == "clipper":
                fdr_method = "Clipper"
            else:
                args.fdr_method = "bh"
                fdr_method = "Benjamini-Hochberg"
                s_logger.info(' ')
                s_logger.info(f"Invalid FDR approach. Defaulting to Benjamini-Hochberg procedure.")

            s_logger.info(' ')
            logger.run(f"Running SICER with a control library "
                       f"using {fdr_method} procedure as the FDR control approach.")

            # Step 7
            if args.fdr_method.lower() == "bh":
                s_logger.info(' ')
                logger.run("Calculate significance of candidate islands using the control library.")
                associate_tags_with_chip_and_control_w_fc_q.main(args=args,
                                                                 chip_library_size=total_treatment_read_count,
                                                                 control_library_size=total_control_read_count,
                                                                 pool=pool)
                # Step 8: Filter out any significant islands whose pvalue is greater than the false discovery rate
                s_logger.info(' ')
                logger.run(f"Identify significant islands using FDR criterion.")
                # 7 represents the ith column we want to filtered by
                significant_read_count, total_island_count = filter_islands_by_significance.main(args, 7, pool)
                s_logger.info(f"Out of the {total_treatment_read_count} reads in {treatment_file_name}, "
                              f"{significant_read_count} reads are in significant islands")
                s_logger.info(f"Given significance {str(args.false_discovery_rate)}, "
                              f"out of the {total_number_islands} candidate islands, "
                              f"there are {total_island_count} significant islands,")
            elif args.fdr_method.lower() == "clipper":
                s_logger.info(' ')
                logger.run("Calculate significance of candidate islands using the control library.")
                associate_tags_with_chip_and_control_w_fc_q.main(args=args,
                                                                 chip_library_size=total_treatment_read_count,
                                                                 control_library_size=total_control_read_count,
                                                                 pool=pool)
                s_logger.info(' ')
                logger.run("Identify significant islands using Clipper FDR control approach.")
                clipper_significant_read_count, total_read_count = clipper_fdr_control.main(
                    args=args,
                    total_treatment_read_count=total_treatment_read_count,
                    total_control_read_count=total_control_read_count,
                    pool=pool)
                s_logger.info(f"Out of the {total_treatment_read_count} reads in {treatment_file_name}, "
                              f"{total_read_count} reads are in significant islands")
                s_logger.info(f"Given significance {str(args.false_discovery_rate)}, "
                              f"out of the {total_number_islands} candidate islands, "
                              f"there are {clipper_significant_read_count} significant islands.")
            else:
                s_logger.error("Invalid FDR approach. Please choose either 'bh' or 'clipper'.")

        # Optional Outputs
        if (args.significant_reads):
            # Step 9: Filter treatment reads by the significant islands found from step 8
            s_logger.info(' ')
            logger.run("Filter reads with identified significant islands.")
            filter_raw_tags_by_islands.main(args, pool, control_lib_exists)

            # Step 10: Produce graph file based on the filtered reads from step 9
            s_logger.info(' ')
            logger.run("Make summary graph with filtered reads.")
            run_make_graph_file_by_chrom.main(args, pool, True)
            # Step 11: Produce Normalized WIG file
            s_logger.info(' ')
            logger.run("Normalize graphs by total island filitered reads per million and generating summary WIG file.")
            output_WIG_name = (treatment_file_name.replace('.bed', '') + "-W" + str(args.window_size) +
                               "-G" + str(args.gap_size) + "-FDR" + str(args.false_discovery_rate) +
                               "-islandfiltered-normalized.wig")
            make_normalized_wig.main(args, output_WIG_name, pool)

        pool.close()
        pool.join()
        # Final Step
        if df_run == True:
            return temp_dir, total_treatment_read_count
        else:
            s_logger.info(' ')
            logger.run("End of SICER 2")
    finally:
        if df_run == False:
            s_logger.info("Remove temporary directory and all files in it.")
            shutil.rmtree(temp_dir)
