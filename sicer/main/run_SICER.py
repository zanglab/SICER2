# Developed by Zang Lab at University of Virginia - 2018

#Author: Mengxue Tian (edited 2025)


import os
import shutil
import sys
import tempfile
import multiprocessing as mp
import subprocess
import warnings

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
from sicer.lib import Utility

''' args: ArgumentParser object formed form command line parameters
    df_run: If df_run is true, then this instance of SICER is called by SICER-DF module.
            Default value is False.'''

def convert_bam_to_bed(bam_path, args, output_dir=None):
    """
    Convert BAM to BED. Automatically detects whether the BAM is paired-end (PE)
    or single-end (SE), sets args.input_type accordingly, and logs information
    if args.verbose is True.

    Returns the output BED file path.
    """
    print("[INFO] Start generating BED file from BAM.\n")
    if output_dir is None:
        output_dir = tempfile.gettempdir()

    bed_file = os.path.join(output_dir, os.path.basename(bam_path).replace(".bam", ".bed"))

    # Detect paired-end or single-end BAM using samtools flag 0x1
    pe_check = subprocess.Popen(
        f"samtools view -c -f 1 {bam_path}",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True
    )
    out, err = pe_check.communicate()

    if pe_check.returncode != 0:
        sys.stderr.write(
            "Error: Failed to determine BAM type using samtools.\n"
            "Make sure samtools is installed correctly and in your PATH.\n"
        )
        sys.exit(1)

    try:
        num_pe = int(out.strip() or 0)
    except ValueError:
        num_pe = 0

    is_paired_end = num_pe > 0
    detected_type = "PE" if is_paired_end else "SE"

    # 记录一下用户原始的 input_type（可能是命令行给的，也可能是默认 SE）
    user_input_type = getattr(args, "input_type", None)

    # 如果用户显式给了 input_type，且和 BAM flags 检测结果不一致，就给个 warning
    if user_input_type is not None and user_input_type != detected_type:
        warnings.warn(
            f"[SICER WARNING] User specified input_type='{user_input_type}' or NA, "
            f"but BAM flags indicate '{detected_type}'. "
            f"SICER will use '{detected_type}' to ensure correct interpretation."
        )

    # 最终还是以 BAM 检测结果为准，保证不会把 PE 当 SE 或相反
    args.input_type = detected_type


    # Verbose logging
    if args.verbose:
        if is_paired_end:
            print(f"[INFO] Detected paired-end BAM: {bam_path}  (Total read count: {num_pe})\n")
        else:
            print(f"[INFO] Detected single-end BAM: {bam_path}\n")
        print(
            f"[INFO] BAM input detected. Setting input_type='{args.input_type}' "
            "based on BAM flags (0x1).\n"
        )

    # Choose conversion method
    if is_paired_end:
        # BEDPE -> fragment BED6
        cmd = (
            "bedtools bamtobed -bedpe -i {bam} | "
            "awk 'BEGIN{{OFS=\"\\t\"}} "
            "{{ if ($1 == $4) {{ "
            "s = ($2 < $5 ? $2 : $5); "
            "e = ($3 > $6 ? $3 : $6); "
            "print $1, s, e, $7, $8, \"+\" }} "
            "}}' > {bed}"
        ).format(bam=bam_path, bed=bed_file)
    else:
        # SE: standard read-level BED6
        cmd = f"bedtools bamtobed -i {bam_path} > {bed_file}"

    # Execute the conversion
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, shell=True)
    process.communicate()

    if process.returncode != 0:
        sys.stderr.write(
            "Error: Failed to convert BAM to BED.\n"
            "Check if bedtools2 and samtools are installed correctly.\n"
        )
        sys.exit(1)

    if args.verbose:
        print(f"[INFO] Created BED file: {bed_file}\n\n")

    return bed_file







def main(args, df_run=False):
    """
    Main SICER pipeline.

    args:
        ArgumentParser object formed from command line parameters.
        Expected to contain:
            - treatment_file / control_file
            - species, window_size, gap_size, fragment_size, redundancy_threshold
            - input_type: 'SE' or 'PE'
    df_run:
        If df_run is True, then this instance of SICER is called by SICER-DF module.
        Default value is False.
    """

    # Check SE/PE mode and print a summary
    input_type = getattr(args, "input_type", "SE")
    if input_type not in ["SE", "PE"]:
        sys.stderr.write(
            f"Error: Unsupported input_type '{input_type}'. Only 'SE' and 'PE' are allowed.\n"
        )
        sys.exit(1)

    if input_type == "SE":
        print(
            "Running SICER in SE mode:\n"
            "  - Input is expected to be single-end read-level BED/BAM with valid +/- strand.\n"
            "  - Tag position will be computed by shifting reads by fragment_size/2 based on strand.\n"
        )
    else:  # PE
        print(
            "Running SICER in PE mode:\n"
            "  - Input is expected to be fragment-level BED/BAM (e.g. from paired-end data).\n"
            "  - Each BED interval is treated as a fragment, and tag position is the fragment midpoint.\n"
            "  - Strand information is ignored in this mode.\n"
        )

    # Checks if there is a control library
    control_lib_exists = True
    if args.control_file is None:
        control_lib_exists = False

    # Creates temporary directory to contain all intermediate files.
    try:
        temp_dir = tempfile.mkdtemp()
        # Change current working directory to temp_dir
        os.chdir(temp_dir)
    except Exception:
        sys.exit(
            "Temporary directory required for SICER cannot be created. "
            f"Check if directories can be created in {curr_path}."
        )

    try:
        # Step 0: create Pool object for parallel processing
        num_chroms = len(GenomeData.species_chroms[args.species])
        pool = mp.Pool(processes=min(args.cpu, num_chroms))

        # Step 1: Remove redundant reads/fragments in treatment according to input threshold
        # Output is the total number of retained tags. Represents size of the library.
        treatment_file_name = os.path.basename(args.treatment_file)
        print(
            "Preprocess the", treatment_file_name,
            "file to remove redundancy with threshold of",
            args.redundancy_threshold, "\n"
        )
        
        if args.treatment_file.endswith(".bam"):
            bed_path = convert_bam_to_bed(args.treatment_file, args, temp_dir)  # auto recognize SE/PE
            args.treatment_file = bed_path
            treatment_file_name = os.path.basename(bed_path)
        
        total_treatment_read_count = remove_redundant_reads.main(args, args.treatment_file, pool)
        args.treatment_file = treatment_file_name
        print('\n')

        # Step 2: Remove redundant reads/fragments in control library according to input threshold
        if control_lib_exists:
            control_file_name = os.path.basename(args.control_file)
            print(
                "Preprocess the", control_file_name,
                "file to remove redundancy with threshold of",
                args.redundancy_threshold, "\n"
            )
        
            if args.control_file.endswith(".bam"):
                bed_path = convert_bam_to_bed(args.control_file, args, temp_dir)
                args.control_file = bed_path
                control_file_name = os.path.basename(bed_path)
        
            total_control_read_count = remove_redundant_reads.main(args, args.control_file, pool)
            args.control_file = control_file_name
            print('\n')

        # Step 3: Partition the genome in windows and generate graph files for each chromosome
        print("Partition the genome in windows and generate summary files... \n")
        total_tag_in_windows = run_make_graph_file_by_chrom.main(args, pool)
        print("\n")

        # Step 4+5: Normalize and generate WIG file
        print("Normalizing graphs by total island filtered reads/fragments per million "
              "and generating summary WIG file...\n")
        output_WIG_name = (
            treatment_file_name.replace('.bed', '') +
            "-W" + str(args.window_size) +
            "-normalized.wig"
        )
        make_normalized_wig.main(args, output_WIG_name, pool)

        # Step 6: Find candidate islands exhibiting clustering
        print("Finding candidate islands exhibiting clustering...\n")
        find_islands_in_pr.main(args, total_tag_in_windows, pool)
        print("\n")

        # Running SICER with a control library
        if control_lib_exists:
            # Step 7
            print("Calculating significance of candidate islands using the control library... \n")
            associate_tags_with_chip_and_control_w_fc_q.main(
                args,
                total_treatment_read_count,
                total_control_read_count,
                pool
            )

            # Step 8: Filter out significant islands by FDR
            print("Identify significant islands using FDR criterion\n")
            # 7 represents the ith column used for filtering in the island file
            significant_read_count = filter_islands_by_significance.main(args, 7, pool)
            print(
                "Out of the", total_treatment_read_count,
                "tags in", treatment_file_name, ",",
                significant_read_count,
                "tags are in significant islands\n"
            )

        # Optional Outputs
        if args.significant_reads:
            # Step 9: Filter treatment reads by the significant islands from step 8
            print("Filtering reads with identified significant islands...\n")
            filter_raw_tags_by_islands.main(args, pool)

            # Step 10: Produce graph file based on the filtered reads
            print("Making summary graph with filtered reads...\n")
            run_make_graph_file_by_chrom.main(args, pool, True)

            # Step 11: Produce Normalized WIG file based on filtered reads
            print("\nNormalizing graphs by total island filtered reads per million "
                  "and generating summary WIG file...\n")
            output_WIG_name = (
                treatment_file_name.replace('.bed', '') +
                "-W" + str(args.window_size) +
                "-G" + str(args.gap_size) +
                "-FDR" + str(args.false_discovery_rate) +
                "-islandfiltered-normalized.wig"
            )
            make_normalized_wig.main(args, output_WIG_name, pool)

        pool.close()
        pool.join()

        # Final Step
        if df_run is True:
            return temp_dir, total_treatment_read_count
        else:
            print("End of SICER")

    finally:
        if df_run is False:
            print("Removing temporary directory and all files in it.")
            shutil.rmtree(temp_dir)
    
    
