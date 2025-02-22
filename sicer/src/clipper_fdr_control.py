import pandas as pd
import numpy as np
import logging

from functools import partial
from sicer.lib import GenomeData
import os


def expand_intervals(data_array):
    # expanded data use ndarray to store
    expanded_data = []

    for row in data_array:
        chrom, start, end, val1, val2, val3 = row
        start, end = int(start), int(end)
        length = end - start

        for i in range(length):
            new_start = start + i
            new_end = new_start + 1
            expanded_data.append([chrom, int(new_start), int(new_end), float(val1), float(val2), float(val3)])

    return np.array(expanded_data, dtype=object)


def index_bp_conversion(bdg_file, window_size=200):
    """
    Scale down the index by window size to reduce the computation burden for Clipper
    """
    return bdg_file.assign(
        start=(bdg_file['start'] / window_size).astype(int),
        end=((bdg_file['end'] + 1) / window_size).astype(int)
    )


def process_chrom(args, t, compare, chrom):
    union_read_file_name = (args.treatment_file.replace('.bed', '') + '_' +
                            args.control_file.replace('.bed', '') + f'_reads_{chrom}_union.npy')
    try:
        union_read_file = np.load(union_read_file_name, allow_pickle=True)
        union_read_file[:, 1:6] = union_read_file[:, 1:6].astype(np.float32)

        if compare == 'lr':
            return np.sum(union_read_file[:, 5] >= t)
        elif compare == 'sl':
            return np.sum(union_read_file[:, 5] <= -t)
        else:
            return 0
    except:
        return 0


def sum_contrast_score(args, chroms, t, compare, pool):
    partial_process_chrom = partial(process_chrom, args, t, compare)
    sum_values = pool.map(partial_process_chrom, chroms)

    return np.sum(sum_values)


def get_contrast_score_partial(args, chrom):
    try:
        union_read_file_name = args.treatment_file.replace('.bed', '') + '_' + \
                                 args.control_file.replace('.bed', '') + f'_reads_{chrom}_union.npy'
        union_read_file = np.load(union_read_file_name, allow_pickle=True)
        contrast_score = np.nan_to_num(union_read_file[:, 5], nan=0)
        return np.sort(np.unique(np.abs(contrast_score[contrast_score != 0])))
    except:
        return np.array([])


def clipper_BC(FDR=0.05, args=None, chroms=None, pool=None):
    '''
    BC procedure calculation taking window-wise contrast score as input to obtain the cutoff
    '''
    contrast_scores = pool.map(partial(get_contrast_score_partial, args), chroms)
    c_abs = np.sort(np.unique(np.concatenate(contrast_scores)))

    emp_fdp = np.full(len(c_abs), np.nan)
    emp_fdp[0] = 1
    for i in range(1, len(c_abs)):
        t = c_abs[i]
        emp_fdp[i] = ((1 + sum_contrast_score(args, chroms, t, 'sl', pool)) /
                      max(sum_contrast_score(args, chroms, t, 'lr', pool), 1))
        emp_fdp[i] = min(emp_fdp[i], emp_fdp[i - 1])

        if emp_fdp[i] <= FDR:
            break

    c_abs = c_abs[~np.isnan(emp_fdp)]
    emp_fdp = emp_fdp[~np.isnan(emp_fdp)]

    indices = np.where(emp_fdp <= FDR)
    if indices[0].size == 0:
        thre = -np.inf
    else:
        thre = c_abs[np.min(indices)]

    return thre


def compute_importance_score(args, factors, chrom):
    (treat_scaling_factor, ctrl_scaling_factor) = factors
    chrom_island_load_name = (args.treatment_file.replace('.bed', '') + '_'
                              + args.control_file.replace('.bed', '') + f'_reads_{chrom}_union.npy')
    try:
        chrom_island_load = np.load(chrom_island_load_name, allow_pickle=True)
    except:
        return
    chrom_island_load = np.c_[chrom_island_load, np.zeros(len(chrom_island_load))]
    chrom_island_load[:, 3] = chrom_island_load[:, 3] / treat_scaling_factor
    chrom_island_load[:, 4] = chrom_island_load[:, 4] / ctrl_scaling_factor
    chrom_island_load[:, 5] = chrom_island_load[:, 3] - chrom_island_load[:, 4]
    # unfold the island to one by one window sized region
    unfold_chrom_island_load = expand_intervals(chrom_island_load)
    # unfold_chrom_island_load = chrom_island_load
    np.save(chrom_island_load_name, unfold_chrom_island_load)


def clipper_threshold(FDR=0.05, args=None, pool=None, chroms=None, factors=None):
    '''
    Compute the cutoff of the island contrast score for Clipper
    '''
    compute_importance_score_partial = partial(compute_importance_score, args, factors)
    pool.map(compute_importance_score_partial, chroms)
    threshold = clipper_BC(FDR=FDR, args=args, chroms=chroms, pool=pool)
    return threshold


def insert_gaps(table):

    table = table.to_numpy()
    new_table = []

    # Insert starting row from 0 to the start of the first row, if it doesn't start at 0
    if table[0][1] != 0:
        new_table.append((table[0][0], 0, table[0][1]-1, 0, 0))

    for i in range(len(table)):
        new_table.append(table[i])
        if i < len(table) - 1 and table[i][0] == table[i + 1][0]:  # same chromosome
            end_current = table[i][2]
            start_next = table[i + 1][1]
            if end_current + 1 < start_next:  # check for gap
                new_row = (table[i][0], end_current+1, start_next-1, 0, 0)
                new_table.append(new_row)
    new_table = pd.DataFrame(new_table, columns=['chrom', 'start', 'end', 'treat', 'ctrl']).astype(
        {'chrom': str, 'start': int, 'end': int, 'treat': int, 'ctrl': int})
    return new_table


def filter_and_expand_rows(base_df, filter_df):
    '''
    Filter the rows in genome-wide sicer windows that are within the range of the rows in sicer candidate island
    '''
    starts = base_df['start'].values
    ends = base_df['end'].values

    expanded_rows = []
    for _, filter_row in filter_df.iterrows():
        mask = (starts >= filter_row['start']) & (ends <= filter_row['end'])
        matched_rows = base_df[mask]
        expanded_rows.extend(matched_rows.values)

    expanded_df = pd.DataFrame(expanded_rows, columns=base_df.columns)
    return expanded_df


def islands_bdg_union_partial(args, chrom):
    treat_file = args.treatment_file.replace('.bed', '') + f'_{chrom}' + '_graph.npy'
    ctrl_file = args.control_file.replace('.bed', '') + f'_{chrom}' + '_graph.npy'
    island_file = args.treatment_file.replace('.bed', '') + f'_{chrom}_island_summary.npy'

    treat_file = 'reads_' + treat_file
    ctrl_file = 'reads_' + ctrl_file

    try:
        treat_bdg = np.load(treat_file, allow_pickle=True)
        ctrl_bdg = np.load(ctrl_file, allow_pickle=True)
        is_island_file = np.load(island_file, allow_pickle=True)[:, :5]
        island_bdg = pd.DataFrame(is_island_file, columns=['chrom', 'start', 'end', 'treat', 'ctrl']).astype(
            {'chrom': str, 'start': int, 'end': int, 'treat': int, 'ctrl': int})
    except:
        return

    if treat_bdg.shape[0] == 0 and ctrl_bdg.shape[0] == 0:
        return
    elif treat_bdg.shape[0] == 0:
        treat_bdg = pd.DataFrame(columns=['chrom', 'start', 'end', 'treat'])
        ctrl_bdg = pd.DataFrame(ctrl_bdg[ctrl_bdg[:, 3] != 0], columns=['chrom', 'start', 'end', 'ctrl'])
    elif ctrl_bdg.shape[0] == 0:
        ctrl_bdg = pd.DataFrame(columns=['chrom', 'start', 'end', 'ctrl'])
        treat_bdg = pd.DataFrame(treat_bdg[treat_bdg[:, 3] != 0], columns=['chrom', 'start', 'end', 'treat'])
    else:
        treat_bdg = pd.DataFrame(treat_bdg[treat_bdg[:, 3] != 0], columns=['chrom', 'start', 'end', 'treat'])
        ctrl_bdg = pd.DataFrame(ctrl_bdg[ctrl_bdg[:, 3] != 0], columns=['chrom', 'start', 'end', 'ctrl'])

    ctrl_bdg = ctrl_bdg.astype({"chrom": str, "start": int, "end": int, "ctrl": int})
    treat_bdg = treat_bdg.astype({"chrom": str, "start": int, "end": int, "treat": int})
    merged_df_outer_join = pd.merge(ctrl_bdg, treat_bdg, how='outer', on=['chrom', 'start', 'end']).fillna(0)
    merged_df_outer_join = merged_df_outer_join.sort_values(by=['start'])
    # reindex
    merged_df_outer_join.index = range(merged_df_outer_join.shape[0])
    merged_df_outer_join = merged_df_outer_join[['chrom', 'start', 'end', 'treat', 'ctrl']]
    # insert gaps and convert index to bp
    merged_df_outer_join = insert_gaps(merged_df_outer_join)
    merged_df_outer_join = index_bp_conversion(merged_df_outer_join, window_size=args.window_size)
    island_bdg = index_bp_conversion(island_bdg, window_size=args.window_size)

    name_for_save_island = (args.treatment_file.replace('.bed', '') + '_' + args.control_file.replace('.bed', '')
                            + f'_read_{chrom}_island_summary_clipper_converted.npy')

    np.save(name_for_save_island, island_bdg.to_numpy())

    save_union_read_candidate_island_wide = (args.treatment_file.replace('.bed', '') + '_'
                                             + args.control_file.replace('.bed', '') + f'_reads_{chrom}_union.npy')

    union_read_candidate_island_wide = filter_and_expand_rows(merged_df_outer_join, island_bdg)
    # save to npy file
    np.save(save_union_read_candidate_island_wide, union_read_candidate_island_wide.to_numpy())


def process_chromosome_partial(args, threshold, chrom):
    chrom_island_load = (args.treatment_file.replace('.bed', '') + '_' + args.control_file.replace('.bed', '')
                         + f'_read_{chrom}_island_summary_clipper_converted.npy')
    chrom_read_union_load = (args.treatment_file.replace('.bed', '') + '_' +
                             args.control_file.replace('.bed', '') + f'_reads_{chrom}_union.npy')

    try:
        chrom_island_list = pd.DataFrame(np.load(chrom_island_load, allow_pickle=True),
                                         columns=['chrom', 'start', 'end', 'treat', 'ctrl']).astype(
            {'chrom': str, 'start': int, 'end': int, 'treat': float, 'ctrl': float})
        read_union = pd.DataFrame(np.load(chrom_read_union_load, allow_pickle=True),
                                  columns=['chrom', 'start', 'end', 'treat', 'ctrl', 'diff']).astype(
            {'chrom': str, 'start': int, 'end': int, 'treat': float, 'ctrl': float, 'diff': float})
    except FileNotFoundError:
        return

    clipper_peak_save = (args.treatment_file.replace('.bed', '') + '_' + args.control_file.replace('.bed', '') +
                         f'_read_{chrom}_island_clipper_filtered.npy')

    filtered_island_list = []
    for i in range(len(chrom_island_list)):
        filtered_reads = filter_rows_within_range(read_union, chrom_island_list.iloc[i]['start'],
                                                  chrom_island_list.iloc[i]['end'])
        # first value of the start and end is the chromosome name
        treat_normalized = filtered_reads['treat']
        ctrl_normalized = filtered_reads['ctrl']
        mean_value = np.mean(treat_normalized - ctrl_normalized)
        if mean_value >= threshold:
            filtered_island_list.append(chrom_island_list.iloc[i])

    filtered_island_list = pd.DataFrame(filtered_island_list)
    np.save(clipper_peak_save, filtered_island_list.to_numpy())


def islands_bdg_union(args, chroms, pool):
    island_bdg_union_pl = partial(islands_bdg_union_partial, args)
    pool.map(island_bdg_union_pl, chroms)


def index_bp_conversion_back(bdg_file, window_size=200):
    return bdg_file.assign(
        start=(bdg_file['start'] * window_size).astype(int),
        end=((bdg_file['end'] * window_size) - 1).astype(int)
    )


def filter_rows_within_range(data_df, start, end):
    filtered_df = data_df[(data_df['start'] < end) & (data_df['end'] > start)]
    return filtered_df


def main(args, total_treatment_read_count, total_control_read_count, pool):
    s_logger = logging.getLogger("s_logger")

    chroms = GenomeData.species_chroms[args.species]
    treat_scaling_factor = total_treatment_read_count / 1000000 * (args.window_size / 1000)
    ctrl_scaling_factor = total_control_read_count / 1000000 * (args.window_size / 1000)

    # preprocess and union reads for clipper procedure
    islands_bdg_union(args, chroms, pool)

    # calculate clipper threshold
    threshold = clipper_threshold(args=args, chroms=chroms, FDR=args.false_discovery_rate, pool=pool,
                                  factors=(treat_scaling_factor, ctrl_scaling_factor))

    if threshold == -np.inf:
        s_logger.info('Cannot find a valid threshold for Clipper, '
                      'please try the Benjamini-Hochberg procedure or stick to the result from Clipper '
                      '(no false discovery).')
    else:
        s_logger.info(f'Clipper contrast score cutoff is {threshold}.')

    process_chrom_func = partial(process_chromosome_partial, args, threshold)
    pool.map(process_chrom_func, chroms)

    # save output filer for clipper filtered island
    island_summary_length = 0
    total_read_count = 0
    try:
        outfile_name = (args.treatment_file.replace('.bed', '') + '-W' + str(args.window_size) + '-G'
                        + str(args.gap_size) + '-FDR' + str(args.false_discovery_rate) + '-island.bed')
    except AttributeError:
        outfile_name = (args.treatment_file.replace('.bed', '') + '-W' + str(args.window_size) +
                        '-FDR' + str(args.false_discovery_rate) + '-island.bed')
    outfile_path = os.path.join(args.output_directory, outfile_name)

    with open(outfile_path, 'w') as outfile:
        for chrom in chroms:
            island_file_name = (args.treatment_file.replace('.bed', '') + '_' + args.control_file.replace('.bed', '') +
                                f'_read_{chrom}_island_clipper_filtered.npy')
            try:
                filtered_island = pd.DataFrame((np.load(island_file_name, allow_pickle=True)),
                                               columns=['chrom', 'start', 'end', 'treat', 'ctrl']).astype(
                    {'chrom': str, 'start': int, 'end': int, 'treat': int, 'ctrl': int})
                # drop the ctrl column
                filtered_island = filtered_island.drop(columns=['ctrl'])
            except:
                continue

            island_summary_length += len(filtered_island)
            filtered_island = index_bp_conversion_back(filtered_island, window_size=args.window_size).to_numpy()

            for island in filtered_island:
                output_line = ''
                for i in range(0, len(island)):
                    output_line += str(island[i]) + '\t'
                output_line += '\n'
                outfile.write(output_line)
                total_read_count += island[3]

    return island_summary_length, total_read_count

