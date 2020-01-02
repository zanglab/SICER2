#Author: Chongzhi Zang, Yiren Wang
#Edits made by Jeffrey Yoo

import multiprocessing as mp
import os
from functools import partial
from math import *

import numpy as np

# From SICER package
from sicer.lib import GenomeData


'''version 8: 3-phase coarse graining, take the phase that has most 1 to next step. '''

def linreg(list X, list Y):
	"from Simple Recipes in Python http://www.phys.uu.nl/~haque/computing/WPark_recipes_in_python.html"
	if len(X) != len(Y):
		raise (ValueError, 'unequal length')
	N = len(X)
	cdef double Sx, Sy, Sxx, Syy, Sxy
	Sx = Sy = Sxx = Syy = Sxy = 0.0
	for i in range(0,len(X)):
		x=X[i]
		y=Y[i]
		Sx = Sx + x
		Sy = Sy + y
		Sxx = Sxx + x*x
		Syy = Syy + y*y
		Sxy = Sxy + x*y
	det = Sxx * N - Sx * Sx
	if det != 0:
		return (Sxy * N - Sy * Sx)/det
	else:
		return 0

def is_list_sorted(List):
		"""
		Check if sorted in ascending order.
		input is a list of pure numbers.
		output: sorted =1 or 0
		"""
		sorted = 1;
		for index in range(0, len(List)-1):
			if List[index] > List[index + 1]:
				sorted = 0;
		return sorted;


cdef float start_list_correlation_r_rev(list List, int win, int r, int chrom_length):
	'''List must be sorted'''
	assert is_list_sorted(List) == 1
	cdef int x, d, n, i, SUMM
	cdef list a
	x = List[0]%win
	d = r // win
	SUMM = 0
	n = (chrom_length - x) // win
	if n - d > 0:
		a = [0] * n
		for island in List:
			i = (island - x) // win
			if i >= 0 and i < n:
				a[i] = 1
		for i in range(0, n - d):
			SUMM += a[i] * a[i + d]
		return SUMM / float(n - d) - ((sum(a) / float(len(a))) ** 2)
	else:
		return 0.0


def start_list_correlation_function(List, win, chrom_length, name):
	xlist = []
	ylist = []
	#file = open("cr_"+name+"_"+str(win)+".txt", 'w')
	for i in range(0, min(3, int(chrom_length/win))):
		r = i * win
		c = start_list_correlation_r_rev(List, win, r, chrom_length)
		xlist.append(i)
		ylist.append(c)
		#file.write(str(i)+'\t'+str(c)+'\n')
	#file.close()
	return (xlist, ylist)


def correlation_length_fit(xlist, ylist):
	assert len(xlist) == len(ylist)
	loglist = []
	for i in range(0, len(ylist)):
		loglist.append(log(max(ylist[i], 0.000000000001)))
	a = linreg(xlist[1:],loglist[1:])
	if abs(a) > 0.000000000001:
		return -1.0/a
	else:
		return 1e12


cdef list graining(list List, int win, int step, int score):
	'''
	1 step coarse graining, phase considered:
	List must be sorted!
	List (list) contains (start) coordinates of positive signals;
	win (int) is the window size in list, coarse graining will start from this resolution;
	step (int) is the number of windows in one graining unit;
	score (int) is the minimum number of positive elements in the graining unit to call the unit positive;
	output is a list of positive unit number in each graining step;
	'''
	result = []
	endlimit = List[-1]
	cdef int i, j, h, k, n
	for p in range(0, step):
		tmp_result = []
		i = List[0] - p * win
		k = 0
		while i <= endlimit and k < len(List):
			j = i + step * win
			h = k
			while h <= (len(List) - 1) and List[h] < j:
				h += 1
			n = h - k
			if n >= score:
				tmp_result.append(i)
			k = h
			i = j
		if len(tmp_result) > len(result):
			result = tmp_result
	return (result)


def coarsegraining(List, win_min, step, score, genome_length):
	if (is_list_sorted(List) != 1):
		List.sort()

	Length_list = []
	Length_list.append(len(List))  #number of eligible windows
	result_list = []
	result_list.append(List)	#list of start positions of eligible windows
	win = win_min
	while len(List) > 0:

		List = graining(List, win, step, score)
		Length_list.append(len(List))
		if len(List) > 0:
			result_list.append(List)
		win = win * step
	return result_list


def union_islands_to_list(islandlist, win):
	'''input islandlist and output list are both lists of BED island objects'''
	islandlist.sort(key=lambda x: x[1]);
	List = []
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current[1] <= compare[1]
		if compare[1] > current[2] + 1 + win:
			List.append(current)
			current = compare
			i += 1
		else:
			current = [current[0], current[1], max(current[2], compare[2])]
			i += 1
	List.append(current)
	return List


def write_islandlist(List, win, chrom):
	'''input a start list and universal island width, output a islandlist of BED objects
	object.start = List[i]
	object.end = List[i] + win - 1'''
	output_list = []
	for start in List:
		island = [chrom, start, int(start + win - 1)]
		output_list.append(island)
	return output_list


def backstep(islandlist, List, win, chrom):
	'''one step trace back
		island[0]=chromosome
		island[1]=start position of island
		island[2]=end position of island'''
	# result_list = []
	# fine_islands = []
	addtional_islands = write_islandlist(List, win, chrom)
	for island in islandlist:
		start_left = (island[1] - win) in List
		start_right = island[1] in List
		if start_left and start_right:
			island[1] = island[1] - win
		elif (not start_left) and (not start_right):
			island[1] = island[1] + win

		end_left = (island[2] + 1 - win) in List
		end_right = (island[2] + 1) in List
		if end_left and end_right:
			island[2] = island[2] + win
		elif (not end_left) and (not end_right):
			island[2] = island[2] - win
		assert island[1] < island[2]

	return union_islands_to_list(islandlist + addtional_islands, win)


def traceback(List, win_min, step, level, genome_length, chrom):
	'''
	Input is a list of lists.
	'''

	win = int(win_min * (step ** (len(List) - 1)))

	islandlist = write_islandlist(List[-1], win, chrom)
	backlist = List[-1]
	(xlist, ylist) = start_list_correlation_function(backlist, win, genome_length, chrom)
	correlation_length = correlation_length_fit(xlist, ylist)

	if len(List) > 1:
		(xlist, ylist) = start_list_correlation_function(List[-2], win // step, genome_length, chrom)
		correlation_length_next = correlation_length_fit(xlist, ylist)

	i = 1
	while i < len(List) - level:
		backlist = List[-i - 1]
		win = win // step

		if correlation_length > 1.0 and correlation_length_next >= correlation_length:
			break;
		else:
			islandlist = write_islandlist(backlist, win, chrom)
			correlation_length = correlation_length_next
			if len(List) > i + 1:
				(xlist, ylist) = start_list_correlation_function(List[-i - 2], win // step, genome_length, chrom)
				correlation_length_next = correlation_length_fit(xlist, ylist)
				# print len(List[-i-2]), correlation_length_next
			else:
				correlation_length_next = 10000
		i += 1
	while i < len(List) - level:
		backlist = List[-i - 1]
		# print len(islandlist)
		islands = islandlist
		islandlist = backstep(islands, backlist, win, chrom)
		win = win // step
		i += 1
	return islandlist

def filter_and_find_islands(args, min_tag_count, chrom):
	graph_file = args.treatment_file.replace('.bed', '') + '_' + chrom + '_graph.npy'
	chrom_windows = np.load(graph_file, allow_pickle=True)
	print_return = ''
	total_count_island = 0
	if (len(chrom_windows) > 0):
		chrom_lengths = GenomeData.species_chrom_lengths[args.species][chrom]
		eligible_start_list = []

		for i in range(0, len(chrom_windows)):
			window = chrom_windows[i]
			read_count = window[3]
			if read_count >= min_tag_count:
				eligible_start_list.append(window[1])
		island_list = coarsegraining(eligible_start_list, args.window_size, args.step_size, args.step_score,
									 chrom_lengths)
		islands = traceback(island_list, args.window_size, args.step_size, 0, chrom_lengths, chrom)

		if not (len(islands) > 0):
			print_return += chrom + " does not have any islands meeting the required significance"
		for i in range(len(islands)):
			island = islands[i]
			islands[i] = (island[0], int(island[1]), int(island[2]))

		np_islands = np.array(islands, dtype=object)
		np.save(graph_file, np_islands)
		total_count_island = len(islands)

	return (total_count_island, print_return)

def main(args, read_count, pool):
	print("Coarse-graining approach to identify ChIP-Seq enriched domains:")
	print("Species: %s" % args.species)
	print("Window_size: %d " % args.window_size)
	print("Coarse graining step: %d" % args.step_size)
	print("Coarse graining score: %d " % args.step_score)

	chroms = GenomeData.species_chroms[args.species]
	total_read_count = read_count
	print("Total read count: %d" % total_read_count)

	genome_length = sum(GenomeData.species_chrom_lengths[args.species].values());
	effective_genome_length = int(args.effective_genome_fraction * genome_length);

	average = float(total_read_count) * args.window_size / effective_genome_length;
	print("Effective genome length: %d " % effective_genome_length)
	print("Window average: %f" % average)

	min_tags_in_window = int(average) + 1
	print("Minimum read count in a qualified window: %d\n" % min_tags_in_window)

	print("Running coarsegraining method... (this might take some time)")
	# read in the summary graph file
	# Use multiprocessing to find islands separately by chromosome
	# pool = mp.Pool(processes=min(args.cpu, len(chroms)))
	filter_and_find_islands_partial = partial(filter_and_find_islands, args, min_tags_in_window)
	filtered_islands_result = pool.map(filter_and_find_islands_partial, chroms)
	# pool.close()

	file_name = args.treatment_file.replace('.bed', '')
	outfile_path = os.path.join(args.output_directory, (file_name + '-W' + str(args.window_size) + '.cgisland'))
	total_number_islands = 0	
	path_to_filtered_graph = []
	with open(outfile_path, 'w') as outfile:
		for i in range(0, len(filtered_islands_result)):
			filtered_chrom_graph = np.load(file_name + '_' + chroms[i] + '_graph.npy', allow_pickle=True)
			total_number_islands += filtered_islands_result[i][0]
			if (filtered_islands_result[i][1] != ""):
				print(filtered_islands_result[i][1])
			for island in filtered_chrom_graph:
				line = (island[0] + '\t' + str(int(island[1])) + '\t' + str(int(island[2]))
						+ '\t' + '1\n')
				outfile.write(line)

	print("Total number of islands: %d" % total_number_islands);
