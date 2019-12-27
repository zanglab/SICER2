from read_bed import read_bed
from sicer.lib import GenomeData

import multiprocessing as mp
import cProfile
import time
import pybedtools

PROFILE = False

if PROFILE:
    pr = cProfile.Profile()
    pr.enable()
'''
bed_reads = read_bed("../../../data/36629_treat_rep1.bed", "hg38")
num_chroms = len(GenomeData.species_chroms["hg38"])
worker_pool = mp.Pool(processes=min(6, num_chroms))

'''
a = read_bed("../../../data/36629_treat_rep1.bed", "hg38")
b = read_bed("../../../data/test.bam", "hg38")

if PROFILE:
    pr.disable()
    pr.print_stats(sort='time')