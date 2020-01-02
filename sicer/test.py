# Python Imports
import multiprocessing as mp
import cProfile
import time
import timeit

# SICER Internal Imports
from sicer.bed_reader import BEDReader

PROFILE = False

if PROFILE:
    pr = cProfile.Profile()
    pr.enable()

bed_reader = BEDReader("../../data/36629_treat_rep1.bed", "hg38", 12, 1)
reads = bed_reader.read_file()

if PROFILE:
    pr.disable()
    pr.print_stats(sort='time')
