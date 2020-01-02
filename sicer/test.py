# Python Imports
import multiprocessing as mp
import cProfile
import time

# SICER Internal Imports
from sicer.bed_reader import BEDReader
from sicer.preprocessing import preprocess_BED_reads

PROFILE = False

print("Start")

if PROFILE:
    pr = cProfile.Profile()
    pr.enable()

bed_reader = BEDReader("../../data/36629_treat_rep1.bed", "hg38")
reads = bed_reader.read_file()
reads.printDataHead()
print("============================================")
pp_reads = preprocess_BED_reads(reads, 6, 1)
pp_reads.printDataHead()

if PROFILE:
    pr.disable()
    pr.print_stats(sort='time')

print("End")