# Python Imports
import os, sys, errno
import pybedtools

# SICER Internal Imports
from sicer.lib import GenomeData
from shared.read import Read

# Cython Imports
from libc.stdlib cimport free, atoi
from libc.string cimport strtok
from libc.stdio cimport fopen, fclose, getline
from libcpp.vector cimport vector
from libcpp.unordered_map cimport unordered_map

cdef createRead(char* line):   
    cdef char* chrom = strtok(&line[0], "\t")
    cdef int start = atoi(strtok(NULL, "\t"))
    cdef int end = atoi(strtok(NULL, "\t"))
    cdef char* name = strtok(NULL, "\t")
    cdef int score = atoi(strtok(NULL, "\t"))
    cdef char strand = strtok(NULL, "\t")[0]

    return Read(chrom, start, end, name, score, strand)

cdef read_bed_file(file, species):
    fp = fopen(file.encode(), "r")
    if fp == NULL:
        raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), file)
        sys.exit(2)

    cdef char *line = NULL
    cdef size_t len = 0

    chrom_list = list(map(lambda x: x.encode(), GenomeData.species_chroms[species]))
    read_lists = {chrm: [] for chrm in chrom_list}
    while getline(&line, &len, fp) != -1:
        read = createRead(line)
        if read.chrom in chrom_list:
            read_lists[read.chrom].append(read)
    
    fclose(fp)
    free(line)

    return read_lists

def read_bed(file, species):
    return read_bed_file(file, species)


