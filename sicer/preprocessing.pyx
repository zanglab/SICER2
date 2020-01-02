# Python Imports
import os, errno

# SICER Internal Imports
from sicer.shared import GenomeData
from sicer.shared.data_classes cimport BEDRead
from sicer.shared.chrom_collections cimport ChromBEDReadCollection

# Cython Imports
cimport cython
from libcpp cimport bool
from cython.parallel import parallel, prange
from libcpp.algorithm cimport sort
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdlib cimport malloc, free
from cython.operator cimport dereference as deref

ctypedef char* cstr
ctypedef bool (*cmp_f)(BEDRead, BEDRead)

cdef bool compare_reads(BEDRead i, BEDRead j) nogil:
    if i.strand != j.strand:
        return i.strand < j.strand
    else:
        if i.start != j.start:
            return i.start < j.start
        else:
            return i.end < j.end

cdef void _preprocess_BED_reads_by_chrom(vector[BEDRead] vec, int threshold) nogil:

    sort[vector[BEDRead].iterator, cmp_f](vec.begin(), vec.end(), compare_reads)

cpdef ChromBEDReadCollection preprocess_BED_reads(
    ChromBEDReadCollection reads,
    int num_cpu,
    int redundancy_threshold
):
    
    # Convert Python list to C-array for no-GIL use
    cdef int num_chroms = len(reads.getChromosomes())
    cdef cstr *chroms = <cstr *>malloc(num_chroms*cython.sizeof(cstr))
    if chroms is NULL:
        raise MemoryError()
    for j in range(num_chroms):
        chroms[j] = reads.getChromosomes()[j]

    cdef mapcpp[string, vector[BEDRead]] data = reads.getData()
    cdef int i
    for i in prange(num_chroms, schedule='guided', num_threads=num_cpu, nogil=True):
        _preprocess_BED_reads_by_chrom(data[chroms[i]], redundancy_threshold)

    free(chroms)
    return reads