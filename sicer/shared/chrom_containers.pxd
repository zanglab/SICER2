# TODO: Using a templated class would be nice, but Cython doesn't support it currently

# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Window, Island, DiffExprIsland

# Cython Imports
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from libc.stdint cimport uint32_t

ctypedef char* cstr

cdef class BEDReadContainer:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vector[BEDRead]] data
        uint32_t read_count

    cdef void insertRead(self, string chrom, BEDRead item)
    cpdef void updateReadCount(self)
    cpdef uint32_t getReadCount(self)
    cpdef list getChromosomes(self)
    cdef mapcpp[string, vector[BEDRead]] getData(self)
    cdef vector[BEDRead]* getVectorPtr(self, string chrom) nogil


cdef class WindowContainer:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vector[Window]] data
        uint32_t window_count
        uint32_t total_tag_count

    cpdef void updateCounts(self)
    cpdef uint32_t getWindowCount(self)
    cpdef list getChromosomes(self)
    cpdef uint32_t getTotalTagCount(self)
    cdef mapcpp[string, vector[Window]] getData(self)
    cdef vector[Window]* getVectorPtr(self, string chrom) nogil


cdef class IslandContainer:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vector[Island]] data
        uint32_t island_count

    cpdef void updateIslandCount(self)
    cpdef uint32_t getIslandCount(self)
    cpdef list getChromosomes(self)
    cdef mapcpp[string, vector[Island]] getData(self)
    cdef vector[Island]* getVectorPtr(self, string chrom) nogil


cdef class DiffExprIslandContainer:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vector[DiffExprIsland]] data
        uint32_t island_count
        vector[double] pvalue_list

    cpdef void updateIslandCount(self)
    cpdef uint32_t getIslandCount(self)
    cpdef list getChromosomes(self)
    cdef mapcpp[string, vector[DiffExprIsland]] getData(self)
    cdef vector[DiffExprIsland]* getVectorPtr(self, string chrom) nogil
