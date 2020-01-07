# TODO: Using a templated class would be nice, but Cython doesn't support it currently

# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Window, Island

# Cython Imports
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef char* cstr

cdef class ChromBEDReadContainer:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vector[BEDRead]] data
        int read_count

    cdef void insertRead(self, string chrom, BEDRead item)
    cpdef void updateReadCount(self)

    cpdef int getReadCount(self)
    cpdef list getChromosomes(self)
    cdef mapcpp[string, vector[BEDRead]] getData(self)
    cdef vector[BEDRead]* getVectorPtr(self, string chrom) nogil
    cdef BEDRead getRead(self, string chrom, int index)


cdef class ChromWindowContainer:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vector[Window]] data
        int window_count
        int total_tag_count

    cdef void insertWindow(self, string chrom, Window item)
    cpdef void updateCounts(self)

    cpdef int getWindowCount(self)
    cpdef list getChromosomes(self)
    cpdef int getTotalTagCount(self)
    cdef mapcpp[string, vector[Window]] getData(self)
    cdef vector[Window]* getVectorPtr(self, string chrom) nogil
    cdef Window getWindow(self, string chrom, int index)


cdef class ChromIslandContainer:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vector[Island]] data
        int island_count

    cdef void insertIsland(self, string chrom, Island item)
    cpdef void updateIslandCount(self)

    cpdef int getIslandCount(self)
    cpdef list getChromosomes(self)
    cdef mapcpp[string, vector[Island]] getData(self)
    cdef vector[Island]* getVectorPtr(self, string chrom) nogil
    cdef Island getIsland(self, string chrom, int index)

