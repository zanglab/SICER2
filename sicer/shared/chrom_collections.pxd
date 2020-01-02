# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead

# Cython Imports
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef char* cstr

cdef class ChromBEDReadCollection:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vector[BEDRead]] data

        void insertRead(self, string chrom, BEDRead input_data)
        mapcpp[string, vector[BEDRead]] getData(self)
        vector[BEDRead] getChromData(self, string chrom)
        BEDRead getRead(self, string chrom, int i)
    
    cpdef list getChromosomes(self)
    cpdef string dataToString(self)
    cpdef void printDataHead(self)