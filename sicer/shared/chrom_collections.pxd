# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead

# Cython Imports
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef char* cstr
ctypedef vector[BEDRead]* vec_ptr

cdef class ChromBEDReadCollection:
    cdef:
        str species
        list chromosomes
        mapcpp[string, vec_ptr] data

        void insertRead(self, string chrom, BEDRead input_data)
        mapcpp[string, vec_ptr] getData(self)
        vec_ptr getChromData(self, string chrom)
        BEDRead getRead(self, string chrom, int i)
    
    cpdef list getChromosomes(self)
    cpdef string dataToString(self)
    cpdef void printDataHead(self)