from libcpp.vector cimport vector
from libcpp.map cimport map as mapcpp
from libcpp.string cimport string

ctypedef char* cstr

# Represents a single line of BED file
cdef extern from "dataStructs.h":
    cdef cppclass BEDRead:
        string chrom
        int start
        int end
        string name
        int score
        char strand

        BEDRead()
        BEDRead(string chrom, int start, int end, string name, int score, char strand)
        string toString()

