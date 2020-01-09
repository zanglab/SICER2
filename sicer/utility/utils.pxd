from sicer.shared.data_classes cimport BEDRead
from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "<string>" namespace "std" nogil:
    string to_string(int val)

cdef extern from "removeAt.cpp" nogil:
    V remove_at[V, I](V first, V last, I ii_first, I ii_last)

# Read interval manipulation functions
cdef int get_tag_pos(BEDRead read, int frag_size) nogil
cdef int bin_tag_in_island(vector[int]& island_starts, vector[int]& island_ends, int tag_pos) nogil

# Math functions
cdef int fact(int n) nogil
cdef double factln(int n) nogil
cpdef double poisson(int n, double avg) nogil