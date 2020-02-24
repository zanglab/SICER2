from sicer.shared.data_classes cimport BEDRead
from libcpp.string cimport string
from libcpp.vector cimport vector
from libc.stdint cimport uint32_t

cdef extern from "<string>" namespace "std" nogil:
    string to_string(int val)
    string to_string(uint32_t val)

cdef extern from "removeAt.h" nogil:
    V remove_at[V, I](V first, V last, I ii_first, I ii_last)

# Read interval manipulation functions
cdef int get_tag_pos(BEDRead& read, int frag_size) nogil
cdef int bin_tag_in_island(
    vector[uint32_t]& island_starts,
    vector[uint32_t]& island_ends, 
    uint32_t tag_pos
) nogil

# Math functions
cdef uint32_t fact(int n) nogil
cdef double factln(int n) nogil
cpdef double poisson(int n, double avg) nogil

cdef extern from "<algorithm>" namespace "std" nogil:
    O merge[I1, I2, O, C] (I1 first1, I1 last1, I2 first2, I2 last2, O result, C compare)