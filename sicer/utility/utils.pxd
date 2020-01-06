from libcpp.string cimport string
from libcpp.vector cimport vector

cdef extern from "<string>" namespace "std" nogil:
    string to_string(int val)

cdef extern from "removeAt.cpp" nogil:
    V remove_at[V, I](V first, V last, I ii_first, I ii_last)

# Math functions
cdef int fact(int n) nogil
cdef double factln(int n) nogil
cdef double poisson(int n, double avg) nogil