from libcpp.string cimport string

cdef extern from "<string>" namespace "std" nogil:
    string to_string(int val)

cdef extern from "removeAt.cpp" nogil:
    V remove_at[V, I](V first, V last, I ii_first, I ii_last)