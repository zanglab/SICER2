from cython.operator cimport predecrement as predec
from libc.math cimport M_PI, log, exp, pow, round
from libcpp.algorithm cimport upper_bound, lower_bound

cdef int get_tag_pos(BEDRead read, int frag_size) nogil:
    cdef int shift = <int> round(frag_size / 2.0)
    if read.strand == b'+':
        return read.start + shift
    else:
        return read.end - shift - 1

cdef int bin_tag_in_island(
    vector[int] &island_starts,
    vector[int] &island_ends,
    int tag_pos
) nogil:
    cdef int start_index = upper_bound(island_starts.begin(), island_starts.end(), tag_pos) - island_starts.begin()
    cdef int end_index = lower_bound(island_ends.begin(), island_ends.end(), tag_pos) - island_ends.begin()
    if (start_index < island_starts.size() 
        and end_index < island_ends.size()
        and start_index - end_index == 1
    ):
        return start_index - 1;
    else:
        return -1;

cdef int fact(int n) nogil:
    cdef int val = 1
    if n != 0:
        while n != 1:
            val = val * n;
            predec(n)
    return val

# Return the log of a factorial, using Srinivasa Ramanujan's approximation when m>=20
cdef double factln(int n) nogil:
    if n < 20:
        return log(fact(n))
    else:
        return n * log(n) - n + log(n * (1 + 4 * n * (1 + 2 * n))) / 6.0 + log(M_PI) / 2

cpdef double poisson(int n, double avg) nogil:
    if n < 20:
        return exp(-avg) * pow(avg, n) / fact(n)
    else:
        return exp(-avg + n * log(avg) - factln(n))