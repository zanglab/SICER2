from cython.operator cimport predecrement as predec
from libc.math cimport M_PI, log, exp, pow

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

cdef double poisson(int n, double avg) nogil:
    if n < 20:
        return exp(-avg) * pow(avg, n) / fact(n)
    else:
        return exp(-avg + n * log(avg) - factln(n))
