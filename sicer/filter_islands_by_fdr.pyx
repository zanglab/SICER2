# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Island
from sicer.utility.utils cimport remove_at
from sicer.shared.chrom_containers cimport ChromBEDReadContainer, ChromWindowContainer, ChromIslandContainer

# Cython Imports
from libc.math cimport log, fmin, isnan
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libc.stdint cimport uint32_t

from scipy.special.cython_special cimport pdtrc as poisson_sf
from scipy.stats import rankdata

cdef void _filter_islands_by_fdr_by_chrom(vector[Island]& islands, double fdr) nogil:

    cdef vector[int] remove
    if islands.size() > 0:
        for i in range(islands.size()):
            if islands[i].alpha_stat > fdr:
                remove.push_back(i)

        islands.erase(
            remove_at(islands.begin(), islands.end(), remove.begin(), remove.end()),
            islands.end()
        )

cdef ChromIslandContainer _filter_islands_by_fdr(
    ChromIslandContainer islands,
    double fdr,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = islands.getChromosomes()
    cdef int i

    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _filter_islands_by_fdr_by_chrom(deref(islands.getVectorPtr(chroms.at(i))), fdr)

    islands.updateIslandCount()

    return islands

cpdef ChromIslandContainer filter_islands_by_fdr(islands, fdr, num_cpu):
    print("Identifying significant islands using FDR criterion...")
    return _filter_islands_by_fdr(islands, fdr, num_cpu)
