# SICER Internal Imports
from sicer.shared.data_classes cimport Island, DiffExprIsland
from sicer.shared.utils cimport remove_at
from sicer.shared.containers cimport IslandContainer, DiffExprIslandContainer

# Cython Imports
from libcpp cimport bool
from cython.operator cimport dereference as deref
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libc.stdint cimport uint32_t

cdef void _filter_islands_by_fdr_by_chrom(vector[Island]& islands, double fdr_threshold) nogil:

    # Items to remove
    cdef vector[size_t] remove
    if islands.size() > 0:
        for i in range(islands.size()):
            if islands[i].alpha_stat > fdr_threshold:
                remove.push_back(i)

        islands.erase(
            remove_at(islands.begin(), islands.end(), remove.begin(), remove.end()),
            islands.end()
        )

cdef void _filter_islands_by_A_vs_B_fdr_by_chrom(
    vector[DiffExprIsland]& islands, 
    double fdr_threshold
) nogil:
    
    # Items to remove 
    cdef vector[size_t] remove
    if islands.size() > 0:
        for i in range(islands.size()):
            if islands[i].fdr_A_vs_B > fdr_threshold:
                remove.push_back(i)

        islands.erase(
            remove_at(islands.begin(), islands.end(), remove.begin(), remove.end()),
            islands.end()
        )

cdef void _filter_islands_by_B_vs_A_fdr_by_chrom(
    vector[DiffExprIsland]& islands, 
    double fdr_threshold
) nogil:
    
    # Items to remove
    cdef vector[size_t] remove
    if islands.size() > 0:
        for i in range(islands.size()):
            if islands[i].fdr_B_vs_A > fdr_threshold:
                remove.push_back(i)

        islands.erase(
            remove_at(islands.begin(), islands.end(), remove.begin(), remove.end()),
            islands.end()
        )

cdef IslandContainer _filter_islands_by_fdr(
    IslandContainer islands,
    double fdr,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = islands.getChromosomes()
    cdef int i

    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _filter_islands_by_fdr_by_chrom(deref(islands.getVectorPtr(chroms.at(i))), fdr)

    islands.updateIslandCount()
    print("Significant island count: ", islands.getIslandCount())

    return islands

cdef DiffExprIslandContainer _filter_diff_expr_islands_by_fdr(
    bool A_vs_B,
    DiffExprIslandContainer islands,
    double fdr,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = islands.getChromosomes()
    cdef int i

    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        if A_vs_B:    
            _filter_islands_by_A_vs_B_fdr_by_chrom(deref(islands.getVectorPtr(chroms.at(i))), fdr)
        else:
            _filter_islands_by_B_vs_A_fdr_by_chrom(deref(islands.getVectorPtr(chroms.at(i))), fdr)

    islands.updateIslandCount()

    return islands

cpdef filter_islands_by_fdr(islands, fdr, num_cpu, diff_expr, A_vs_B=False):
    print("Identifying significant islands using FDR criterion...")
    if diff_expr:
        return _filter_diff_expr_islands_by_fdr(A_vs_B, islands, fdr, num_cpu)
    else:
        return _filter_islands_by_fdr(islands, fdr, num_cpu)
