# SICER Internal Imports
from sicer.shared.utils cimport get_tag_pos, bin_tag_in_island, remove_at
from sicer.shared.data_classes cimport BEDRead, Island, DiffExprIsland
from sicer.shared.containers cimport BEDReadContainer, IslandContainer, DiffExprIslandContainer

# Cython Imports
from libc.math cimport log, fmin, isnan
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libc.stdint cimport uint32_t

from scipy.special.cython_special cimport pdtrc as poisson_sf
import scipy.stats

cdef struct ReturnItem:
    vector[double] pvalues_A_vs_B
    vector[double] pvalues_B_vs_A

cdef double _calc_pvalue(
    uint32_t count_A,
    uint32_t count_B,
    double scaling_factor,
) nogil:
    """
    Currently using poisson distribution

    scaling_factor: the factor that accounts for the differences of control library and ChIP library. effective control read count
    is control_read_count * scaling factor
    pseudocount: when control_read_count is zero, replace zero with pseudocount to alleviate the impact of statistical fluctuation

    output: pvalue
    """
    cdef double average
    cdef double pvalue
    if count_B > 0:
        average = count_B * scaling_factor
    else:
        average = scaling_factor
    if count_A > average:
        pvalue = poisson_sf(count_A, average)
    else:
        pvalue = 1
    return pvalue

cdef ReturnItem _associate_tag_count_to_regions_by_chrom (
    vector[DiffExprIsland]& islands,
    vector[BEDRead]& reads_A,
    vector[BEDRead]& reads_B,
    double scaling_factor,
    int frag_size
) nogil:
    cdef ReturnItem ret

    if islands.size() == 0:
        return ret

    cdef uint32_t pos
    cdef int index
    cdef vector[uint32_t] island_starts = vector[uint32_t](islands.size())
    cdef vector[uint32_t] island_ends = vector[uint32_t](islands.size())
    for i in range(islands.size()):
        island_starts[i] = islands[i].start
        island_ends[i] = islands[i].end

    for i in range(reads_A.size()):
        pos = get_tag_pos(reads_A[i], frag_size)
        index = bin_tag_in_island(island_starts, island_ends, pos)
        if index >= 0:
            preinc(islands[index].count_A)

    for i in range(reads_B.size()):
        pos = get_tag_pos(reads_B[i], frag_size)
        index = bin_tag_in_island(island_starts, island_ends, pos)
        if index >= 0:
            preinc(islands[index].count_B)

    # We store pvalues again in vectors to easily compute fdr later
    cdef vector[double] pvalues_A_vs_B = vector[double](islands.size())
    cdef vector[double] pvalues_B_vs_A = vector[double](islands.size())
    cdef double pvalue_A_vs_B, pvalue_B_vs_A 

    for i in range(islands.size()):
        pvalue_A_vs_B = _calc_pvalue(islands[i].count_A, islands[i].count_B, scaling_factor)
        pvalue_B_vs_A = _calc_pvalue(islands[i].count_B, islands[i].count_A, 1/scaling_factor)
        pvalues_A_vs_B[i] = pvalue_A_vs_B
        pvalues_B_vs_A[i] = pvalue_B_vs_A
        islands[i].pvalue_A_vs_B = pvalue_A_vs_B
        islands[i].pvalue_B_vs_A = pvalue_B_vs_A

    ret.pvalues_A_vs_B.swap(pvalues_A_vs_B)
    ret.pvalues_B_vs_A.swap(pvalues_B_vs_A)

    return ret

cpdef DiffExprIslandContainer compare_two_libraries(
    object genome_data,
    BEDReadContainer reads_A, 
    BEDReadContainer reads_B,
    DiffExprIslandContainer union_islands,
    int frag_size,
    int num_cpu
):
    print("Comparing two treatment libraries...")
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = reads_A.getChromosomes()
    cdef uint32_t lib_size_A = reads_A.getTotalReadCount()
    cdef uint32_t lib_size_B = reads_B.getTotalReadCount()
    cdef double lib_scaling_factor = (<double> lib_size_A) / lib_size_B

    cdef int i
    cdef ReturnItem result
    cdef vector[double] pvalues_A_vs_B
    pvalues_A_vs_B.reserve(union_islands.getIslandCount())
    cdef vector[double] pvalues_B_vs_A
    pvalues_B_vs_A.reserve(union_islands.getIslandCount())

    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        result = _associate_tag_count_to_regions_by_chrom(
                            deref(union_islands.getVectorPtr(chroms.at(i))),
                            deref(reads_A.getVectorPtr(chroms.at(i))),
                            deref(reads_B.getVectorPtr(chroms.at(i))),
                            lib_scaling_factor,
                            frag_size
                        )
        pvalues_A_vs_B.insert(pvalues_A_vs_B.end(), result.pvalues_A_vs_B.begin(), result.pvalues_A_vs_B.end())
        pvalues_B_vs_A.insert(pvalues_B_vs_A.end(), result.pvalues_B_vs_A.begin(), result.pvalues_B_vs_A.end())

    union_islands.updateIslandCount()

    pvalues_rank_A_vs_B = scipy.stats.rankdata(pvalues_A_vs_B)
    pvalues_rank_B_vs_A = scipy.stats.rankdata(pvalues_B_vs_A)

    cdef uint32_t total_count = pvalues_A_vs_B.size()
    cdef vector[double] norm_counts_A = vector[double](total_count)
    cdef vector[double] norm_counts_B = vector[double](total_count)

    cdef uint32_t k = 0
    cdef vector[DiffExprIsland] vec
    cdef double norm_count_A, norm_count_B

    for chrom in union_islands.getChromosomes():
        vecptr = union_islands.getVectorPtr(chrom)
        for j in range(deref(vecptr).size()):
            deref(vecptr)[j].fdr_A_vs_B = fmin(1.0, pvalues_A_vs_B[k] * total_count / pvalues_rank_A_vs_B[k])
            deref(vecptr)[j].fdr_B_vs_A = fmin(1.0, pvalues_B_vs_A[k] * total_count / pvalues_rank_B_vs_A[k])

            norm_count_A = (<double> deref(vecptr)[j].count_A) / (<double> lib_size_A) * 1000000
            norm_count_B = (<double> deref(vecptr)[j].count_B) / (<double> lib_size_B) * 1000000

            deref(vecptr)[j].norm_count_A = norm_count_A
            deref(vecptr)[j].norm_count_B = norm_count_B
            norm_counts_A[k] = norm_count_A
            norm_counts_B[k] = norm_count_B
            deref(vecptr)[j].fc_A_vs_B = ((<double> (deref(vecptr)[j].count_A + 1)) / (deref(vecptr)[j].count_B + 1)) / lib_scaling_factor
            deref(vecptr)[j].fc_B_vs_A = ((<double> (deref(vecptr)[j].count_B + 1)) / (deref(vecptr)[j].count_A + 1)) * lib_scaling_factor
            preinc(k)

    pearson = scipy.stats.pearsonr(norm_counts_A, norm_counts_B)
    print("Pearson's correlation is: ", pearson[0], " with p-value ", pearson[1])
    spearman = scipy.stats.spearmanr(norm_counts_A, norm_counts_B)
    print("Spearman's correlation is: ", spearman[0], " with p-value ", spearman[1])

    return union_islands
