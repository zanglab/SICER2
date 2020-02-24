# SICER Internal Imports
from sicer.shared.utils cimport get_tag_pos, bin_tag_in_island
from sicer.shared.data_classes cimport BEDRead, Island
from sicer.shared.containers cimport BEDReadContainer, IslandContainer

# Cython Imports
from libc.math cimport log, fmin
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libc.stdint cimport uint32_t

from scipy.special.cython_special cimport pdtrc as poisson_sf
from scipy.stats import rankdata

cdef vector[double] _associate_tag_count_to_regions_by_chrom (
    vector[Island]& islands,
    vector[BEDRead]& treatment_reads,
    vector[BEDRead]& control_reads,
    double genome_size,
    double scaling_factor,
    int frag_size,
    int ctrl_lib_size
) nogil:
    if islands.size() == 0:
        # Return empty vector
        return vector[double]()

    cdef uint32_t pos
    cdef int index
    cdef vector[uint32_t] island_starts = vector[uint32_t](islands.size())
    cdef vector[uint32_t] island_ends = vector[uint32_t](islands.size())
    for i in range(islands.size()):
        island_starts[i] = islands[i].start
        island_ends[i] = islands[i].end

    for i in range(treatment_reads.size()):
        pos = get_tag_pos(treatment_reads[i], frag_size)
        index = bin_tag_in_island(island_starts, island_ends, pos)
        if index >= 0:
            preinc(islands[index].obs_count)

    for i in range(control_reads.size()):
        pos = get_tag_pos(control_reads[i], frag_size)
        index = bin_tag_in_island(island_starts, island_ends, pos)
        if index >= 0:
            preinc(islands[index].control_count)

    cdef uint32_t length
    cdef double avg, fc, pvalue
    cdef vector[double] pvalue_vec = vector[double](islands.size())
    for i in range(islands.size()):
        if (islands[i].control_count > 0):
            avg = islands[i].control_count * scaling_factor
        else:
            length = (islands[i].end - islands[i].start + 1) * ctrl_lib_size
            avg = fmin(0.25, (length / genome_size)) * scaling_factor
        fc = islands[i].obs_count / avg
        if islands[i].obs_count > avg:
            pvalue = poisson_sf(islands[i].obs_count, avg)
        else:
            pvalue = 1

        pvalue_vec[i] = pvalue
        islands[i].pvalue = pvalue
        islands[i].fc = fc

    return pvalue_vec

cdef IslandContainer _associate_tag_count_to_region(
    IslandContainer islands,
    BEDReadContainer treatment_reads,
    BEDReadContainer control_reads,
    double genome_size,
    double scaling_factor,
    int frag_size,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = islands.getChromosomes()
    cdef int ctrl_lib_size = control_reads.getReadCount()
    cdef int i
    cdef vector[double] pvalues
    cdef vector[double] returned_vec
    pvalues.reserve(islands.getIslandCount())
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        returned_vec = _associate_tag_count_to_regions_by_chrom(
                            deref(islands.getVectorPtr(chroms.at(i))),
                            deref(treatment_reads.getVectorPtr(chroms.at(i))),
                            deref(control_reads.getVectorPtr(chroms.at(i))),
                            genome_size,
                            scaling_factor,
                            frag_size,
                            ctrl_lib_size
                        )
        # Concatenate returned vectors
        pvalues.insert(pvalues.end(), returned_vec.begin(), returned_vec.end())

    cdef list pvalues_list = pvalues
    cdef vector[double] pvalue_rank = rankdata(pvalues_list)
    cdef vector[Island]* vptr
    cdef int k = 0
    cdef double alpha_stat

    for i in range(chroms.size()):
        vptr = islands.getVectorPtr(chroms[i])
        for j in range(deref(vptr).size()):
            alpha_stat = pvalues[k] * pvalues.size() / pvalue_rank[k]
            if alpha_stat > 1:
                alpha_stat = 1
            deref(vptr)[j].alpha_stat = alpha_stat
            preinc(k)

    islands.updateIslandCount()

    return islands

cpdef IslandContainer associate_tags_with_control(
    islands,
    treatment_reads,
    control_reads, 
    genome_size, 
    scaling_factor, 
    frag_size,
    num_cpu
):
    print("Calculating significance of candidate islands using the control library...")
    return _associate_tag_count_to_region(
        islands,
        treatment_reads,
        control_reads,
        genome_size,
        scaling_factor,
        frag_size,
        num_cpu
    )
