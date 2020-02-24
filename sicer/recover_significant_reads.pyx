# SICER Internal Imports
from sicer.utility.utils cimport remove_at, get_tag_pos, bin_tag_in_island
from sicer.shared.data_classes cimport BEDRead, Island
from sicer.shared.chrom_containers cimport BEDReadContainer, IslandContainer

# Cython Imports
from libc.math cimport log, fmin
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libc.stdint cimport uint32_t

cdef void _recover_significant_reads_by_chrom(
    vector[Island]& islands,
    vector[BEDRead]& treatment_reads,
    int frag_size
) nogil:
    if islands.size() == 0 or treatment_reads.size() == 0:
        return

    cdef uint32_t pos
    cdef int index
    cdef vector[uint32_t] island_starts = vector[uint32_t](islands.size())
    cdef vector[uint32_t] island_ends = vector[uint32_t](islands.size())
    cdef vector[uint32_t] remove
    for i in range(islands.size()):
        island_starts[i] = islands[i].start
        island_ends[i] = islands[i].end

    for i in range(treatment_reads.size()):
        pos = get_tag_pos(treatment_reads[i], frag_size)
        index = bin_tag_in_island(island_starts, island_ends, pos)
        if index >= 0:
            remove.push_back(i)

    treatment_reads.erase(
        remove_at(treatment_reads.begin(), treatment_reads.end(), remove.begin(), remove.end()), 
        treatment_reads.end()
    )

cdef BEDReadContainer _recover_significant_reads(
    IslandContainer islands,
    BEDReadContainer treatment_reads,
    int frag_size,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = islands.getChromosomes()
    cdef int i
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _recover_significant_reads_by_chrom(
            deref(islands.getVectorPtr(chroms.at(i))),
            deref(treatment_reads.getVectorPtr(chroms.at(i))),
            frag_size
        )

    treatment_reads.updateReadCount()

    return treatment_reads

cpdef BEDReadContainer recover_significant_reads(
    islands,
    treatment_reads,
    frag_size,
    num_cpu
):
    print("Filtering reads with identified significant islands...")
    return _recover_significant_reads(
        islands,
        treatment_reads,
        frag_size,
        num_cpu
    )


