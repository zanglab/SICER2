# SICER Internal Imports
from sicer.shared.data_classes cimport Island
from sicer.shared.chrom_containers cimport ChromIslandContainer
from sicer.utility.utils cimport merge

# Cython Imports
from libcpp cimport bool
from libc.stdint cimport uint32_t
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange

ctypedef bool (*cmp_f)(Island, Island)
ctypedef vector[Island].iterator vi_itr

cdef bool compare_islands(Island i, Island j) nogil:
    return i.start < j.start

cdef void _find_union_islands_by_chrom(
    vector[Island]& union_islands,
    vector[Island]& islands_1,
    vector[Island]& islands_2
) nogil:

    if islands_1.size() == 0:
        union_islands.swap(islands_2)
        return
    if islands_2.size() == 0:
        union_islands.swap(islands_1)
        return

    cdef vector[Island] merged_islands
    merge[vi_itr, vi_itr, vi_itr, cmp_f](
        islands_1.begin(), 
        islands_1.end(),
        islands_2.begin(),
        islands_2.end(),
        merged_islands.begin(),
        compare_islands
    )

    cdef Island current = merged_islands[0]
    cdef Island next;
    cdef uint32_t i = 1

    while i < merged_islands.size():
        next = merged_islands[i]
        if next.start > current.end:
            union_islands.push_back(current)
            current = next
        else:
            if next.end > current.end:
                current.end = next.end
        preinc(i)

    union_islands.push_back(current)

cdef ChromIslandContainer _find_union_islands(
    ChromIslandContainer islands_1,
    ChromIslandContainer islands_2,
    object genome_data,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = islands_1.getChromosomes()

    cdef ChromIslandContainer union_islands = ChromIslandContainer(genome_data)

    cdef int i
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _find_union_islands_by_chrom(
            deref(union_islands.getVectorPtr(chroms[i])),
            deref(islands_1.getVectorPtr(chroms[i])),
            deref(islands_2.getVectorPtr(chroms[i]))
        )

    union_islands.updateIslandCount()
    print("Union of islands count: ", union_islands.getIslandCount())

    return union_islands

cpdef ChromIslandContainer find_union_islands(
    islands_1,
    islands_2,
    genome_data,
    num_cpu
):
    print("Finding the union islands between two treatment libraries...")
    return _find_union_islands(
        islands_1,
        islands_2,
        genome_data,
        num_cpu
    )


