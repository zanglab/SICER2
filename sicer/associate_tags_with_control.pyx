# SICER Internal Imports
from sicer.utility.utils cimport poisson
from sicer.shared.data_classes cimport BEDRead, Window, Island
from sicer.shared.chrom_containers cimport ChromBEDReadContainer, ChromWindowContainer, ChromIslandContainer

# Cython Imports
from libc.math cimport log
from libcpp cimport bool
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libcpp.algorithm cimport sort

# Typedefs
ctypedef char* cstr
ctypedef vector[BEDRead]* bed_vec_ptr
ctypedef vector[Window]* win_vec_ptr
ctypedef vector[Island]* isl_vec_ptr

cdef char PLUS = b'+'


cdef void _associate_tag_count_to_regions_by_chrom()
    isl_vec_ptr islands,
    bed_vec_ptr treatment_reads,
    bed_vec_ptr control_reads,
    double genome_size,
    double scaling_factor,
    int frag_size
) nogil:
    
    if deref(islands).size() > 0:
        cdef vector[int] island_starts(deref(islands).size())
        cdef vector[int] island_ends(deref(islands).size())
        for i in range(deref(islands).size()):
            islands_starts[i] = deref(islands)[i].start
            islands_ends[i] = deref(islands)[i].end

        for i in range(deref(treatment_reads.size())):


cdef ChromIslandContainer _associate_tag_count_to_region(
    ChromIslandContainer islands,
    ChromBEDReadContainer treatment_reads,
    ChromBEDReadContainer control_reads,
    double genome_size
    double scaling_factor,
    int frag_size,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = windows.getChromosomes()

    cdef int i
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _find_islands_by_chrom(
            islands.getChromVector(chroms.at(i)),
            treatment_reads.getChromVector(chroms.at(i)),
            control_reads.getChromVector(chroms.at(i)),
            genome_size,
            scaling_factor
            frag_size
        )

    islands.updateIslandCount()

    return islands

cpdef ChromIslandContainer associate_tags_with_control(
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
        windows,
        genome_data,
        min_tag_threshold, 
        score_threshold,
        gap_size,
        avg_tag_count,
        num_cpu
    )


