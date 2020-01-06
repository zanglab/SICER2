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
ctypedef vector[Window]* win_vec_ptr
ctypedef vector[Island]* isl_vec_ptr

cdef char PLUS = b'+'
cdef int WINDOW_SIZE_BUFFER = 2

cdef void _filter_by_threshold(isl_vec_ptr islands, double score_threshold) nogil:
    cdef vector[Island] filtered_islands
    if deref(islands).size() > 0:
        for i in range(0, deref(islands).size()):
            if deref(islands)[i].score >= (score_threshold - .0000000001):
                filtered_islands.push_back(deref(islands)[i])

        deref(islands).swap(filtered_islands)

cdef void _combine_proximal_islands(isl_vec_ptr islands, int gap_size) nogil:
    cdef int proximal_island_dist = gap_size + WINDOW_SIZE_BUFFER
    # Create new vector of islands b/c we generally throw away majority of islands
    cdef vector[Island] final_islands
    cdef string chrom
    cdef int curr_start
    cdef int curr_end
    cdef double curr_score

    if deref(islands).size() > 0:
        if deref(islands).size() == 1:
            return 

        chrom = deref(islands)[0].chrom
        curr_start = deref(islands)[0].start
        curr_end = deref(islands)[0].end
        curr_score = deref(islands)[0].score
        for i in range(1, deref(islands).size()):
            dist = deref(islands)[i].start - curr_end
            if dist <= proximal_island_dist:
                curr_end = deref(islands)[i].end
                curr_score += deref(islands)[i].score
            else:
                final_islands.push_back(Island(chrom, curr_start, curr_end, curr_score))
                curr_start = deref(islands)[i].start
                curr_end = deref(islands)[i].end
                curr_score = deref(islands)[i].score
        
        final_islands.push_back(Island(chrom, curr_start, curr_end, curr_score))

        # Swap with our final islands
        deref(islands).swap(final_islands)

cdef void _generate_islands(
    win_vec_ptr windows,
    isl_vec_ptr islands, 
    int min_tag_threshold, 
    double avg_tag_count
) nogil:
    cdef Window window
    cdef double score
    for i in range(deref(windows).size()):
        window = deref(windows)[i]
        score = -1
        if window.count >= min_tag_threshold:
            prob = poisson(window.count, avg_tag_count)

            if prob < 1e-250:
                # Outside of the scale, so take an arbitrary value
                score = 1000
            else:
                score = -log(prob)
            if score > 0:
                # Create Island
                deref(islands).push_back(Island(window.chrom, window.start, window.end, score))

cdef void _find_islands_by_chrom(
    win_vec_ptr windows,
    isl_vec_ptr islands,
    int min_tag_threshold,
    double score_threshold,
    int gap_size,
    double avg_tag_count
) nogil:
    
    if deref(windows).size() > 0:
        # First create base islands from windows
        _generate_islands(windows, islands, min_tag_threshold, avg_tag_count)
        _combine_proximal_islands(islands, gap_size)
        _filter_by_threshold(islands, score_threshold)

cdef ChromIslandContainer _find_islands(
    ChromWindowContainer windows,
    object genome_data,
    int min_tag_threshold, 
    double score_threshold,
    int gap_size,
    double avg_tag_count,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = windows.getChromosomes()

    cdef ChromIslandContainer islands = ChromIslandContainer(genome_data)

    cdef int i
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _find_islands_by_chrom(
            windows.getChromVector(chroms.at(i)),
            islands.getChromVector(chroms.at(i)),
            min_tag_threshold,
            score_threshold,
            gap_size,
            avg_tag_count
        )

    islands.updateIslandCount()
    print("Island count: ", islands.getIslandCount())

    return islands

cpdef ChromIslandContainer find_islands(
    windows,
    genome_data,
    min_tag_threshold, 
    score_threshold, 
    gap_size, 
    avg_tag_count,
    num_cpu
):
    print("Finding candidate islands exhibiting clustering...")
    return _find_islands(
        windows,
        genome_data,
        min_tag_threshold, 
        score_threshold,
        gap_size,
        avg_tag_count,
        num_cpu
    )


