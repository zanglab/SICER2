# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Window, Island
from sicer.shared.containers cimport BEDReadContainer, WindowContainer, IslandContainer

# Cython Imports
from libc.math cimport log, exp, M_PI
from libcpp cimport bool
from libc.stdint cimport uint32_t
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libcpp.algorithm cimport sort


cdef int WINDOW_SIZE_BUFFER = 2


cdef double fact(int m) nogil:
    cdef double value = 1.0
    if m != 0:
        while m != 1:
            value = value * m
            m = m - 1
    return value


# Return the log of a factorial, using Srinivasa Ramanujan's approximation when m>=20
cdef double factln(int m) nogil:
    cdef double value = 1.0
    if m < 20:
        if m != 0:
            while m != 1:
                value = value * m;
                m = m - 1;
        return log(value);
    else:
        return m * log(m) - m + log(m * (1 + 4 * m * (1 + 2 * m))) / 6.0 + log(M_PI) / 2


cdef double poisson(int i, double average) nogil:
    if i < 20:
        return exp(-average) * average ** i / fact(i)
    else:
        return exp(-average + i * log(average) - factln(i))


cdef void _filter_by_threshold(vector[Island]& islands, double score_threshold) nogil:
    cdef vector[Island] filtered_islands
    if islands.size() > 0:
        for i in range(0, islands.size()):
            if islands[i].score >= (score_threshold - .0000000001):
                filtered_islands.push_back(islands[i])

        islands.swap(filtered_islands)


cdef void _combine_proximal_islands(vector[Island]& islands, int gap_size) nogil:
    cdef uint32_t proximal_island_dist = gap_size + WINDOW_SIZE_BUFFER
    # Create new vector of islands b/c we generally throw away majority of islands
    cdef vector[Island] final_islands
    cdef string chrom
    cdef uint32_t curr_start
    cdef uint32_t curr_end
    cdef double curr_score
    cdef int dist

    if islands.size() > 0:
        if islands.size() == 1:
            return 

        chrom = islands[0].chrom
        curr_start = islands[0].start
        curr_end = islands[0].end
        curr_score = islands[0].score
        for i in range(1, islands.size()):
            dist = islands[i].start - curr_end
            if chrom == b"chr2" and curr_start == 32916000:
                print(islands[i].score, dist)
            if dist <= proximal_island_dist:
                curr_end = islands[i].end
                curr_score += islands[i].score
            else:
                #if chrom == b"chr2" and curr_start == 32916000 and curr_end == 32916799:
                    #print(curr_score)
                final_islands.push_back(Island(chrom, curr_start, curr_end, curr_score))
                curr_start = islands[i].start
                curr_end = islands[i].end
                curr_score = islands[i].score
        
        final_islands.push_back(Island(chrom, curr_start, curr_end, curr_score))

        # Swap with our final islands
        islands.swap(final_islands)


cdef void _generate_islands(
    vector[Window]& windows,
    vector[Island]& islands,
    int min_tag_threshold, 
    double avg_tag_count
) nogil:
    cdef double score
    for i in range(windows.size()):
        score = -1
        if windows[i].count >= min_tag_threshold:
            prob = poisson(windows[i].count, avg_tag_count)

            if prob < 1e-250:
                # Outside of the scale, so take an arbitrary value
                score = 1000
            else:
                score = -log(prob)
            if score > 0:
                # Create Island
                islands.push_back(Island(windows[i].chrom, windows[i].start, windows[i].end, score))


cdef void _find_islands_by_chrom(
    vector[Window]& windows,
    vector[Island]& islands,
    int min_tag_threshold,
    double score_threshold,
    int gap_size,
    double avg_tag_count
) nogil:
    
    if windows.size() > 0:
        # First create base islands from windows
        _generate_islands(windows, islands, min_tag_threshold, avg_tag_count)
        _combine_proximal_islands(islands, gap_size)
        _filter_by_threshold(islands, score_threshold)


cdef IslandContainer _find_islands(
    WindowContainer windows,
    object genome_data,
    int min_tag_threshold, 
    double score_threshold,
    int gap_size,
    double avg_tag_count,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = windows.getChromosomes()

    cdef IslandContainer islands = IslandContainer(genome_data)

    cdef int i
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
    #for i in range(chroms.size()):
        _find_islands_by_chrom(
            deref(windows.getVectorPtr(chroms[i])),
            deref(islands.getVectorPtr(chroms[i])),
            min_tag_threshold,
            score_threshold,
            gap_size,
            avg_tag_count
        )

    islands.updateIslandCount()
    print("Total island count: ", islands.getIslandCount())

    return islands


cpdef IslandContainer find_islands(
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


