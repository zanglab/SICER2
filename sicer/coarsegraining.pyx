# SICER Internal Imports
from sicer.shared.data_classes cimport Window, Island
from sicer.shared.chrom_containers cimport ChromWindowContainer, ChromIslandContainer
from sicer.utility.utils cimport merge

# Cython Imports
from libc.math cimport log, fmax, fabs, pow
from libcpp cimport bool
from libc.stdint cimport uint32_t
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from cython.parallel import parallel, prange
from libcpp.algorithm cimport binary_search, sort

ctypedef pair[vector[uint32_t], vector[double]] corr_pair
ctypedef bool (*cmp_f)(Island, Island)
ctypedef vector[Island].iterator vi_itr
ctypedef vector[uint32_t]* v_uint_ptr

cdef bool compare_islands(Island i, Island j) nogil:
    return i.start < j.start

cdef vector[Island] _generate_islands_from_start_pos(
    vector[uint32_t]& start_list,
    int window_size,
    string chrom
) nogil:

    cdef vector[Island] islands
    for i in range(start_list.size()):
        islands.push_back(Island(chrom, start_list[i], start_list[i] + window_size - 1, 1))

    return islands

cdef vector[Island] _backstep(
    vector[Island]& islands,
    vector[uint32_t] start_list,
    int window_size,
    string chrom
) nogil:
    cdef bool start_left, start_right, end_left, end_right
    cdef vector[Island] merged_islands, additional_islands, kept_islands
    additional_islands = _generate_islands_from_start_pos(start_list, window_size, chrom)

    for i in range(islands.size()):
        start_left = binary_search(start_list.begin(), start_list.end(), (islands[i].start - window_size))
        start_right = binary_search(start_list.begin(), start_list.end(), islands[i].start)

        if start_left and start_right:
            islands[i].start = islands[i].start - window_size
        elif not start_left and not start_right:
            islands[i].start = islands[i].start + window_size

        end_left = binary_search(start_list.begin(), start_list.end(), islands[i].end + 1 - window_size)
        end_right = binary_search(start_list.begin(), start_list.end(), islands[i].end + 1)

        if end_left and end_right:
            islands[i].end = islands[i].end + window_size
        elif not end_left and not end_right:
            islands[i].end = islands[i].end - window_size

    sort[vi_itr, cmp_f](islands.begin(), islands.end(), compare_islands)

    merged_islands = vector[Island](islands.size() + additional_islands.size())
    merge[vi_itr, vi_itr, vi_itr, cmp_f](
        islands.begin(), 
        islands.end(),
        additional_islands.begin(),
        additional_islands.end(),
        merged_islands.begin(),
        compare_islands
    )

    kept_islands = vector[Island]()

    cdef uint32_t current_start = merged_islands[0].start
    cdef uint32_t current_end = merged_islands[0].end
    cdef uint32_t j = 1
    cdef Island next_island
    while j < merged_islands.size():
        next_island = merged_islands[j]
        if next_island.start > current_end + window_size + 1:
            kept_islands.push_back(Island(chrom, current_start, current_end, 1))
            current_start = next_island.start
            current_end = next_island.end
        else:
            if current_end < next_island.end: 
                current_end = next_island.end

        preinc(j)
    kept_islands.push_back(Island(chrom, current_start, current_end, 1))

    return kept_islands

cdef double _linreg(vector[uint32_t]& x, vector[double]& y) nogil:
    
    cdef int N = x.size()
    cdef double Sx, Sy, Sxx, Syy, Sxy
    Sx = Sy = Sxx = Syy = Sxy = 0.0

    for i in range(N):
        Sx += x[i]
        Sy += y[i]
        Sxx += x[i] * x[i]
        Syy += y[i] * y[i]
        Sxy += x[i] * y[i]
    det = Sxx * N - Sx * Sx
    if det != 0:
        return (Sxy * N - Sy * Sx)/det
    else:
        return 0.0

cdef double _start_list_correlation_r_rev(
    vector[uint32_t]& start_list,
    int win, 
    int r, 
    int chrom_length
) nogil:
    cdef int x, d, n, i, SUMM
    cdef double s
    cdef vector[uint32_t] a
    x = start_list[0] % win
    d = r // win
    SUMM = 0
    n = (chrom_length - x) // win
    if n - d > 0:
        a = vector[uint32_t](n)
        for j in range(start_list.size()):
            i = (start_list[j] - x) // win
            if i >= 0 and i < n:
                a[i] = 1
        for i in range(0, n - d):
            SUMM += a[i] * a[i + d]
        s = 0.0
        for i in range(a.size()):
            s += a[i]
        s = s / a.size()
        return SUMM / (n - d) - (s * s)
    else:
        return 0.0

cdef corr_pair _start_list_correlation_function(
    vector[uint32_t]& start_list, 
    int win, 
    int chrom_length
) nogil:
    cdef vector[uint32_t] xlist
    cdef vector[double] ylist
    cdef int max_iter = chrom_length // win
    cdef uint32_t r
    cdef double c
    if max_iter > 3:
        max_iter = 3

    cdef uint32_t i = 0
    for i in range(max_iter):
        r = i * win
        c = _start_list_correlation_r_rev(start_list, win, r, chrom_length)
        xlist.push_back(i)
        ylist.push_back(c)

    return corr_pair(xlist, ylist)

cdef double _correlation_length_fit(corr_pair pair) nogil:
    cdef vector[uint32_t] xlist = pair.first
    cdef vector[double] ylist = pair.second

    cdef vector[double] loglist
    for i in range(ylist.size()):
        loglist.push_back(log(fmax(ylist[i], 1e-12)))
    
    if xlist.size() > 0:
        xlist.erase(xlist.begin())
        loglist.erase(loglist.begin())
    a = _linreg(xlist,loglist)
    if fabs(a) > 1e-12:
        return -1.0/a
    else:
        return 1e12

cdef vector[Island] _traceback(
    vector[v_uint_ptr] graining_results,
    int window_size,
    int step_size,
    uint32_t chrom_length,
    string chrom
) nogil:
    window_size = window_size * <int> pow(step_size, (graining_results.size() - 1))
    cdef uint32_t end = graining_results.size() - 1
    cdef vector[uint32_t] backlist = deref(graining_results[end])
    cdef corr_pair correlation_pair = _start_list_correlation_function(backlist, window_size, chrom_length)
    cdef double correlation_length = _correlation_length_fit(correlation_pair)
    cdef double correlation_length_next

    cdef int i = 1
    while i < graining_results.size():
        backlist = deref(graining_results[end - i])
        window_size = window_size // step_size
        correlation_pair = _start_list_correlation_function(backlist, window_size, chrom_length)
        correlation_length_next = _correlation_length_fit(correlation_pair)

        if correlation_length > 1.0 and correlation_length_next >= correlation_length:
            break
        else:
            correlation_length = correlation_length_next

        preinc(i)
    cdef vector[Island] islands = _generate_islands_from_start_pos(backlist, window_size, chrom)
    
    while i < graining_results.size():
        backlist = deref(graining_results[end - i])
        islands = _backstep(islands, backlist, window_size, chrom)
        window_size = window_size // step_size
        preinc(i)

    return islands

cdef v_uint_ptr _graining(vector[uint32_t]& starts, int window_size, int step_size, int step_score) nogil:

    cdef int i, j, h, k, n  # variables used by legacy code
    cdef uint32_t end_limit = starts[starts.size()-1]
    cdef v_uint_ptr new_starts = new vector[uint32_t]()

    cdef int step
    for step in range(step_size):
        i = starts[0] - step * window_size
        k = 0
        while i <= end_limit and k < starts.size():
            j = i + step_size * window_size
            h = k
            while h <= starts.size() - 1 and starts[h] < j:
                preinc(h)
            n = h - k
            if n >= step_score:
                deref(new_starts).push_back(i)
            k = h
            i = j

    return new_starts

cdef vector[v_uint_ptr] _coarsegraining(
    v_uint_ptr eligible_starts,
    int window_size,
    int step_size,
    int step_score
) nogil:
    
    cdef vector[v_uint_ptr] graining_results
    graining_results.push_back(eligible_starts)

    while deref(eligible_starts).size() > 0:
        graining_results.push_back(eligible_starts)

        # `graining` modifies `eligible_starts` inplace
        eligible_starts = _graining(deref(eligible_starts), window_size, step_size, step_score)

        window_size = window_size * step_size

    return graining_results

cdef void _coarsegraining_islands_by_chrom(
    vector[Window]& windows,
    vector[Island]& islands,
    uint32_t chrom_length,
    int min_tag_threshold,
    int window_size,
    int step_size,
    int step_score
) nogil:  
    cdef vector[uint32_t] eligible_starts
    cdef vector[Island] tracebacked_islands
    cdef string chrom
    cdef vector[v_uint_ptr] graining_results

    if windows.size() > 0:
        chrom = windows[0].chrom
        for i in range(windows.size()):
            if windows[i].count >= min_tag_threshold:
                eligible_starts.push_back(windows[i].start)

        # In-place modify islands vector
        if eligible_starts.size() > 0:
            graining_results = _coarsegraining(&eligible_starts, window_size, step_size, step_score)
            tracebacked_islands = _traceback(graining_results, window_size, step_size, chrom_length, chrom)
        # Clear vector[uint32_t] that are allocated in heap
        for i in range(graining_results.size()):
            graining_results[i].clear()
            graining_results[i].shrink_to_fit()

    islands.swap(tracebacked_islands)

cdef ChromIslandContainer _find_islands_by_coarsegraining(
    ChromWindowContainer windows,
    object genome_data,
    int min_tag_threshold,
    int window_size,
    int step_size,
    int step_score,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = windows.getChromosomes()

    cdef ChromIslandContainer islands = ChromIslandContainer(genome_data)
    cdef vector[uint32_t] chrom_lengths
    for c in genome_data.chrom:
        chrom_lengths.push_back(genome_data.chrom_length[c])

    cdef int i
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _coarsegraining_islands_by_chrom(
            deref(windows.getVectorPtr(chroms[i])),
            deref(islands.getVectorPtr(chroms[i])),
            chrom_lengths[i],
            min_tag_threshold,
            window_size,
            step_size,
            step_score
        )

    islands.updateIslandCount()
    print("Island count: ", islands.getIslandCount())

    return islands

cpdef ChromIslandContainer find_islands_by_coarsegraining(
    windows,
    genome_data,
    min_tag_threshold,
    window_size,
    step_size,
    step_score,
    num_cpu
):
    print("Finding candidate islands exhibiting clustering by coarsegraining...")
    return _find_islands_by_coarsegraining(
        windows,
        genome_data,
        min_tag_threshold,
        window_size,
        step_size,
        step_score,
        num_cpu
    )