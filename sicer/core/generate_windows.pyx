# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Window
from sicer.shared.containers cimport BEDReadContainer, WindowContainer

# Cython Imports
from libcpp cimport bool
from libc.math cimport round
from libc.stdint cimport uint32_t
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libcpp.algorithm cimport sort


cdef vector[uint32_t] _get_tag_list(vector[BEDRead]& reads, uint32_t chrom_length, int frag_size) nogil:
    cdef vector[uint32_t] tag_list
    cdef int shift = <int> round(frag_size / 2.0)
    cdef uint32_t pos

    for i in range(reads.size()):
        read = reads[i]
        if read.start >= 0 and read.end < chrom_length:
            if read.strand == b'+':
                pos = read.start + shift
                if pos >= chrom_length:
                    pos = chrom_length - 1
            else:
                pos = read.end - shift - 1
                if pos < 0:
                    pos = 0

            tag_list.push_back(pos)

    sort(tag_list.begin(), tag_list.end())

    return tag_list

cdef void _generate_window_from_tags(
    vector[Window]& windows,
    vector[uint32_t] tag_list,
    string chrom,
    uint32_t chrom_length,
    int window_size
) nogil:
    # Modify window vector in place
    cdef uint32_t curr_win_start
    cdef uint32_t curr_win_end
    cdef uint32_t curr_tag_count
    cdef uint32_t adjusted_tag_pos
    
    if tag_list.size() > 0:
        curr_win_start = (tag_list[0] // window_size) * window_size
        curr_tag_count = 1

        if tag_list.size() > 1:
            for i in range(1, tag_list.size()):
                adjusted_tag_pos = (tag_list[i] // window_size) * window_size
                if adjusted_tag_pos == curr_win_start:
                    preinc(curr_tag_count)
                elif adjusted_tag_pos > curr_win_start:
                    curr_win_end = curr_win_start + window_size - 1

                    if curr_win_end < chrom_length:
                        # Create new window
                        windows.push_back(Window(
                            chrom, 
                            curr_win_start, 
                            curr_win_end, 
                            curr_tag_count
                        ))

                    curr_win_start = adjusted_tag_pos
                    curr_tag_count = 1
        
        curr_win_end = curr_win_start + window_size - 1
        if curr_win_end < chrom_length:
            windows.push_back(Window(
                chrom, 
                curr_win_start, 
                curr_win_end, 
                curr_tag_count
            ))
    
cdef void _generate_windows_by_chrom(
    vector[BEDRead]& reads, 
    vector[Window]& windows,
    string chrom, 
    uint32_t chrom_length,
    int frag_size,
    int window_size
) nogil:
    cdef vector[uint32_t] tag_list
    if reads.size() > 0:
        tag_list = _get_tag_list(reads, chrom_length, frag_size)
        _generate_window_from_tags(windows, tag_list, chrom, chrom_length, window_size)
    
cdef WindowContainer _generate_windows(
    BEDReadContainer reads, 
    object genome_data, 
    int frag_size,
    int window_size,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = reads.getChromosomes()
    cdef vector[uint32_t] chrom_lengths
    for c in genome_data.chrom:
        chrom_lengths.push_back(genome_data.chrom_length[c])

    cdef WindowContainer windows = WindowContainer(genome_data)

    cdef int i
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _generate_windows_by_chrom(
            deref(reads.getVectorPtr(chroms.at(i))),
            deref(windows.getVectorPtr(chroms.at(i))),
            chroms[i],
            chrom_lengths[i],
            frag_size,
            window_size
        )

    windows.updateCounts()
    print("Window count: ", windows.getWindowCount())

    return windows

cpdef WindowContainer generate_windows(reads, genome_data, frag_size, window_size, num_cpu):
    print("Generating windows from treatement reads...")
    return _generate_windows(reads, genome_data, frag_size, window_size, num_cpu)

