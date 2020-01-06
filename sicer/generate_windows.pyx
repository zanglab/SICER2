# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Window
from sicer.shared.chrom_containers cimport ChromBEDReadContainer, ChromWindowContainer

# Cython Imports
from libc.math cimport round as roundcpp
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

cdef char PLUS = b'+'

cdef vector[int] _get_tag_list(bed_vec_ptr reads, int chrom_length, int frag_size) nogil:
    cdef vector[int] tag_list
    cdef int shift = <int> roundcpp(frag_size / 2.0)
    cdef int pos

    for i in range(deref(reads).size()):
        read = deref(reads)[i]
        if read.start >= 0 and read.end < chrom_length:
            if read.strand == PLUS:
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
    win_vec_ptr windows,
    vector[int] tag_list,
    string chrom,
    int chrom_length,
    int window_size
) nogil:
    # Modify window vector in place
    cdef int curr_win_start
    cdef int curr_win_end
    cdef int curr_tag_count
    cdef int adjusted_tag_pos
    
    if tag_list.size() > 0:
        curr_win_start = (tag_list.at(0) // window_size) * window_size
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
                        deref(windows).push_back(Window(
                            chrom, 
                            curr_win_start, 
                            curr_win_end, 
                            curr_tag_count
                        ))

                    curr_win_start = adjusted_tag_pos
                    curr_tag_count = 1
        
        curr_win_end = curr_win_start + window_size - 1
        if curr_win_end < chrom_length:
            deref(windows).push_back(Window(
                chrom, 
                curr_win_start, 
                curr_win_end, 
                curr_tag_count
            ))
    
cdef void _generate_windows_by_chrom(
    bed_vec_ptr reads, 
    win_vec_ptr windows,
    string chrom, 
    int chrom_length,
    int frag_size,
    int window_size
) nogil:
    cdef vector[int] tag_list
    if deref(reads).size() > 0:
        tag_list = _get_tag_list(reads, chrom_length, frag_size)
        _generate_window_from_tags(windows, tag_list, chrom, chrom_length, window_size)
    
cdef ChromWindowContainer _generate_windows(
    ChromBEDReadContainer reads, 
    object genome_data, 
    int frag_size,
    int window_size,
    int num_cpu
):
    # Convert Python list to vector for no-GIL use
    cdef vector[string] chroms = reads.getChromosomes()
    cdef vector[int] chrom_lengths
    for c in genome_data.chrom:
        chrom_lengths.push_back(genome_data.chrom_length[c])

    cdef ChromWindowContainer windows = ChromWindowContainer(genome_data)

    cdef int i
    for i in prange(chroms.size(), schedule='guided', num_threads=num_cpu, nogil=True):
        _generate_windows_by_chrom(
            reads.getChromVector(chroms.at(i)),
            windows.getChromVector(chroms.at(i)),
            chroms.at(i),
            chrom_lengths.at(i),
            frag_size,
            window_size
        )

    windows.updateCounts()
    print("Window count: ", windows.getWindowCount())

    return windows

cpdef ChromWindowContainer generate_windows(reads, genome_data, frag_size, window_size, num_cpu):
    print("Generating windows from treatement reads...")
    return _generate_windows(reads, genome_data, frag_size, window_size, num_cpu)


