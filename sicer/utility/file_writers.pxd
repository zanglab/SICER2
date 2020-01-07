# SICER Internal Imports
from sicer.shared.chrom_containers cimport ChromWindowContainer

ctypedef char* cstr

cdef class WigFileWriter:
    # Writes .wig file of windows
    cdef:
        str file_name
        str output_dir
        ChromWindowContainer windows
        int window_size
        bint filtered

        cstr format_line(self, int pos, double count)

        void c_write(self, 
            cstr outfile_path,
            cstr header, 
            int window_size,
            double scaling_factor
        )

    cpdef void write(self)

cdef class IslandSummaryFileWriter:
    # Writes .wig file of windows
    cdef:
        str file_name
        str output_dir
        ChromWindowContainer windows
        int window_size
        bint filtered

        cstr format_line(self, int pos, double count)

        void c_write(self, 
            cstr outfile_path,
            cstr header, 
            int window_size,
            double scaling_factor
        )

    cpdef void write(self)