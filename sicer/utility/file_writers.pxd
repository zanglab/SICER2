# SICER Internal Imports
from sicer.shared.data_classes cimport Island, BEDRead
from sicer.shared.chrom_containers cimport ChromBEDReadContainer, ChromWindowContainer, ChromIslandContainer

from libc.stdint cimport uint32_t

ctypedef char* cstr

cdef class WigFileWriter:
    # Writes .wig file of windows
    cdef:
        str file_name
        str output_dir
        ChromWindowContainer windows
        int window_size
        bint filtered
        object fdr

        cstr format_line(self, uint32_t pos, double count)

        void c_write(self, 
            cstr outfile_path,
            cstr header, 
            int window_size,
            double scaling_factor
        )

    cpdef void write(self)


cdef class IslandFileWriter:
    cdef:
        str file_name
        str output_dir
        object file_type
        ChromIslandContainer islands
        int window_size
        object gap_size
        object fdr

        cstr format_summary_line(self, Island island)
        cstr format_bed_line(self, Island island)
        cstr format_scoreisland_line(self, Island island)

        void c_write(self, cstr outfile_path)

    cpdef void write(self)


cdef class BEDFileWriter:
    # Writes BED files
    cdef:
        str file_name
        str output_dir
        ChromBEDReadContainer reads
        int window_size
        object gap_size
        object fdr

        cstr format_read(self, BEDRead read)

        void c_write(self, cstr outfile_path)

    cpdef void write(self)