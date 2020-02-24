# SICER Internal Imports
from sicer.shared.data_classes cimport Island, BEDRead
from sicer.shared.containers cimport BEDReadContainer, WindowContainer, IslandContainer, DiffExprIslandContainer

from libc.stdint cimport uint32_t

ctypedef char* cstr

cdef class WigFileWriter:
    # Writes .wig file of windows
    cdef:
        str file_name
        str output_dir
        WindowContainer windows
        int window_size
        bint filtered
        object fdr

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
        IslandContainer islands
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
        BEDReadContainer reads
        int window_size
        object gap_size
        object fdr

        void c_write(self, cstr outfile_path)

    cpdef void write(self)


cdef class DiffExprIslandWriter:
    # Writes islands produced from differential expression analysis
    cdef public:
        str file_name_1
        str file_name_2
        str output_dir
        DiffExprIslandContainer islands
        int window_size
        bint fdr_filtered
        bint increased
        object fdr
        object gap_size
        cstr header
        cstr format

        void c_write(self, cstr outfile_path)

    cpdef void write(self)