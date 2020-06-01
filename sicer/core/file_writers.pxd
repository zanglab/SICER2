# SICER Internal Imports
from sicer.shared.data_classes cimport Island, BEDRead
from sicer.shared.containers cimport BEDReadContainer, WindowContainer, IslandContainer, DiffExprIslandContainer

from libc.stdint cimport uint32_t
from libc.stdio cimport FILE

ctypedef char* cstr

cdef class WigFileWriter:
    # Writes .wig file of windows
    cdef:
        str file_name
        str output_dir
        WindowContainer windows
        int window_size
        bint filtered
        object gap_size
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
        str file_type
        IslandContainer islands
        int window_size
        object gap_size
        object fdr

        void write_summary_line(self, Island island, FILE *fp)
        void write_bed_line(self, Island island, FILE *fp)
        void write_scoreisland_line(self, Island island, FILE *fp)

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


cdef class DiffExprIslandFileWriter:
    # Writes islands produced from differential expression analysis
    cdef public:
        str file_name_1
        str file_name_2
        str output_dir
        str file_type
        DiffExprIslandContainer islands
        int window_size
        object e_value
        object fdr
        object gap_size
        cstr header
        cstr format

        void c_write_all(self, cstr outfile_path)
        void c_write_basic(self, cstr outfile_path)

    cpdef void write(self)