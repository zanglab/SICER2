# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Window
from sicer.utility.utils cimport to_string
from sicer.shared.chrom_containers cimport ChromWindowContainer

#Cython Imports
from libc.stdio cimport FILE, fopen, snprintf, fprintf, fclose
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.operator cimport dereference as deref

ctypedef vector[Window]* win_vec_ptr

cdef class WigFileWriter:
    # Writes .wig file of windows

    def __cinit__(self, 
        str file_name, 
        str output_dir, 
        ChromWindowContainer windows,
        int window_size,
        bint filtered
    ):
        self.file_name = file_name
        self.output_dir = output_dir
        self.windows = windows
        self.window_size = window_size
        self.filtered = filtered

    cdef cstr format_line(self, int pos, double count):
        cdef char buffer[50]
        snprintf (buffer, 44, "%d\t%.2f\n", pos, count)
        return buffer

    cdef void c_write(self,
        cstr outfile_path,
        cstr header,
        int window_size,
        double scaling_factor
    ):
        cdef FILE *fp = fopen(outfile_path, "w")
        fprintf(fp, header)

        cdef vector[string] chroms = self.windows.getChromosomes()
        cdef string chrom_header
        cdef int pos
        cdef double norm_count
        cdef cstr line
        cdef win_vec_ptr vptr

        for i in range(chroms.size()):
            vptr = self.windows.getChromVector(chroms[i])
            if deref(vptr).size() > 0:
                chrom_header = b"variableStep chrom=" + chroms[i] + b" span=" + to_string(window_size) + b"\n"
                fprintf(fp, chrom_header.c_str())
                for j in range(deref(vptr).size()):
                    line = self.format_line(deref(vptr)[j].start + 1, deref(vptr)[j].count / scaling_factor)
                    fprintf(fp, line)

    cpdef void write(self):
        print("Normalizing graphs by total island filitered reads per million and generating summary WIG file...\n")
        # We first need to normalize
        cdef int count = self.windows.getTotalTagCount()
        cdef double scaling_factor =  count / 1000000.0 * (self.window_size / 1000.0)

        # Format final file_name
        self.file_name += "-W" + str(self.window_size)
        if self.filtered:
            self.file_name += "-islandfiltered"

        cdef bytes wig_header = ("track type=wiggle_0 name=" + self.file_name + "\n").encode("UTF-8")
        self.file_name += "-normalized.wig"

        cdef bytes outfile_path = (self.output_dir + "/" + self.file_name).encode("UTF-8")

        self.c_write(outfile_path, wig_header, self.window_size, scaling_factor)
