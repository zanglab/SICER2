import os

# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Window, DiffExprIsland
from sicer.shared.utils cimport to_string

#Cython Imports
from libc.stdio cimport fopen, snprintf, fprintf, fclose, printf
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.operator cimport dereference as deref

ctypedef vector[Window]* win_vec_ptr
ctypedef void (*format_f)(IslandFileWriter, Island, FILE*)

cdef class WigFileWriter:
    # Writes .wig file of windows

    def __cinit__(self, 
        str file_name, 
        str output_dir, 
        WindowContainer windows,
        int window_size,
        bint filtered,
        object fdr = None,
        object gap_size = None
    ):
        self.file_name = file_name
        self.output_dir = output_dir
        self.windows = windows
        self.window_size = window_size
        self.filtered = filtered
        self.gap_size = gap_size
        self.fdr = fdr

        if filtered and fdr is None:
            raise ValueError("Missing FDR value")

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
        cdef uint32_t pos
        cdef double norm_count
        cdef win_vec_ptr vptr

        for i in range(chroms.size()):
            vptr = self.windows.getVectorPtr(chroms[i])
            if deref(vptr).size() > 0:
                chrom_header = b"variableStep chrom=" + chroms[i] + b" span=" + to_string(window_size) + b"\n"
                fprintf(fp, chrom_header.c_str())
                for j in range(deref(vptr).size()):
                    fprintf(fp, "%u\t%.2f\n", deref(vptr)[j].start + 1, deref(vptr)[j].count / scaling_factor)

    cpdef void write(self):
        print("Normalizing graphs by total island filitered reads per million and generating summary WIG file...\n")
        # We first need to normalize
        cdef int count = self.windows.getTotalTagCount()
        cdef double scaling_factor =  count / 1000000.0 * (self.window_size / 1000.0)

        cdef str base = self.file_name
        if self.filtered:
            base += '-islandfiltered'
        cdef bytes wig_header = ("track type=wiggle_0 name=" + base + "\n").encode("UTF-8")

        # Format final file_name
        self.file_name += "-W" + str(self.window_size)
        if self.filtered:
            if self.gap_size:
                self.file_name += "-G" + str(self.gap_size) 
            self.file_name += "-FDR" + str(self.fdr) + "-islandfiltered"

        self.file_name += "-normalized.wig"

        cdef bytes outfile_path = (self.output_dir + "/" + self.file_name).encode("UTF-8")

        self.c_write(outfile_path, wig_header, self.window_size, scaling_factor)


cdef class IslandFileWriter:
    # Writes island files

    def __cinit__(self, 
        str file_name,
        str output_dir,
        str file_type,
        IslandContainer islands,
        int window_size,
        object gap_size = None,
        object fdr = None
    ):
        self.file_name = file_name
        self.file_type = file_type
        self.output_dir = output_dir
        self.islands = islands
        self.window_size = window_size
        self.gap_size = gap_size
        self.fdr = fdr

        if file_type == "fdr-filtered" and fdr is None:
            raise ValueError("Missing FDR value")

    cdef void write_summary_line(self, Island island, FILE *fp):
        fprintf(fp, b"%s\t%u\t%u\t%u\t%u\t%.10e\t%.10f\t%.10e\n",
                island.chrom.c_str(), island.start, island.end, island.obs_count,
                island.control_count, island.pvalue, island.fc, island.alpha_stat
                )

    cdef void write_bed_line(self, Island island, FILE *fp):
        fprintf(fp, b"%s\t%u\t%u\t%u\n",
                island.chrom.c_str(), island.start, island.end, island.obs_count
                )

    cdef void write_scoreisland_line(self, Island island, FILE *fp):
        fprintf(fp, b"%s\t%u\t%u\t%.10f\n",
                island.chrom.c_str(), island.start, island.end, island.score
                )

    cdef void c_write(self, cstr outfile_path):
        cdef FILE *fp = fopen(outfile_path, "w")

        cdef vector[string] chroms = self.islands.getChromosomes()
        cdef cstr line
        cdef vector[Island]* vptr

        cdef format_f write
        if self.file_type == "summary":
            write = self.write_summary_line
        elif self.file_type == "fdr-filtered":
            write = self.write_bed_line
        elif self.file_type == "scoreisland" or self.file_type == "cgisland":
            write = self.write_scoreisland_line

        for i in range(chroms.size()):
            vptr = self.islands.getVectorPtr(chroms[i])
            for j in range(deref(vptr).size()):
                write(self, deref(vptr)[j], fp)

    cpdef void write(self):
        self.file_name += "-W" + str(self.window_size)
        if self.gap_size is not None:
            self.file_name += "-G" + str(self.gap_size)

        if self.file_type == "scoreisland":
            self.file_name += ".scoreisland"
        elif self.file_type == "summary":
            self.file_name += "-islands-summary"
        elif self.file_type == "fdr-filtered":
            self.file_name += "-FDR" + str(self.fdr) + "-island.bed"
        elif self.file_type == "cgisland":
            self.file_name += ".cgisland"

        cdef bytes outfile_path = (self.output_dir + "/" + self.file_name).encode("UTF-8")

        self.c_write(outfile_path)


cdef class BEDFileWriter:
    def __cinit__(self, 
        str file_name, 
        str output_dir, 
        BEDReadContainer reads,
        int window_size,
        object fdr,
        object gap_size = None
    ):
        self.file_name = file_name
        self.output_dir = output_dir
        self.reads = reads
        self.window_size = window_size
        self.fdr = fdr
        self.gap_size = gap_size

    cdef void c_write(self, cstr outfile_path):
        cdef FILE *fp = fopen(outfile_path, "w")

        cdef vector[string] chroms = self.reads.getChromosomes()
        cdef vector[BEDRead]* vptr
        cdef BEDRead read

        for i in range(chroms.size()):
            vptr = self.reads.getVectorPtr(chroms[i])
            for j in range(deref(vptr).size()):
                read = deref(vptr)[j]
                fprintf(fp, "%s\t%u\t%u\t%s\t%d\t%c\n", 
                    read.chrom.c_str(), read.start, read.end, 
                    read.name.c_str(), read.score, read.strand
                )

    cpdef void write(self):
        self.file_name += "-W" + str(self.window_size)
        if self.gap_size is not None:
            self.file_name += "-G" + str(self.gap_size)
        self.file_name += "-FDR" + str(self.fdr) + "-islandfiltered.bed"

        cdef bytes outfile_path = (self.output_dir + "/" + self.file_name).encode("UTF-8")

        self.c_write(outfile_path)


cdef class DiffExprIslandFileWriter:
    def __cinit__(self, 
        str file_name_1, 
        str file_name_2,
        str output_dir,
        str file_type,
        DiffExprIslandContainer islands,
        int window_size,
        object e_value = None,
        object fdr = None,
        object gap_size = None
    ):
        self.file_name_1 = file_name_1
        self.file_name_2 = file_name_2
        self.file_type = file_type
        self.output_dir = output_dir
        self.islands = islands
        self.window_size = window_size
        self.e_value = e_value
        self.fdr = fdr
        self.gap_size = gap_size
        self.header = b"#chrom\tstart\tend\tReadcount_A\tNormalized_Readcount_A\tReadcountB\tNormalized_Readcount_B\tFc_A_vs_B\tpvalue_A_vs_B\tFDR_A_vs_B\tFc_B_vs_A\tpvalue_B_vs_A\tFDR_B_vs_A\n"
        self.format = b"%s\t%u\t%u\t%u\t%.10f\t%u\t%.10f\t%.10f\t%.10e\t%.10e\t%.10f\t%.10e\t%.10e\n"

        if "fdr-filtered" in self.file_type and fdr is None:
            raise ValueError("Missing FDR value")

    cdef void c_write_all(self, cstr outfile_path):
        cdef FILE *fp = fopen(outfile_path, "w")

        if "summary" in self.file_type:
            fprintf(fp, self.header)

        cdef vector[string] chroms = self.islands.getChromosomes()
        cdef vector[DiffExprIsland]* vptr
        cdef DiffExprIsland island

        for i in range(chroms.size()):
            vptr = self.islands.getVectorPtr(chroms[i])
            for j in range(deref(vptr).size()):
                island = deref(vptr)[j]
                fprintf(fp, self.format, 
                    island.chrom.c_str(), island.start, island.end, 
                    island.count_A, island.norm_count_A,
                    island.count_B, island.norm_count_B,
                    island.fc_A_vs_B, island.pvalue_A_vs_B, island.fdr_A_vs_B,
                    island.fc_B_vs_A, island.pvalue_B_vs_A, island.fdr_B_vs_A
                )

    cdef void c_write_basic(self, cstr outfile_path):
        cdef FILE *fp = fopen(outfile_path, "w")

        cdef vector[string] chroms = self.islands.getChromosomes()
        cdef vector[DiffExprIsland]* vptr
        cdef DiffExprIsland island
        
        for i in range(chroms.size()):
            vptr = self.islands.getVectorPtr(chroms[i])
            for j in range(deref(vptr).size()):
                island = deref(vptr)[j]
                fprintf(fp, "%s\t%u\t%u\n", 
                    island.chrom.c_str(), island.start, island.end
                )

    cpdef void write(self):
        if self.file_type == "union-island":
            file_name = self.file_name_1 + "-vs-" + self.file_name_2
        elif self.file_type == "summary":
            file_name = self.file_name_1 + "-and-" + self.file_name_2
        else:
            file_name = self.file_name_1

        file_name += "-W" + str(self.window_size)

        if self.gap_size is not None:
            file_name += "-G" + str(self.gap_size)

        if self.file_type == "union-island":
            if self.e_value is not None:
                file_name += "-E" + str(self.e_value)
            file_name += "-union.island"

        if self.file_type == "fdr-filtered-increased":
            file_name += '-increased-islands-summary-FDR' + str(self.fdr)
        elif self.file_type == "fdr-filtered-decreased":
            file_name += '-decreased-islands-summary-FDR' + str(self.fdr)
        elif self.file_type == "summary":
            file_name += "-summary"

        cdef bytes outfile_path = (self.output_dir + "/" + file_name).encode("UTF-8")

        if self.file_type == "union-island":
            self.c_write_basic(outfile_path)
        else:
            self.c_write_all(outfile_path)

