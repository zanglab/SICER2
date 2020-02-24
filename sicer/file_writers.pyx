# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead, Window, DiffExprIsland
from sicer.shared.utils cimport to_string

#Cython Imports
from libc.stdio cimport FILE, fopen, snprintf, fprintf, fclose, printf
from cython.operator cimport preincrement as preinc
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.operator cimport dereference as deref

ctypedef vector[Window]* win_vec_ptr
ctypedef cstr (*format_f)(IslandFileWriter, Island)

cdef class WigFileWriter:
    # Writes .wig file of windows

    def __cinit__(self, 
        str file_name, 
        str output_dir, 
        WindowContainer windows,
        int window_size,
        bint filtered,
        object fdr = None
    ):
        self.file_name = file_name
        self.output_dir = output_dir
        self.windows = windows
        self.window_size = window_size
        self.filtered = filtered

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
                    fprintf(fp, "%d\t%.2f\n", deref(vptr)[j].start + 1, deref(vptr)[j].count / scaling_factor)

    cpdef void write(self):
        print("Normalizing graphs by total island filitered reads per million and generating summary WIG file...\n")
        # We first need to normalize
        cdef int count = self.windows.getTotalTagCount()
        cdef double scaling_factor =  count / 1000000.0 * (self.window_size / 1000.0)

        # Format final file_name
        self.file_name += "-W" + str(self.window_size)
        if self.filtered:
            self.file_name += "-islandfiltered" + "-FDR" + str(self.fdr)

        cdef bytes wig_header = ("track type=wiggle_0 name=" + self.file_name + "\n").encode("UTF-8")
        self.file_name += "-normalized.wig"

        cdef bytes outfile_path = (self.output_dir + "/" + self.file_name).encode("UTF-8")

        self.c_write(outfile_path, wig_header, self.window_size, scaling_factor)


cdef class IslandFileWriter:
    # Writes island files

    def __cinit__(self, 
        str file_name, 
        str output_dir,
        object file_type,
        IslandContainer islands,
        int window_size,
        object gap_size = None,
        object fdr = None
    ):
        self.file_name = file_name
        self.output_dir = output_dir
        self.file_type = file_type
        self.islands = islands
        self.window_size = window_size
        self.gap_size = gap_size
        self.fdr = fdr

        if file_type == "fdr-filtered" and fdr is None:
            raise ValueError("Missing FDR value")

    cdef cstr format_summary_line(self, Island island):
        cdef char buffer[128]
        snprintf(buffer, 128, "%s\t%d\t%d\t%d\t%d\t%.10e\t%.10f\t%.10e\n",
                island.chrom.c_str(), island.start, island.end, island.obs_count,
                island.control_count, island.pvalue, island.fc, island.alpha_stat
                )
        return buffer

    cdef cstr format_bed_line(self, Island island):
        cdef char buffer[64]
        snprintf(buffer, 64, "%s\t%d\t%d\t%d\n",
                island.chrom.c_str(), island.start, island.end, island.obs_count
                )
        return buffer

    cdef cstr format_scoreisland_line(self, Island island):
        cdef char buffer[64]
        snprintf(buffer, 64, "%s\t%d\t%d\t%.10f\n",
                island.chrom.c_str(), island.start, island.end, island.score
                )
        return buffer

    cdef void c_write(self, cstr outfile_path):
        cdef FILE *fp = fopen(outfile_path, "w")

        cdef vector[string] chroms = self.islands.getChromosomes()
        cdef cstr line
        cdef vector[Island]* vptr

        cdef format_f func
        if self.file_type == "summary":
            func = self.format_summary_line
        elif self.file_type == "fdr-filtered":
            func = self.format_bed_line
        elif self.file_type == "scoreisland":
            func = self.format_scoreisland_line
        elif self.file_type == "cgisland":
            func = self.format_bed_line

        for i in range(chroms.size()):
            vptr = self.islands.getVectorPtr(chroms[i])
            for j in range(deref(vptr).size()):
                line = func(self, deref(vptr)[j])
                fprintf(fp, line)

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
                fprintf(fp, "%s\t%d\t%d\t%s\t%d\t%c\n", 
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


cdef class DiffExprIslandWriter:
    def __cinit__(self, 
        str file_name_1, 
        str file_name_2,
        str output_dir, 
        DiffExprIslandContainer islands,
        int window_size,
        bint fdr_filtered,
        bint increased,
        object fdr = None,
        object gap_size = None
    ):
        self.file_name_1 = file_name_1
        self.file_name_2 = file_name_2
        self.output_dir = output_dir
        self.islands = islands
        self.window_size = window_size
        self.fdr_filtered = fdr_filtered
        self.increased = increased
        self.fdr = fdr
        self.gap_size = gap_size
        self.header = b"#chrom\tstart\tend\tReadcount_A\tNormalized_Readcount_A\tReadcountB\tNormalized_Readcount_B\tFc_A_vs_B\tpvalue_A_vs_B\tFDR_A_vs_B\tFc_B_vs_A\tpvalue_B_vs_A\tFDR_B_vs_A"
        self.format = b"%s\t%d\t%d\t%d\t%.10f\t%d\t%.10f\t%.10f\t%.10e\t%.10e\t%.10f\t%.10e\t%.10e\n"

        if fdr_filtered and fdr is None:
            raise ValueError("Missing FDR value")

    cdef void c_write(self, cstr outfile_path):
        cdef FILE *fp = fopen(outfile_path, "w")

        if not self.fdr_filtered:
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

    cpdef void write(self):
        file_name = self.file_name_1 + "-and-" + self.file_name_2
        file_name += "-W" + str(self.window_size)

        if self.gap_size is not None:
            file_name += "-G" + str(self.gap_size)

        if self.fdr_filtered:
            if self.increased:
                file_name += '-increased-islands-summary-FDR' + str(self.fdr)
            else:
                file_name += '-decreased-islands-summary-FDR' + str(self.fdr)
        else:
            file_name += "-summary"

        cdef bytes outfile_path = (self.output_dir + "/" + file_name).encode("UTF-8")

        self.c_write(outfile_path)