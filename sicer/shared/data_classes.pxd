from libcpp.string cimport string
from libc.stdint cimport uint32_t

ctypedef char* cstr

cdef extern from "data_objects.h" nogil:
    # Represents a single line of BED file
    cdef cppclass BEDRead:
        string chrom
        uint32_t start
        uint32_t end
        string name
        int score
        char strand

        BEDRead()
        BEDRead(string chrom, int start, int end, string name, int score, char strand)
        string toString()

    cdef cppclass Window:
        string chrom
        uint32_t start
        uint32_t end
        uint32_t count

        Window()
        Window(string chrom, uint32_t start, uint32_t end, uint32_t count)
        string toString()

    cdef cppclass Island:
        string chrom
        uint32_t start
        uint32_t end
        double score
        uint32_t obs_count;
        uint32_t control_count;
        double pvalue;
        double fc;
        double alpha_stat;

        Island()
        Island(string chrom, uint32_t start, uint32_t end, double score)
        string toString()

    cdef cppclass DiffExprIsland:
        string chrom
        uint32_t start
        uint32_t end
        uint32_t count_A
        double norm_count_A
        uint32_t count_B
        double norm_count_B
        double fc_A_vs_B
        double pvalue_A_vs_B
        double fdr_A_vs_B
        double fc_B_vs_A
        double pvalue_B_vs_A
        double fdr_B_vs_A

        DiffExprIsland()
        DiffExprIsland(string chrom, uint32_t start, uint32_t end)
        DiffExprIsland(Island island)
        string toString()