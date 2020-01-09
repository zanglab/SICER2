from libcpp.string cimport string

ctypedef char* cstr

cdef extern from "data_objects.h" nogil:
    # Represents a single line of BED file
    cdef cppclass BEDRead:
        string chrom
        int start
        int end
        string name
        int score
        char strand

        BEDRead()
        BEDRead(string chrom, int start, int end, string name, int score, char strand)
        string toString()

    cdef cppclass Window:
        string chrom
        int start
        int end
        int count

        Window()
        Window(string chrom, int start, int end, int count)
        string toString()

    cdef cppclass Island:
        string chrom
        int start
        int end
        double score
        int obs_count;
        int control_count;
        double pvalue;
        double fc;
        double alpha_stat;

        Island()
        Island(string chrom, int start, int end, double score)
        string toString()