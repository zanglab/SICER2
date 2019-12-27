# Represents a single line of BED file

cdef class Read:
    cdef public:
        char *chrom
        int start
        int end
        char * name
        int score
        char strand