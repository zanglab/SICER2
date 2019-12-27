cimport read

cdef class Read:

    def __cinit__(self, chrom, start, end, name, score, strand):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.name = name
        self.score = score
        self.strand = strand