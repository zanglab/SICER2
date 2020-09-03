# Python Imports
import os, errno

# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead
from sicer.shared.containers cimport BEDReadContainer
from sicer.shared.utils cimport to_string, remove_at

# Cython Imports
from libcpp cimport bool
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libc.stdlib cimport free, strtoul, malloc
from libc.string cimport strtok
from libc.stdio cimport fopen, fclose, getline, printf
from libc.stdint cimport uint32_t
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libcpp.algorithm cimport sort

# Typedefs
ctypedef char* cstr
ctypedef bool (*cmp_f)(BEDRead, BEDRead)

cdef bool compare_reads(BEDRead i, BEDRead j) nogil:
        if i.strand != j.strand:
            return i.strand < j.strand
        else:
            if i.start != j.start:
                return i.start < j.start
            else:
                return i.end < j.end

cdef class BEDReader:
    cdef:
        str file_name
        object genome_data
        int num_cpu
        int redundancy_threshold

    def __cinit__(
        self, 
        str file_name, 
        object genome_data,
        int num_cpu,
        int redundancy_threshold
    ):
        self.file_name = file_name
        self.genome_data = genome_data
        self.num_cpu = num_cpu
        self.redundancy_threshold = redundancy_threshold

    cdef void _remove_redudant_reads(
        self, 
        vector[BEDRead]& reads, 
        int threshold
    ) nogil:
        # Indices to delete
        cdef vector[uint32_t] shouldDelete

        cdef uint32_t start = -1
        cdef uint32_t end = -1
        cdef int redund_count = 0

        for i in range(reads.size()):
            read = reads[i]

            if read.start != start or read.end != end:
                # Retain read
                start = read.start
                end = read.end
                redund_count = 1
            else:
                preinc(redund_count)
                if redund_count > threshold:
                    # Delete read
                    shouldDelete.push_back(i)

        # Remove duplicate reads
        reads.erase(remove_at(reads.begin(), reads.end(), shouldDelete.begin(), shouldDelete.end()), reads.end())

    cdef void _preprocess_BED_reads_by_chrom(self, vector[BEDRead]& reads) nogil:
        # First, sort reads by strand, start pos, and end pos order
        sort[vector[BEDRead].iterator, cmp_f](reads.begin(), reads.end(), compare_reads)

        self._remove_redudant_reads(reads, self.redundancy_threshold)

    cdef BEDReadContainer _preprocess_BED_reads(self, BEDReadContainer reads):
        print("Preprocessing reads...")

        # Convert Python list to vector for no-GIL use
        cdef vector[string] chroms = reads.getChromosomes()
        cdef vector[BEDRead]* read_vec
        cdef int i
        for i in prange(chroms.size(), schedule='guided', num_threads=self.num_cpu, nogil=True):
            self._preprocess_BED_reads_by_chrom(deref(reads.getVectorPtr(chroms.at(i))))

        reads.updateReadCount()

        print("Retained read count: ", reads.getReadCount())

        return reads

    cdef BEDRead _parseBEDLine(self, cstr line):
        cdef cstr read[6]
        cdef int count = 0

        cdef cstr token = strtok(line, "\t")

        while token != NULL:
            read[count] = token
            preinc(count)
            token = strtok(NULL, "\t")

        if count != 6:
            print("Error: Not a valid BED6 line: %s " % line.decode("UTF-8"))
            raise ValueError("Not a valid BED6 line: %s " % line.decode("UTF-8"))

        return BEDRead(string(read[0]), strtoul(read[1],NULL,10), strtoul(read[2],NULL,10), 
                string(read[3]), strtoul(read[4],NULL,10), read[5][0])

    cdef BEDReadContainer _read_file(self):
        print("Reading file \"" + self.file_name + "\" ...")
        fp = fopen(self.file_name.encode("UTF-8"), "r")
 
        if fp == NULL:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.file_name)

        cdef BEDReadContainer reads = BEDReadContainer(self.genome_data)

        cdef set chromosomes = set(map(lambda x: x.encode("UTF-8"), self.genome_data.chrom))
        cdef cstr line = NULL
        cdef size_t len = 0
        cdef BEDRead read

        while getline(&line, &len, fp) != -1:
            read = self._parseBEDLine(line)
            if read.chrom.c_str() in chromosomes:
                reads.insertRead(read.chrom, read)

        fclose(fp)
        free(line)

        return reads
   
    cpdef BEDReadContainer read_file(self):
        return self._preprocess_BED_reads(self._read_file())
