# Python Imports
import os, errno

# SICER Internal Imports
from sicer.shared.data_classes cimport BEDRead
from sicer.shared.chrom_containers cimport ChromBEDReadContainer
from sicer.utility.utils cimport to_string, remove_at

# Cython Imports
from libcpp cimport bool
from cython.operator cimport dereference as deref
from cython.operator cimport preincrement as preinc
from libc.stdlib cimport free, atoi, malloc
from libc.string cimport strtok
from libc.stdio cimport fopen, fclose, getline, printf
from libcpp.map cimport map as mapcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libcpp.algorithm cimport sort

# Typedefs
ctypedef char* cstr
ctypedef bool (*cmp_f)(BEDRead, BEDRead)
ctypedef vector[BEDRead]* vec_ptr

cdef char PLUS = b'+'

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
        int line_count

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
        self.line_count = 0

    cdef void _remove_redudant_reads(self, vec_ptr vec, int threshold) nogil:
        # Indices to delete
        cdef vector[int] shouldDelete

        cdef int start = -1
        cdef int end = -1
        cdef int redund_count = 0

        for i in range(deref(vec).size()):
            read = deref(vec)[i]

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
        deref(vec).erase(remove_at(deref(vec).begin(), deref(vec).end(), shouldDelete.begin(), shouldDelete.end()), deref(vec).end())

    cdef void _preprocess_BED_reads_by_chrom(self, vec_ptr vec) nogil:
        # First, sort reads by strand, start pos, and end pos order
        sort[vector[BEDRead].iterator, cmp_f](deref(vec).begin(), deref(vec).end(), compare_reads)

        self._remove_redudant_reads(vec, self.redundancy_threshold)

    cdef ChromBEDReadContainer _preprocess_BED_reads(self, ChromBEDReadContainer reads):
        print("Preprocessing reads...")

        # Convert Python list to vector for no-GIL use
        cdef int num_chroms = len(reads.getChromosomes())
        cdef vector[string] chroms = reads.getChromosomes()

        cdef mapcpp[string, vec_ptr] data = reads.getData()
        cdef int i

        for i in prange(num_chroms, schedule='guided', num_threads=self.num_cpu, nogil=True):
            self._preprocess_BED_reads_by_chrom(data.at(chroms.at(i)))

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
            raise ValueError("Not a valid BED6 line: %s " % line.decode("UTF-8"))

        return BEDRead(string(read[0]), atoi(read[1]), atoi(read[2]), string(read[3]), atoi(read[4]), read[5][0])

    cdef ChromBEDReadContainer _read_file(self):
        print("Reading file \"" + self.file_name + "\" ...")
        fp = fopen(self.file_name.encode("UTF-8"), "r")
 
        if fp == NULL:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.file_name)

        cdef ChromBEDReadContainer reads = ChromBEDReadContainer(self.genome_data)

        cdef set chromosomes = set(map(lambda x: x.encode("UTF-8"), self.genome_data.chrom))
        cdef cstr line = NULL
        cdef size_t len = 0
        cdef BEDRead read

        while getline(&line, &len, fp) != -1:
            self.line_count+=1
            read = self._parseBEDLine(line)
            if read.chrom.c_str() in chromosomes: 
                reads.insertRead(read.chrom, read)

        fclose(fp)
        free(line)

        return reads
   
    cpdef ChromBEDReadContainer read_file(self):
        return self._preprocess_BED_reads(self._read_file())