# Python Imports
import os, errno

# SICER Internal Imports
from sicer.shared import GenomeData
from sicer.shared.data_classes cimport BEDRead
from sicer.shared.chrom_collections cimport ChromBEDReadCollection
from sicer.utils_cpp cimport to_string, remove_at

# Cython Imports
cimport cython
from libcpp cimport bool
from cython.operator cimport dereference as deref
from cython.operator cimport postincrement as postinc
from libc.stdlib cimport free, atoi, malloc
from libc.string cimport strtok
from libc.stdio cimport fopen, fclose, getline, printf
from libcpp.map cimport map as mapcpp
from libcpp.set cimport set as setcpp
from libcpp.vector cimport vector
from libcpp.string cimport string
from cython.parallel import parallel, prange
from libcpp.algorithm cimport sort

ctypedef char* cstr
ctypedef bool (*cmp_f)(BEDRead, BEDRead)
ctypedef vector[BEDRead]* vec_ptr

cdef char PLUS = b'+'
cdef string TAB = string(b'\t')
cdef string NEWLINE = string(b'\n')

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
        str species
        set chromosomes
        int num_cpu
        int redundancy_threshold
        int line_count

    def __cinit__(
        self, 
        str file_name, 
        str species,
        int num_cpu,
        int redundancy_threshold
    ):
        self.file_name = file_name
        self.species = species
        self.chromosomes = set(map(lambda x: x.encode("UTF-8"), GenomeData.species_chroms[species]))
        self.num_cpu = num_cpu
        self.redundancy_threshold = redundancy_threshold
        self.line_count = 0

    cdef (int, int, int, int) _remove_redudant_reads(self, vec_ptr vec, int threshold) nogil:
        # Indices to delete
        cdef vector[int] shouldDelete

        cdef int start = -1
        cdef int end = -1
        cdef int redund_count = 0
        cdef bool retain = 0

        cdef int plus_total_count = 0
        cdef int plus_retained_count = 0
        cdef int minus_total_count = 0
        cdef int minus_retained_count = 0

        for i in range(deref(vec).size()):
            read = deref(vec).at(i)

            if read.start != start or read.end != end:
                start = read.start
                end = read.end
                redund_count = 1
                retain = 1
            else:
                postinc(redund_count)
                if redund_count <= threshold:
                    retain = 1
                else:
                    # Delete read
                    retain = 0
                    shouldDelete.push_back(i)

            if read.strand == PLUS:
                postinc(plus_total_count)
                if retain:
                    postinc(plus_retained_count)
            else:
                postinc(minus_total_count)
                if retain:
                    postinc(minus_retained_count)

        # Remove duplicate reads
        deref(vec).erase(remove_at(deref(vec).begin(), deref(vec).end(), shouldDelete.begin(), shouldDelete.end()), deref(vec).end())

        return plus_total_count, plus_retained_count, minus_total_count, minus_retained_count

    cdef string _preprocess_BED_reads_by_chrom(self, vec_ptr vec) nogil:
        # First, sort reads by strand, start pos, and end pos order
        sort[vector[BEDRead].iterator, cmp_f](deref(vec).begin(), deref(vec).end(), compare_reads)

        cdef int pt
        cdef int pr
        cdef int mt
        cdef int mr

        pt, pr, mt, mr = self._remove_redudant_reads(vec, self.redundancy_threshold)

        return to_string(pt) + TAB + to_string(pr) + TAB + to_string(mt) + TAB + to_string(mr)

    cdef ChromBEDReadCollection _preprocess_BED_reads(self, ChromBEDReadCollection reads):
        print("Preprocessing reads...")

        # Convert Python list to C-array for no-GIL use
        cdef int num_chroms = len(reads.getChromosomes())
        cdef cstr *chroms = <cstr *>malloc(num_chroms*cython.sizeof(cstr))
        if chroms is NULL:
            raise MemoryError()
        for j in range(num_chroms):
            chroms[j] = reads.getChromosomes()[j]

        cdef mapcpp[string, vec_ptr] data = reads.getData()
        cdef int i
        cdef string s

        for i in prange(num_chroms, schedule='guided', num_threads=self.num_cpu, nogil=True):
            s += self._preprocess_BED_reads_by_chrom(data[chroms[i]]) + NEWLINE

        free(chroms)

        # Print summary of preprocessing
        print("Chromsome\t Total Plus Reads \t Retained Plus Reads \t Total Minus Reads \t Retained Minus Reads")
        print(s.decode("UTF-8"))

        return reads

    cdef BEDRead * _parseBEDLine(self, cstr line):
        cdef cstr read[6]
        cdef int count = 0

        cdef cstr token = strtok(line, "\t")

        while token != NULL:
            read[count] = token
            count += 1
            token = strtok(NULL, "\t")

        if count != 6:
            raise ValueError("Not a valid BED6 line: %s " % line.decode("UTF-8"))

        return new BEDRead(string(read[0]), atoi(read[1]), atoi(read[2]), string(read[3]), atoi(read[4]), read[5][0])

    cdef ChromBEDReadCollection _read_file(self):
        print("Reading file \"" + self.file_name + "\" ...")
        fp = fopen(self.file_name.encode("UTF-8"), "r")
 
        if fp == NULL:
            raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), self.file_name)

        cdef ChromBEDReadCollection read_collection = ChromBEDReadCollection(self.species)

        cdef cstr line = NULL
        cdef size_t len = 0
        cdef BEDRead read

        while getline(&line, &len, fp) != -1:
            self.line_count+=1
            read = deref(self._parseBEDLine(line))
            if read.chrom.c_str() in self.chromosomes: 
                read_collection.insertRead(read.chrom, read)

        fclose(fp)
        free(line)

        return read_collection
   
    cpdef ChromBEDReadCollection read_file(self):
        return self._preprocess_BED_reads(self._read_file())