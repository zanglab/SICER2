# Python Imports
import os, errno

# SICER Internal Imports
from sicer.shared import GenomeData
from sicer.shared.data_classes cimport BEDRead
from sicer.shared.chrom_collections cimport ChromBEDReadCollection

# Cython Imports
from cython.operator cimport dereference as deref
from libc.stdlib cimport free, atoi
from libc.string cimport strtok
from libc.stdio cimport fopen, fclose, getline, printf
from libcpp.map cimport map as mapcpp
from libcpp.set cimport set as setcpp
from libcpp.vector cimport vector
from libcpp.string cimport string

ctypedef char* cstr

cdef class BEDReader:
    cdef:
        str file_name
        str species
        set chromosomes
        int line_count

    def __cinit__(self, str file_name, str species):
        self.file_name = file_name
        self.species = species
        self.chromosomes = set(map(lambda x: x.encode("UTF-8"), GenomeData.species_chroms[species]))
        self.line_count = 0

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
        return self._read_file()