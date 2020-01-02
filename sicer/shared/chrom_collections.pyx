from sicer.shared import GenomeData
from cython.operator cimport dereference as deref
from cython.operator cimport postincrement
from libcpp.utility cimport pair
import sys

cdef class ChromBEDReadCollection:

    def __cinit__(self, str species):
        self.species = species
        self.chromosomes = list(map(lambda x: x.encode("UTF-8"), GenomeData.species_chroms[species]))

        for chrom in self.chromosomes:
            self.data.insert(pair[string, vec_ptr](chrom, new vector[BEDRead](0)))

    cdef void insertRead(self, string chrom, BEDRead input_data):
        deref(self.data.at(chrom)).push_back(input_data)

    cdef mapcpp[string, vec_ptr] getData(self):
        return self.data

    cdef vec_ptr getChromData(self, string chrom):
        return self.data.at(chrom)

    cdef BEDRead getRead(self, string chrom, int i):
        return deref(self.data.at(chrom)).at(i)

    cpdef list getChromosomes(self):
        return self.chromosomes

    cpdef string dataToString(self):
        cdef string s = "".encode("UTF-8")
        for chrom in self.chromosomes:
            for i in range(deref(self.data.at(chrom)).size()):
                s = s + deref(self.data.at(chrom)).at(i).toString().decode("UTF-8") + "\n"
        return s

    cpdef void printDataHead(self):
        head = 10;
        cdef mapcpp[string, vec_ptr].iterator it = self.data.begin()
        while(it != self.data.end()):
            if deref(deref(it).second).size() > 0:
                for i in range(head):
                    print(deref(deref(it).second).at(i).toString().decode("UTF-8"))
            postincrement(it)

    def __dealloc__(self):
        for chrom in self.chromosomes:
            self.data[chrom].clear()
            self.data[chrom].shrink_to_fit()
