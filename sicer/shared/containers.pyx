from sicer.shared.utils cimport to_string

from cython.operator cimport dereference as deref
from cython.operator cimport postincrement
from libc.stdlib cimport free
from libc.stdio cimport printf
from libcpp.utility cimport pair

cdef class BEDReadContainer:

    def __cinit__(self, object genome_data):
        self.species = genome_data.species
        self.chromosomes = list(map(lambda x: x.encode("UTF-8"), genome_data.chrom))
        self.read_count = 0
        self.total_count = 0
        for chrom in self.chromosomes:
            self.data.insert(pair[string, vector[BEDRead]](chrom, vector[BEDRead]()))

    cdef void insertRead(self, string chrom, BEDRead item):
        self.data.at(chrom).push_back(item)

    cpdef void updateReadCount(self):
        cdef int new_count = 0
        for chrom in self.chromosomes:
            new_count = new_count + self.data[chrom].size()

        self.read_count = new_count

    cpdef void setTotalReadCount(self):
        cdef int new_count = 0
        for chrom in self.chromosomes:
            new_count = new_count + self.data[chrom].size()

        self.total_count = new_count

    cdef mapcpp[string, vector[BEDRead]] getData(self):
        return self.data

    cdef vector[BEDRead]* getVectorPtr(self, string chrom) nogil:
        return &self.data.at(chrom)

    cpdef list getChromosomes(self):
        return self.chromosomes

    cpdef uint32_t getReadCount(self):
        return self.read_count

    cpdef uint32_t getTotalReadCount(self):
        return self.total_count

    def __dealloc__(self):
        for chrom in self.chromosomes:
            self.data.at(chrom).clear()
            self.data.at(chrom).shrink_to_fit()


cdef class WindowContainer:

    def __cinit__(self, object genome_data):
        self.species = genome_data.species
        self.chromosomes = list(map(lambda x: x.encode("UTF-8"), genome_data.chrom))
        self.window_count = 0
        for chrom in self.chromosomes:
            self.data.insert(pair[string, vector[Window]](chrom, vector[Window]()))

    cpdef void updateCounts(self):
        # Update two counters
        cdef int new_window_count = 0
        cdef int new_tag_count = 0
        for chrom in self.chromosomes:
            new_window_count += self.data.at(chrom).size()
            for i in range(self.data[chrom].size()):
                new_tag_count += self.data[chrom][i].count

        self.window_count = new_window_count
        self.total_tag_count = new_tag_count

    cdef mapcpp[string, vector[Window]] getData(self):
        return self.data

    cdef vector[Window]* getVectorPtr(self, string chrom) nogil:
        return &self.data.at(chrom)

    cpdef list getChromosomes(self):
        return self.chromosomes

    cpdef uint32_t getWindowCount(self):
        return self.window_count

    cpdef uint32_t getTotalTagCount(self):
        return self.total_tag_count

    def __dealloc__(self):
        for chrom in self.chromosomes:
            self.data.at(chrom).clear()
            self.data.at(chrom).shrink_to_fit()


cdef class IslandContainer:

    def __cinit__(self, object genome_data):
        self.species = genome_data.species
        self.chromosomes = list(map(lambda x: x.encode("UTF-8"), genome_data.chrom))
        self.island_count = 0
        for chrom in self.chromosomes:
            self.data.insert(pair[string, vector[Island]](chrom, vector[Island]()))

    cpdef void updateIslandCount(self):
        cdef int new_count = 0
        for chrom in self.chromosomes:
            new_count += self.data.at(chrom).size()

        self.island_count = new_count

    cdef mapcpp[string, vector[Island]] getData(self):
        return self.data

    cdef vector[Island]* getVectorPtr(self, string chrom) nogil:
        return &self.data.at(chrom)

    cpdef list getChromosomes(self):
        return self.chromosomes

    cpdef uint32_t getIslandCount(self):
        return self.island_count

    def __dealloc__(self):
        for chrom in self.chromosomes:
            self.data.at(chrom).clear()
            self.data.at(chrom).shrink_to_fit()

cdef class DiffExprIslandContainer:

    def __cinit__(self, object genome_data):
        self.species = genome_data.species
        self.chromosomes = list(map(lambda x: x.encode("UTF-8"), genome_data.chrom))
        self.island_count = 0
        for chrom in self.chromosomes:
            self.data.insert(pair[string, vector[DiffExprIsland]](chrom, vector[DiffExprIsland]()))

    cpdef void updateIslandCount(self):
        cdef int new_count = 0
        for chrom in self.chromosomes:
            new_count += self.data.at(chrom).size()

        self.island_count = new_count

    cdef mapcpp[string, vector[DiffExprIsland]] getData(self):
        return self.data

    cdef vector[DiffExprIsland]* getVectorPtr(self, string chrom) nogil:
        return &self.data.at(chrom)

    cpdef list getChromosomes(self):
        return self.chromosomes

    cpdef uint32_t getIslandCount(self):
        return self.island_count

    def __dealloc__(self):
        for chrom in self.chromosomes:
            self.data.at(chrom).clear()
            self.data.at(chrom).shrink_to_fit()
