from sicer.utility.utils cimport to_string

from cython.operator cimport dereference as deref
from cython.operator cimport postincrement
from libc.stdlib cimport free
from libc.stdio cimport printf
from libcpp.utility cimport pair

cdef class ChromBEDReadContainer:

    def __cinit__(self, object genome_data):
        self.species = genome_data.species
        self.chromosomes = list(map(lambda x: x.encode("UTF-8"), genome_data.chrom))
        self.read_count = 0
        for chrom in self.chromosomes:
            self.data.insert(pair[string, bed_vec_ptr](chrom, new vector[BEDRead](0)))

    cdef void insertRead(self, string chrom, BEDRead item):
        deref(self.data.at(chrom)).push_back(item)

    cpdef void updateReadCount(self):
        cdef int new_count = 0
        for chrom in self.chromosomes:
            new_count = new_count + deref(self.data.at(chrom)).size()

        self.read_count = new_count

    cdef mapcpp[string, bed_vec_ptr] getData(self):
        return self.data

    cdef bed_vec_ptr getChromVector(self, string chrom) nogil:
        return self.data.at(chrom)

    cdef BEDRead getRead(self, string chrom, int i):
        return deref(self.data.at(chrom)).at(i)

    cpdef list getChromosomes(self):
        return self.chromosomes

    cpdef int getReadCount(self):
        return self.read_count

    cpdef void printDataHead(self):
        head = 10;
        cdef mapcpp[string, bed_vec_ptr].iterator it = self.data.begin()
        while(it != self.data.end()):
            if deref(deref(it).second).size() > 0:
                for i in range(head):
                    print(deref(deref(it).second).at(i).toString().decode("UTF-8"))
            postincrement(it)

    def __dealloc__(self):
        for chrom in self.chromosomes:
            deref(self.data.at(chrom)).clear()
            deref(self.data.at(chrom)).shrink_to_fit()
            free(self.data.at(chrom))


cdef class ChromWindowContainer:

    def __cinit__(self, object genome_data):
        self.species = genome_data.species
        self.chromosomes = list(map(lambda x: x.encode("UTF-8"), genome_data.chrom))
        self.window_count = 0
        for chrom in self.chromosomes:
            self.data.insert(pair[string, win_vec_ptr](chrom, new vector[Window](0)))

    cdef void insertWindow(self, string chrom, Window item):
        deref(self.data.at(chrom)).push_back(item)

    cpdef void updateCounts(self):
        # Update two counters
        cdef int new_window_count = 0
        cdef int new_tag_count = 0
        for chrom in self.chromosomes:
            new_window_count += deref(self.data.at(chrom)).size()
            for i in range(deref(self.data[chrom]).size()):
                new_tag_count += deref(self.data[chrom])[i].count

        self.window_count = new_window_count
        self.total_tag_count = new_tag_count

    cdef mapcpp[string, win_vec_ptr] getData(self):
        return self.data

    cdef win_vec_ptr getChromVector(self, string chrom) nogil:
        return self.data.at(chrom)

    cdef Window getWindow(self, string chrom, int index):
        return deref(self.data.at(chrom)).at(index)

    cpdef list getChromosomes(self):
        return self.chromosomes

    cpdef int getWindowCount(self):
        return self.window_count

    cpdef int getTotalTagCount(self):
        return self.total_tag_count

    cpdef void printDataHead(self):
        head = 10;
        cdef mapcpp[string, win_vec_ptr].iterator it = self.data.begin()
        while(it != self.data.end()):
            if deref(deref(it).second).size() > 0:
                for i in range(head):
                    print(deref(deref(it).second).at(i).toString().decode("UTF-8"))
            postincrement(it)

    cpdef void printChrCount(self):
        cdef mapcpp[string, win_vec_ptr].iterator it = self.data.begin()
        cdef string line
        for c in self.chromosomes:
            line = c + b": " +  to_string(deref(self.data[c]).size()) + b"\n"
            printf(line.c_str())
            postincrement(it)

    def __dealloc__(self):
        for chrom in self.chromosomes:
            deref(self.data.at(chrom)).clear()
            deref(self.data.at(chrom)).shrink_to_fit()
            free(self.data.at(chrom))


cdef class ChromIslandContainer:

    def __cinit__(self, object genome_data):
        self.species = genome_data.species
        self.chromosomes = list(map(lambda x: x.encode("UTF-8"), genome_data.chrom))
        self.island_count = 0
        for chrom in self.chromosomes:
            self.data.insert(pair[string, isl_vec_ptr](chrom, new vector[Island](0)))

    cdef void insertIsland(self, string chrom, Island item):
        deref(self.data.at(chrom)).push_back(item)

    cpdef void updateIslandCount(self):
        cdef int new_count = 0
        for chrom in self.chromosomes:
            new_count += deref(self.data.at(chrom)).size()

        self.island_count = new_count

    cdef mapcpp[string, isl_vec_ptr] getData(self):
        return self.data

    cdef isl_vec_ptr getChromVector(self, string chrom) nogil:
        return self.data.at(chrom)

    cdef Island getIsland(self, string chrom, int index):
        return deref(self.data.at(chrom)).at(index)

    cpdef list getChromosomes(self):
        return self.chromosomes

    cpdef int getIslandCount(self):
        return self.island_count

    cpdef void printDataHead(self):
        head = 10;
        cdef mapcpp[string, isl_vec_ptr].iterator it = self.data.begin()
        while(it != self.data.end()):
            if deref(deref(it).second).size() > 0:
                for i in range(head):
                    print(deref(deref(it).second).at(i).toString().decode("UTF-8"))
            postincrement(it)

    cpdef void printChrCount(self):
        cdef mapcpp[string, isl_vec_ptr].iterator it = self.data.begin()
        cdef string line
        for c in self.chromosomes:
            line = c + b": " +  to_string(deref(self.data[c]).size()) + b"\n"
            printf(line.c_str())
            postincrement(it)

    def __dealloc__(self):
        for chrom in self.chromosomes:
            deref(self.data.at(chrom)).clear()
            deref(self.data.at(chrom)).shrink_to_fit()
            free(self.data.at(chrom))