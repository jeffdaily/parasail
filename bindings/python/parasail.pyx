cdef extern from "parasail.h":
    ctypedef struct parasail_result_t:
        pass
    ctypedef struct parasail_matrix_t:
        pass
    parasail_result_t* nw(const char * s1, const int s1Len, const char * s2, const int s2Len, const int open, const int gap, const parasail_matrix_t* matrix)

