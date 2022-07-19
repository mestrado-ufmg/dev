import numpy
cimport numpy

cpdef double[:] func(double[:] a, double[:] b):

    cdef int n
    cdef int i
    cdef double[:] out

    n = len(a)
    out = numpy.zeros(n)

    for i in range(n):
        out[i] = a[i] / b[i]

    return out