from numpy import log, arctan, pi, zeros, ndarray
from numpy.linalg import solve
cimport numpy

cpdef double division(double a, double b):
    if -1e-8 < b < 1e-8:
        sign = 1 if b > 0 else -1
        return sign * a / 1e-8
    else:
        return a / b

# cpdef double[:] source_velocity(double sigma, double[:] p, double[:] p1, double[:] p2, double[:] p3, double[:] x, double[:] y, double[:] z):

#     # Params
#     cdef double r1 = ((p[0] - p1[0]) ** 2 + (p[1] - p1[1]) ** 2 + p[2] ** 2) ** 0.5
#     cdef double r2 = ((p[0] - p2[0]) ** 2 + (p[1] - p2[1]) ** 2 + p[2] ** 2) ** 0.5
#     cdef double r3 = ((p[0] - p3[0]) ** 2 + (p[1] - p3[1]) ** 2 + p[2] ** 2) ** 0.5
#     cdef double e1 = (p[0] - p1[0]) ** 2 + p[2] ** 2
#     cdef double e2 = (p[0] - p2[0]) ** 2 + p[2] ** 2
#     cdef double e3 = (p[0] - p3[0]) ** 2 + p[2] ** 2
#     cdef double h1 = (p[0] - p1[0]) * (p[1] - p1[1])
#     cdef double h2 = (p[0] - p2[0]) * (p[1] - p2[1])
#     cdef double h3 = (p[0] - p3[0]) * (p[1] - p3[1])

#     # Point 1 and 2
#     cdef double d12 = ((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2) ** 0.5
#     cdef double m12 = division((p2[1] - p1[1]), (p2[0] - p1[0]))

#     # Point 2 and 3
#     cdef double d23 = ((p3[0] - p2[0]) ** 2 + (p3[1] - p2[1]) ** 2) ** 0.5
#     cdef double m23 = division((p3[1] - p2[1]), (p3[0] - p2[0]))

#     # Point 3 and 1
#     cdef double d31 = ((p1[0] - p3[0]) ** 2 + (p1[1] - p3[1]) ** 2) ** 0.5
#     cdef double m31 = division((p1[1] - p3[1]), (p1[0] - p3[0]))

#     # Velocities
#     cdef double ln12 = log( division((r1 + r2 - d12), (r1 + r2 + d12)) )
#     cdef double ln23 = log( division((r2 + r3 - d23), (r2 + r3 + d23)) )
#     cdef double ln31 = log( division((r3 + r1 - d31), (r3 + r1 + d31)) )

#     cdef double mult = (sigma / (4 * pi))

#     cdef double u = -mult * ( division((p2[1] - p1[1]), d12) * ln12 + division((p3[1] - p2[1]), d23) * ln23 + division((p1[1] - p3[1]), d31) * ln31 )
#     cdef double v = mult * ( division((p2[0] - p1[0]), d12) * ln12 + division((p3[0] - p2[0]), d23) * ln23 + division((p1[0] - p3[0]), d31) * ln31 )
#     cdef double w = -mult * ( arctan(division(m12 * e1 - h1, p[2] * r1)) - arctan(division(m12 * e2 - h2, p[2] * r2)) + arctan(division(m23 * e2 - h2, p[2] * r2)) - arctan(division(m23 * e3 - h3, p[2] * r3)) + arctan(division(m31 * e3 - h3, p[2] * r3)) - arctan(division(m31 * e1 - h1, p[2] * r1)) )
    
#     cdef double[:] vel = zeros(3)

#     vel[0] = x[0] * u + y[0] * v + z[0] * w
#     vel[1] = x[1] * u + y[1] * v + z[1] * w
#     vel[2] = x[2] * u + y[2] * v + z[2] * w

#     return vel

cpdef ndarray[float32, ndim=2000] source_distribution(ndarray[float32, ndim=2000] U,
                                    ndarray[float32, ndim=2000] transpiration,
                                    ndarray[float32, ndim=2000] controlPoints,
                                    ndarray[float32, ndim=2000] p1Local,
                                    ndarray[float32, ndim=2000] p2Local,
                                    ndarray[float32, ndim=2000] p3Local,
                                    ndarray[float32, ndim=2000] e1Local,
                                    ndarray[float32, ndim=2000] e2Local,
                                    ndarray[float32, ndim=2000] e3Local):

    # Parameters
    cdef int n
    cdef ndarray p
    cdef ndarray p1
    cdef ndarray p2
    cdef ndarray p3
    cdef ndarray e1
    cdef ndarray e2
    cdef ndarray e3
    cdef ndarray vel
    cdef ndarray a
    cdef ndarray b
    cdef int i
    cdef int j
    cdef double r1
    cdef double r2
    cdef double r3
    cdef double l1
    cdef double l2
    cdef double l3
    cdef double h1
    cdef double h2
    cdef double h3
    cdef double d12
    cdef double m12
    cdef double d23
    cdef double m23
    cdef double d31
    cdef double m31
    cdef double ln12
    cdef double ln23
    cdef double ln31
    cdef double u
    cdef double v
    cdef double w

    cdef double mult

    n = len(transpiration)
    a = zeros((n, n))
    b = zeros(transpiration)
    vel = zeros(3)

    for i in range(n):

        p1 = p1Local[i, :]
        p2 = p2Local[i, :]
        p3 = p3Local[i, :]
        e1 = e1Local[i, :]
        e2 = e2Local[i, :]
        e3 = e3Local[i, :]
        b[i] = transpiration[i] - (U[0] * e3[0] + U[1] * e3[1] + U[2] * e3[2])

        for j in range(n):

            p = controlPoints[i, j, :]

            # Source
            r1 = ((p[0] - p1[0]) ** 2 + (p[1] - p1[1]) ** 2 + p[2] ** 2) ** 0.5
            r2 = ((p[0] - p2[0]) ** 2 + (p[1] - p2[1]) ** 2 + p[2] ** 2) ** 0.5
            r3 = ((p[0] - p3[0]) ** 2 + (p[1] - p3[1]) ** 2 + p[2] ** 2) ** 0.5
            l1 = (p[0] - p1[0]) ** 2 + p[2] ** 2
            l2 = (p[0] - p2[0]) ** 2 + p[2] ** 2
            l3 = (p[0] - p3[0]) ** 2 + p[2] ** 2
            h1 = (p[0] - p1[0]) * (p[1] - p1[1])
            h2 = (p[0] - p2[0]) * (p[1] - p2[1])
            h3 = (p[0] - p3[0]) * (p[1] - p3[1])

            # Point 1 and 2
            d12 = ((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2) ** 0.5
            m12 = division((p2[1] - p1[1]), (p2[0] - p1[0]))

            # Point 2 and 3
            d23 = ((p3[0] - p2[0]) ** 2 + (p3[1] - p2[1]) ** 2) ** 0.5
            m23 = division((p3[1] - p2[1]), (p3[0] - p2[0]))

            # Point 3 and 1
            d31 = ((p1[0] - p3[0]) ** 2 + (p1[1] - p3[1]) ** 2) ** 0.5
            m31 = division((p1[1] - p3[1]), (p1[0] - p3[0]))

            # Velocities
            ln12 = log( division((r1 + r2 - d12), (r1 + r2 + d12)) )
            ln23 = log( division((r2 + r3 - d23), (r2 + r3 + d23)) )
            ln31 = log( division((r3 + r1 - d31), (r3 + r1 + d31)) )

            mult = (1 / (4 * pi))

            u = -mult * ( division((p2[1] - p1[1]), d12) * ln12 + division((p3[1] - p2[1]), d23) * ln23 + division((p1[1] - p3[1]), d31) * ln31 )
            v = mult * ( division((p2[0] - p1[0]), d12) * ln12 + division((p3[0] - p2[0]), d23) * ln23 + division((p1[0] - p3[0]), d31) * ln31 )
            w = -mult * ( arctan(division(m12 * l1 - h1, p[2] * r1)) - arctan(division(m12 * l2 - h2, p[2] * r2)) + arctan(division(m23 * l2 - h2, p[2] * r2)) - arctan(division(m23 * l3 - h3, p[2] * r3)) + arctan(division(m31 * l3 - h3, p[2] * r3)) - arctan(division(m31 * l1 - h1, p[2] * r1)) )
            
            vel[0] = u * e1[0] + v * e2[0] + w * e3[0]
            vel[1] = u * e1[1] + v * e2[1] + w * e3[1]
            vel[2] = u * e1[2] + v * e2[2] + w * e3[2]

            a[j, i] = (vel[0] + U[0]) * e3[0] + (vel[1] + U[1]) * e3[1] + (vel[2] + U[2]) * e3[2]
    
    sol = solve(a, b)

    return sol