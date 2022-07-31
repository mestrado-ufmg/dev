import ctypes
from dataclasses import dataclass
from math import pi
import numpy.ctypeslib as npct
from numpy import array, concatenate, double, linspace, meshgrid, ndarray, zeros
import matplotlib.pyplot as plt


# Types
INT = ctypes.c_int
DOUBLE = ctypes.c_double
PDOUBLE = ctypes.POINTER(DOUBLE)
PPDOUBLE = ctypes.POINTER(PDOUBLE)
PPPDOUBLE = ctypes.POINTER(PPDOUBLE)

# Define the structures
class Point3DStruct(ctypes.Structure):
    _fields_ = [
        ('x', DOUBLE),
        ('y', DOUBLE),
        ('z', DOUBLE),
    ]

class Point2DStruct(ctypes.Structure):
    _fields_ = [
        ('x', DOUBLE),
        ('y', DOUBLE),
    ]

@dataclass
class Point3D:
    x: float
    y: float
    z: float

@dataclass
class Point2D:
    x: float
    y: float

def doublePointerToArray(ptr, n):

    arr = zeros(n)

    for i in range(n):
        arr[i] = ptr[i]

    return arr

def source(p: Point3D, p1: Point2D, p2: Point2D, p3: Point2D, e1: Point3D, e2: Point3D, e3: Point3D, area: float, maxDistance: float) -> ndarray:
    
    lib = npct.load_library('lib', './')

    # Define the return type of the C function
    lib.source.restype = ctypes.c_void_p

    # Define arguments of the C function
    lib.source.argtypes = [
        Point3DStruct,
        Point2DStruct,
        Point2DStruct,
        Point2DStruct,
        Point3DStruct,
        Point3DStruct,
        Point3DStruct,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
    ]

    pStruct = Point3DStruct()
    pStruct.x, pStruct.y, pStruct.z = p[0], p[1], p[2]

    p1Struct = Point2DStruct()
    p1Struct.x, p1Struct.y = p1[0], p1[1]

    p2Struct = Point2DStruct()
    p2Struct.x, p2Struct.y = p2[0], p2[1]

    p3Struct = Point2DStruct()
    p3Struct.x, p3Struct.y = p3[0], p3[1]

    e1Struct = Point3DStruct()
    e1Struct.x, e1Struct.y, e1Struct.z = e1[0], e1[1], e1[2]

    e2Struct = Point3DStruct()
    e2Struct.x, e2Struct.y, e2Struct.z = e2[0], e2[1], e2[2]

    e3Struct = Point3DStruct()
    e3Struct.x, e3Struct.y, e3Struct.z = e3[0], e3[1], e3[2]

    out = npct.as_ctypes(zeros(3).astype(double))

    lib.source(pStruct, p1Struct, p2Struct, p3Struct, e1Struct, e2Struct, e3Struct, area, maxDistance, out)
    
    return doublePointerToArray(out, 3)

def doublet(p: Point3D, p1: Point3D, p2: Point3D, p3: Point3D, e1: Point3D, e2: Point3D, e3: Point3D, area: float, maxDistance: float) -> ndarray:
    
    lib = npct.load_library('lib', './')

    # Define the return type of the C function
    lib.doublet.restype = ctypes.c_void_p

    # Define arguments of the C function
    lib.doublet.argtypes = [
        Point3DStruct,
        Point3DStruct,
        Point3DStruct,
        Point3DStruct,
        Point3DStruct,
        Point3DStruct,
        Point3DStruct,
        ctypes.c_double,
        ctypes.c_double,
        ctypes.POINTER(ctypes.c_double),
    ]

    pStruct = Point3DStruct()
    pStruct.x, pStruct.y, pStruct.z = p[0], p[1], p[2]

    p1Struct = Point3DStruct()
    p1Struct.x, p1Struct.y, p1Struct.z = p1[0], p1[1], p1[2]

    p2Struct = Point3DStruct()
    p2Struct.x, p2Struct.y, p2Struct.z = p2[0], p2[1], p2[2]

    p3Struct = Point3DStruct()
    p3Struct.x, p3Struct.y, p3Struct.z = p3[0], p3[1], p3[2]

    e1Struct = Point3DStruct()
    e1Struct.x, e1Struct.y, e1Struct.z = e1[0], e1[1], e1[2]

    e2Struct = Point3DStruct()
    e2Struct.x, e2Struct.y, e2Struct.z = e2[0], e2[1], e2[2]

    e3Struct = Point3DStruct()
    e3Struct.x, e3Struct.y, e3Struct.z = e3[0], e3[1], e3[2]

    out = npct.as_ctypes(zeros(3).astype(double))

    lib.doublet(pStruct, p1Struct, p2Struct, p3Struct, e1Struct, e2Struct, e3Struct, area, maxDistance, out)
    
    return doublePointerToArray(out, 3)

def line(p: Point3D, p1: Point3D, p2: Point3D) -> ndarray:
    
    lib = npct.load_library('lib', './')

    # Define the return type of the C function
    lib.line.restype = ctypes.c_void_p

    # Define arguments of the C function
    lib.line.argtypes = [
        Point3DStruct,
        Point3DStruct,
        Point3DStruct,
        ctypes.POINTER(ctypes.c_double),
    ]

    pStruct = Point3DStruct()
    pStruct.x, pStruct.y, pStruct.z = p[0], p[1], p[2]

    p1Struct = Point3DStruct()
    p1Struct.x, p1Struct.y, p1Struct.z = p1[0], p1[1], p1[2]

    p2Struct = Point3DStruct()
    p2Struct.x, p2Struct.y, p2Struct.z = p2[0], p2[1], p2[2]

    out = npct.as_ctypes(zeros(3).astype(double))

    lib.line(pStruct, p1Struct, p2Struct, out)
    
    return doublePointerToArray(out, 3)

def source_test():

    e1 = array([1, 0, 0])
    e2 = array([0, 1, 0])
    e3 = array([0, 0, 1])

    p1 = array([0, 0])
    p2 = array([1, 0])
    p3 = array([0, 1])
    center = (1 / 3) * (p1 + p2 + p3)
    p1 = p1 - center
    p2 = p2 - center
    p3 = p3 - center

    area = 0.1
    diameter = (4 * area / pi) ** 0.5
    maxDistance = 5 * diameter

    n = 50

    x = linspace(-2, 2, num=n)
    y = linspace(-2, 2, num=n)
    z = -0.05

    X, Y = meshgrid(x, y)

    vel = zeros((n, n, 3))

    for i in range(n):
        for j in range(n):
            p = array([X[i, j], Y[i, j], z])
            vel[i, j, :] = source(p, p1, p2, p3, e1, e2, e3, area, maxDistance)
    
    plt.figure()
    plt.contourf(X, Y, vel[:, :, 2])
    plt.axis('equal')
    plt.colorbar()

    plt.figure()
    plt.quiver(X, Y, vel[:, :, 0], vel[:, :, 1])
    plt.axis('equal')

    plt.show()

    return

def line_test():

    p1 = array([0, 0, 0])
    p2 = array([1, 0, 0])

    n = 100

    x = linspace(-2, 2, num=n)
    y = linspace(-2, 2, num=n)
    z = 0.05

    X, Y = meshgrid(x, y)

    vel = zeros((n, n, 3))

    for i in range(n):
        for j in range(n):
            p = array([X[i, j], Y[i, j], z])
            vel[i, j, :] = line(p, p1, p2)
    
    plt.figure()
    plt.contourf(X, Y, vel[:, :, 2])
    plt.axis('equal')
    plt.colorbar()

    plt.figure()
    plt.quiver(X, Y, vel[:, :, 0], vel[:, :, 1])
    plt.axis('equal')

    plt.show()

    return

def doublet_test():

    e1 = array([1, 0, 0])
    e2 = array([0, 1, 0])
    e3 = array([0, 0, 1])

    p1 = array([0, 0, 0])
    p2 = array([2, 0, 0])
    p3 = array([0, 1, 0])
    center = (1 / 3) * (p1 + p2 + p3)

    p1 = p1 - center
    p2 = p2 - center
    p3 = p3 - center

    area = 0.1
    diameter = (4 * area / pi) ** 0.5
    maxDistance = 5 * diameter

    n = 100

    x = linspace(-2, 2, num=n)
    y = linspace(-2, 2, num=n)
    z = 0.1

    X, Y = meshgrid(x, y)

    vel = zeros((n, n, 3))

    for i in range(n):
        for j in range(n):
            p = array([X[i, j], Y[i, j], z])
            vel[i, j, :] = doublet(p, p1, p2, p3, e1, e2, e3, area, maxDistance)
    
    plt.figure()
    plt.contourf(X, Y, vel[:, :, 2])
    plt.axis('equal')
    plt.colorbar()

    plt.figure()
    plt.quiver(X, Y, vel[:, :, 0] + 1e-8, vel[:, :, 1] + 1e-8)
    plt.axis('equal')

    plt.show()

    return

if __name__ == '__main__':
    doublet_test()