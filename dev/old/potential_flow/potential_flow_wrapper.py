import ctypes
import os
from typing import List
import numpy.ctypeslib as npct
from numpy import double, ndarray, zeros, random

# Types
INT = ctypes.c_int
DOUBLE = ctypes.c_double
PDOUBLE = ctypes.POINTER(DOUBLE)
PPDOUBLE = ctypes.POINTER(PDOUBLE)

# Define the input and output
class InputStruct(ctypes.Structure):
    _fields_ = [
        ('n', INT),
        ('controlPoints', PPDOUBLE),
        ('e1', PPDOUBLE),
        ('e2', PPDOUBLE),
        ('e3', PPDOUBLE),
        ('p1', PPDOUBLE),
        ('p2', PPDOUBLE),
        ('p3', PPDOUBLE),
        ('freestream', PDOUBLE),
        ('density', DOUBLE),
        ('pressure', DOUBLE),
    ]

class OutputStruct(ctypes.Structure):
    _fields_ = [
        ('sigma', PDOUBLE),
        ('doublet', PDOUBLE),
        ('velocityNorm', PDOUBLE),
        ('velocityField', PPDOUBLE),
        ('pressure', PDOUBLE),
    ]

# Converters
def double2ArrayToPointer(arr):
    """ Converts a 2D numpy to ctypes 2D array. 
    
    Arguments:
        arr: [ndarray] 2D numpy float64 array

    Return:
        arr_ptr: [ctypes double pointer]

    """

    # Init needed data types
    ARR_DIMX = DOUBLE*arr.shape[1]
    ARR_DIMY = PDOUBLE*arr.shape[0]

    # Init pointer
    arr_ptr = ARR_DIMY()

    # Fill the 2D ctypes array with values
    for i, row in enumerate(arr):
        arr_ptr[i] = ARR_DIMX()

        for j, val in enumerate(row):
            arr_ptr[i][j] = val


    return arr_ptr


def double2pointerToArray(ptr, n, m):
    """ Converts ctypes 2D array into a 2D numpy array. 
    
    Arguments:
        arr_ptr: [ctypes double pointer]

    Return:
        arr: [ndarray] 2D numpy float64 array
        
    """

    arr = zeros(shape=(n, m))

    for i in range(n):
        for j in range(m):
            arr[i,j] = ptr[i][j]

    return arr

# Wrapper function
def solver(controlPoints: ndarray,
           e1: ndarray,
           e2: ndarray,
           e3: ndarray,
           p1: ndarray,
           p2: ndarray,
           p3: ndarray,
           freestream: ndarray,
           density: float,
           pressure: float,) -> List[ndarray]:
    
    # Create structures
    inputStruct = InputStruct()
    inputStruct.n = len(controlPoints[:, 0])
    inputStruct.controlPoints = double2ArrayToPointer(controlPoints.astype(double))
    inputStruct.e1 = double2ArrayToPointer(e1.astype(double))
    inputStruct.e2 = double2ArrayToPointer(e2.astype(double))
    inputStruct.e3 = double2ArrayToPointer(e3.astype(double))
    inputStruct.p1 = double2ArrayToPointer(p1.astype(double))
    inputStruct.p2 = double2ArrayToPointer(p2.astype(double))
    inputStruct.p3 = double2ArrayToPointer(p3.astype(double))
    inputStruct.freestream = npct.as_ctypes(freestream.astype(double))
    inputStruct.density = density
    inputStruct.pressure = pressure

    outputStruct = OutputStruct()
    outputStruct.sigma = npct.as_ctypes(random.random(inputStruct.n).astype(double))
    outputStruct.doublet = npct.as_ctypes(random.random(inputStruct.n).astype(double))
    outputStruct.velocityNorm = npct.as_ctypes(zeros(inputStruct.n).astype(double))
    outputStruct.velocityField = double2ArrayToPointer(zeros((inputStruct.n, 3)).astype(double))
    outputStruct.pressure = npct.as_ctypes(zeros(inputStruct.n).astype(double))

    # Load the compiled library
    lib = npct.load_library('lib', './pybird/solver/utils/')

    # Define the return type of the C function
    lib.myfunc.restype = ctypes.c_void_p

    # Define arguments of the C function
    lib.myfunc.argtypes = [
        ctypes.POINTER(InputStruct),
        ctypes.POINTER(OutputStruct),
    ]

    lib.myfunc(inputStruct, outputStruct)

    return [
        outputStruct.sigma,
        outputStruct.doublet,
        outputStruct.velocityNorm,
        double2pointerToArray(outputStruct.velocityField, n=len(controlPoints[:, 0]), m=3),
        outputStruct.pressure,
    ]