import ctypes
from typing import List
import numpy.ctypeslib as npct
from numpy import double, ndarray, zeros, random

# Types
INT = ctypes.c_int
DOUBLE = ctypes.c_double
PDOUBLE = ctypes.POINTER(DOUBLE)
PPDOUBLE = ctypes.POINTER(PDOUBLE)
PPPDOUBLE = ctypes.POINTER(PPDOUBLE)

# Define the input and output
class InputStruct(ctypes.Structure):
    _fields_ = [
        ('n', INT),
        ('facesCenters', PPDOUBLE),
        ('controlPoints', PPDOUBLE),
        ('e1', PPDOUBLE),
        ('e2', PPDOUBLE),
        ('e3', PPDOUBLE),
        ('p1', PPDOUBLE),
        ('p2', PPDOUBLE),
        ('p3', PPDOUBLE),
        ('freestream', PDOUBLE),
    ]

class OutputStruct(ctypes.Structure):
    _fields_ = [
        ('matrix', PPDOUBLE),
        ('array', PDOUBLE),
        ('velocityCoefs', PPPDOUBLE),
    ]

class OutputVelStruct(ctypes.Structure):
    _fields_ = [
        ('vel', PPDOUBLE),
        ('velNorm', PDOUBLE),
        ('pressure', PDOUBLE),
    ]

# Converters
def double1spointerToArray(ptr, n, m):
    """ Converts ctypes 1D array into a 1D numpy array. 
    
    Arguments:
        arr_ptr: [ctypes double pointer]

    Return:
        arr: [ndarray] 1D numpy float64 array
        
    """

    arr = zeros(shape=(n, m))

    for i in range(n):
        for j in range(m):
            arr[i,j] = ptr[i][j]

    return arr

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

def double3ArrayToPointer(arr):
    """ Converts a 3D numpy to ctypes 3D array. 
    
    Arguments:
        arr: [ndarray] 3D numpy float64 array

    Return:
        arr_ptr: [ctypes double pointer]

    """

    # Init needed data types
    ARR_DIMX = DOUBLE*arr.shape[2]
    ARR_DIMY = PDOUBLE*arr.shape[1]
    ARR_DIMZ = PPDOUBLE*arr.shape[0]

    # Init pointer
    arr_ptr = ARR_DIMZ()

    # Fill the 2D ctypes array with values
    for i, row in enumerate(arr):
        arr_ptr[i] = ARR_DIMY()

        for j, col in enumerate(row):
            arr_ptr[i][j] = ARR_DIMX()

            for k, val in enumerate(col):
                arr_ptr[i][j][k] = val

    return arr_ptr

def double3pointerToArray(ptr, n, m, o):
    """ Converts ctypes 3D array into a 3D numpy array. 
    
    Arguments:
        arr_ptr: [ctypes double pointer]

    Return:
        arr: [ndarray] 3D numpy float64 array
        
    """

    arr = zeros(shape=(n, m, o))

    for i in range(n):
        for j in range(m):
            for k in range(o):
                arr[i,j,k] = ptr[i][j][k]

    return arr

# Wrapper function
def create(controlPoints: ndarray,
           e1: ndarray,
           e2: ndarray,
           e3: ndarray,
           p1: ndarray,
           p2: ndarray,
           p3: ndarray,
           freestream: ndarray) -> List[ndarray]:
    
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

    outputStruct = OutputStruct()
    outputStruct.matrix = double2ArrayToPointer(zeros((inputStruct.n, inputStruct.n)).astype(double))
    outputStruct.array = npct.as_ctypes(zeros(inputStruct.n).astype(double))
    outputStruct.velocityCoefs = double3ArrayToPointer(zeros((inputStruct.n, inputStruct.n, 3)).astype(double))

    # Load the compiled library
    lib = npct.load_library('source', './pybird/solver/utils/')

    # Define the return type of the C function
    lib.create.restype = ctypes.c_void_p

    # Define arguments of the C function
    lib.create.argtypes = [
        ctypes.POINTER(InputStruct),
        ctypes.POINTER(OutputStruct),
    ]

    lib.create(inputStruct, outputStruct)

    return [
        double2pointerToArray(outputStruct.matrix),
        double1spointerToArray(outputStruct.array),
        double3pointerToArray(outputStruct.velocityCoefs),
    ]

# Wrapper function
def parameters(controlPoints: ndarray,
               e1: ndarray,
               e2: ndarray,
               e3: ndarray,
               p1: ndarray,
               p2: ndarray,
               p3: ndarray,
               freestream: ndarray) -> List[ndarray]:
    
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

    outputStruct = OutputVelStruct()
    outputStruct.vel = double2ArrayToPointer(zeros((inputStruct.n, 3)).astype(double))
    outputStruct.velNorm = npct.as_ctypes(zeros(inputStruct.n).astype(double))
    outputStruct.pressure = npct.as_ctypes(zeros(inputStruct.n).astype(double))

    # Load the compiled library
    lib = npct.load_library('source', './pybird/solver/utils/')

    # Define the return type of the C function
    lib.parameters.restype = ctypes.c_void_p

    # Define arguments of the C function
    lib.parameters.argtypes = [
        ctypes.POINTER(InputStruct),
        ctypes.POINTER(OutputVelStruct),
    ]

    lib.parameters(inputStruct, outputStruct)

    return [
        double2pointerToArray(outputStruct.vel),
        double1spointerToArray(outputStruct.velNorm),
        double1spointerToArray(outputStruct.pressure),
    ]