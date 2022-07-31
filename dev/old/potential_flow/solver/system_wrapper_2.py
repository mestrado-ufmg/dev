import ctypes
from typing import List, Tuple
import numpy.ctypeslib as npct
from numpy import double, int32, ndarray, zeros
from numpy import ctypeslib

# Types
INT = ctypes.c_int
PINT = ctypes.POINTER(INT)
PPINT = ctypes.POINTER(PINT)
DOUBLE = ctypes.c_double
PDOUBLE = ctypes.POINTER(DOUBLE)
PPDOUBLE = ctypes.POINTER(PDOUBLE)
PPPDOUBLE = ctypes.POINTER(PPDOUBLE)

# Define the input and output
class InputStruct(ctypes.Structure):
    _fields_ = [
        ('nf', INT),
        ('vertices', PPDOUBLE),
        ('faces', PPINT),
        ('facesAreas', PDOUBLE),
        ('facesMaxDistance', PDOUBLE),
        ('facesCenter', PPDOUBLE),
        ('controlPoints', PPDOUBLE),
        ('e1', PPDOUBLE),
        ('e2', PPDOUBLE),
        ('e3', PPDOUBLE),
        ('freestream', PDOUBLE),
        ('sigma', PDOUBLE),
    ]

class OutputStruct(ctypes.Structure):
    _fields_ = [
        ('matrix', PPDOUBLE),
        ('array', PDOUBLE),
        ('velMatrix', PPPDOUBLE),
        ('velArray', PPDOUBLE),
    ]

# Converters
def convertPointerToArray(ptr: ctypes.POINTER, shape: Tuple[int]):

    arr = zeros(shape=tuple(shape))

    if len(shape) == 1:

        for i in range(shape[0]):
            arr[i] = ptr[i]

    elif len(shape) == 2:

        for i in range(shape[0]):
            for j in range(shape[1]):
                arr[i, j] = ptr[i][j]

    elif len(shape) == 3:

        for i in range(shape[0]):
            for j in range(shape[1]):
                for k in range(shape[2]):
                    arr[i,j,k] = ptr[i][j][k]

    return arr

def convertArrayToPointer(arr: ndarray):

    if len(arr.shape) == 1:

        return npct.as_ctypes(arr)

    elif len(arr.shape) == 2:

        if arr.dtype == int32:
            ARR_DIMX = INT*arr.shape[1]
            ARR_DIMY = PINT*arr.shape[0]
        else:
            ARR_DIMX = DOUBLE*arr.shape[1]
            ARR_DIMY = PDOUBLE*arr.shape[0]
        arr_ptr = ARR_DIMY()
        for i, row in enumerate(arr):
            arr_ptr[i] = ARR_DIMX()
            for j, val in enumerate(row):
                arr_ptr[i][j] = val
        return arr_ptr
    
    elif len(arr.shape) == 3:

        if arr.dtype == int32:
            ARR_DIMX = INT*arr.shape[2]
            ARR_DIMY = PINT*arr.shape[1]
            ARR_DIMZ = PPINT*arr.shape[0]
        else:
            ARR_DIMX = DOUBLE*arr.shape[2]
            ARR_DIMY = PDOUBLE*arr.shape[1]
            ARR_DIMZ = PPDOUBLE*arr.shape[0]
        arr_ptr = ARR_DIMZ()
        for i, row in enumerate(arr):
            arr_ptr[i] = ARR_DIMY()
            for j, col in enumerate(row):
                arr_ptr[i][j] = ARR_DIMX()
                for k, val in enumerate(col):
                    arr_ptr[i][j][k] = val
        return arr_ptr

# Wrapper function
def create(vertices: ndarray,
           faces: ndarray,
           facesAreas: ndarray,
           facesMaxDistance: ndarray,
           facesCenter: ndarray,
           controlPoints: ndarray,
           e1: ndarray,
           e2: ndarray,
           e3: ndarray,
           freestream: ndarray,
           sigma: ndarray) -> List[ndarray]:

    print('0')
    
    # Create structures
    nf = len(controlPoints[:, 0])
    # inputStruct = InputStruct()
    # inputStruct.nf = nf
    # inputStruct.vertices = convertArrayToPointer(vertices.astype(double))
    # inputStruct.faces = convertArrayToPointer(faces.astype(int32))
    # inputStruct.facesAreas = convertArrayToPointer(facesAreas.astype(double))
    # inputStruct.facesMaxDistance = convertArrayToPointer(facesMaxDistance.astype(double))
    # inputStruct.facesCenter = convertArrayToPointer(facesCenter.astype(double))
    # inputStruct.controlPoints = convertArrayToPointer(controlPoints.astype(double))
    # inputStruct.e1 = convertArrayToPointer(e1.astype(double))
    # inputStruct.e2 = convertArrayToPointer(e2.astype(double))
    # inputStruct.e3 = convertArrayToPointer(e3.astype(double))
    # inputStruct.freestream = convertArrayToPointer(freestream.astype(double))
    # inputStruct.sigma = convertArrayToPointer(sigma.astype(double))

    # outputStruct = OutputStruct()
    # outputStruct.matrix = convertArrayToPointer(zeros((nf, nf)).astype(double))
    # outputStruct.array = convertArrayToPointer(zeros(nf).astype(double))
    # outputStruct.velMatrix = convertArrayToPointer(zeros((nf, nf, 3)).astype(double))
    # outputStruct.velArray = convertArrayToPointer(zeros((nf, 3)).astype(double))

    matrix = ctypeslib.as_ctypes(zeros((nf, nf)).astype(double))
    array = ctypeslib.as_ctypes(zeros(nf).astype(double))
    velMatrix = ctypeslib.as_ctypes(zeros((nf, nf, 3)).astype(double))
    velArray = ctypeslib.as_ctypes(zeros((nf, 3)).astype(double))

    # Load the compiled library
    lib = npct.load_library('lib', './')

    # Define the return type of the C function
    lib.create.restype = ctypes.c_void_p

    # Define arguments of the C function
    # lib.create.argtypes = [
    #     ctypes.c_int,
    #     ctypeslib.ndpointer(dtype=double, ndim=2, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=int32, ndim=2, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=1, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=1, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=2, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=2, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=2, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=2, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=2, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=1, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=1, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=2, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=1, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=3, flags='C_CONTIGUOUS'),
    #     ctypeslib.ndpointer(dtype=double, ndim=2, flags='C_CONTIGUOUS'),
    # ]

    print('1')

    lib.create(
        nf,
        ctypeslib.as_ctypes(vertices.astype(double)),
        ctypeslib.as_ctypes(faces.astype(int32)),
        ctypeslib.as_ctypes(facesAreas.astype(double)),
        ctypeslib.as_ctypes(facesMaxDistance.astype(double)),
        ctypeslib.as_ctypes(facesCenter.astype(double)),
        ctypeslib.as_ctypes(controlPoints.astype(double)),
        ctypeslib.as_ctypes(e1.astype(double)),
        ctypeslib.as_ctypes(e2.astype(double)),
        ctypeslib.as_ctypes(e3.astype(double)),
        ctypeslib.as_ctypes(freestream.astype(double)),
        ctypeslib.as_ctypes(sigma.astype(double)),
        matrix,
        array,
        velMatrix,
        velArray,
    )

    # lib.create(
    #     nf,
    #     ctypeslib.as_ctypes(vertices.astype(double)),
    #     faces.astype(int32),
    #     facesAreas.astype(double),
    #     facesMaxDistance.astype(double),
    #     facesCenter.astype(double),
    #     controlPoints.astype(double),
    #     e1.astype(double),
    #     e2.astype(double),
    #     e3.astype(double),
    #     freestream.astype(double),
    #     sigma.astype(double),
    #     matrix,
    #     array,
    #     velMatrix,
    #     velArray,
    # )

    print('2')

    return [
        matrix,
        array,
        velMatrix,
        velArray,
    ]

    # return [
    #     convertPointerToArray(outputStruct.matrix, [nf, nf]),
    #     convertPointerToArray(outputStruct.array, [nf]),
    #     convertPointerToArray(outputStruct.velMatrix, [nf, nf, 3]),
    #     convertPointerToArray(outputStruct.velArray, [nf, 3]),
    # ]