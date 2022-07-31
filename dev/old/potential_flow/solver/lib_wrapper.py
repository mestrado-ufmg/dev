import ctypes
from typing import List, Tuple
from numpy import ctypeslib, float32, float64
import numpy.ctypeslib as npct
from numpy import asarray, double, empty, int32, ndarray, zeros, zeros_like

# Types
INT = ctypes.c_int
PINT = ctypes.POINTER(INT)
PPINT = ctypes.POINTER(PINT)
DOUBLE = ctypes.c_double
PDOUBLE = ctypes.POINTER(DOUBLE)
PPDOUBLE = ctypes.POINTER(PDOUBLE)
PPPDOUBLE = ctypes.POINTER(PPDOUBLE)

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

        if arr.dtype == int32:
            ARR_DIMX = INT*arr.shape[0]
        else:
            ARR_DIMX = DOUBLE*arr.shape[0]

        arr_ptr = ARR_DIMX()
        for i, val in enumerate(arr):
            arr_ptr[i] = val
        return arr_ptr

        # return npct.as_ctypes(arr)

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
def func(input: ndarray, type:int=1) -> ndarray:
    
    # Load the compiled library
    lib = ctypes.CDLL('./lib3.so')
    
    if type == 1:

        lib.func1.restype = ctypes.c_void_p

        lib.func1.argtypes = [
            INT,
            INT,
            PDOUBLE,
            PDOUBLE,
        ]

        size = input.shape[0] * input.shape[1]

        inputPointer = npct.as_ctypes(input.reshape(size))
        outputPointer = npct.as_ctypes(zeros(size).reshape(size))

        lib.func1(input.shape[0], input.shape[1], inputPointer, outputPointer)

        out = asarray(outputPointer).reshape((input.shape[0], input.shape[1]))

    elif type == 2:

        lib.func2.restype = ctypes.c_void_p

        lib.func2.argtypes = [
            INT,
            INT,
            PPDOUBLE,
            PPDOUBLE,
        ]

        lib.func2.restype = None
        
        inputPointer = convertArrayToPointer(input)
        outputPointer = convertArrayToPointer(zeros_like(input))

        lib.func2(input.shape[0], input.shape[1], inputPointer, outputPointer)

        out = convertPointerToArray(outputPointer, [input.shape[0], input.shape[1]])
    
    elif type == 3:
        
        ND_POINTER_1 = ctypeslib.ndpointer(dtype=float64, ndim=1, flags="C")
        
        lib.func3.restype = None

        lib.func3.argtypes = [
            INT,
            INT,        
            ND_POINTER_1,
            ND_POINTER_1,
        ]

        lib.func3.restype = None
        
        shape = input.shape
        input = input.reshape(-1)
        out = empty(input.shape)

        lib.func3(shape[0], shape[1], input, out)

        out = out.reshape(shape)
    
    return out


# Wrapper function
def solve(a: ndarray, b: ndarray, type: int = 1) -> ndarray:
    
    # Load the compiled library
    lib = ctypes.CDLL('./lib3.so')

    if type == 1:
        # Define the return type of the C function
        lib.solve_system.restype = None

        # Define arguments of the C function
        lib.solve_system.solve_system = [
            INT,
            PPDOUBLE,
            PDOUBLE,
            PDOUBLE,
        ]

        n = len(a[:, 0])

        aPointer = convertArrayToPointer(a)
        bPointer = convertArrayToPointer(b)
        solPointer = convertArrayToPointer(empty(n))

        lib.solve_system(n, aPointer, bPointer, solPointer)

        out = convertPointerToArray(solPointer, [n])

        return out
    
    elif type == 2:

        ND_POINTER_1 = ctypeslib.ndpointer(dtype=float32, ndim=1, flags="C")

        # Define the return type of the C function
        lib.solve_system_1D.restype = None

        # Define arguments of the C function
        lib.solve_system_1D.argtypes = [
            INT,
            ND_POINTER_1,
            ND_POINTER_1,
            ND_POINTER_1,
        ]

        n = len(b)

        # a = a.reshape(-1)
        out = empty(b.shape).astype(float32)

        lib.solve_system_1D(n, a.astype(float32), b.astype(float32), out)

        return out