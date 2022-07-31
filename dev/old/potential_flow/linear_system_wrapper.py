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

def double1pointerToArray(ptr, n):
    """ Converts ctypes 2D array into a 1D numpy array. 
    
    Arguments:
        arr_ptr: [ctypes double pointer]

    Return:
        arr: [ndarray] 1D numpy float64 array
        
    """

    arr = zeros(n)

    for i in range(n):
        arr[i] = ptr[i]

    return arr

# Wrapper function
def solver(a: ndarray,
           b: ndarray,
           x0: ndarray) -> List[ndarray]:
    
    # Inputs
    n = len(b)
    apt = double2ArrayToPointer(a.astype(double))
    bpt = npct.as_ctypes(b.astype(double))
    x0pt = npct.as_ctypes(x0.astype(double))

    # Load the compiled library
    lib = npct.load_library('lib', './dev/potential_flow/')

    # Define the return type of the C function
    lib.ls_wrapper.restype = ctypes.c_void_p

    # Define arguments of the C function
    lib.ls_wrapper.argtypes = [
        INT,
        PPDOUBLE,
        PDOUBLE,
    ]

    lib.ls_wrapper(n, apt, bpt, x0pt)

    return double1pointerToArray(x0pt, n)