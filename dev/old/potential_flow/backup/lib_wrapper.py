import ctypes
import numpy as np

def func(a: int, b: np.ndarray) -> np.ndarray:

    result = np.zeros(a).astype(np.float64)
    b = b.astype(np.float64)

    lib = ctypes.CDLL('./dev/potential_flow/lib.so')
    lib.func(ctypes.c_int(a), b.ctypes.data, result.ctypes.data)

    return result

def func2(a: int, b: np.ndarray) -> np.ndarray:

    result = np.zeros_like(b).astype(np.double)
    b = b.astype(np.double)

    lib = ctypes.CDLL('./dev/potential_flow/lib.so')
    lib.func2.argtypes = [
        ctypes.c_int,
        np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ]
    lib.func2.restype = None
    lib.func2(a * a, b, result)

    return result

def func3(a: int, b: np.ndarray) -> np.ndarray:

    result = np.zeros_like(b).astype(np.double)
    b = b.astype(np.double)

    lib = ctypes.CDLL('./dev/potential_flow/lib.so')
    lib.func2.argtypes = [
        ctypes.c_int,
        np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
        np.ctypeslib.ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
    ]

    lib.func2.restype = None
    lib.func2(a, b, result)

    return result