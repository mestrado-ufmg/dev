from subprocess import check_output
import ctypes

_sum = ctypes.CDLL('./dev/c_from_python/libsum.so')

def library(num_numbers):
    global _sum
    result = _sum.our_function(ctypes.c_int(num_numbers))
    return float(result)

def exe(num_numbers):
    res = check_output('./dev/c_from_python/exesum {}'.format(num_numbers), shell=True)
    return int(str(res).replace('b', '').replace('\'', ''))