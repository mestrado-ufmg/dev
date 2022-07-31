from time import time
import sum_wrapper

def our_function(num_numbers):
    sum = 0
    for i in range(num_numbers):
        sum += i
    return sum

if __name__ == '__main__':
    n = 10000
    t1 = time()
    res1 = sum_wrapper.library(n)
    t2 = time()
    res2 = sum_wrapper.exe(n)
    t3 = time()
    res3 = our_function(n)
    t4 = time()

    print(res1, res2, res3)
    print('1: {} | 2: {} | 3: {}'.format(t2 - t1, t3 - t2, t4 - t3))