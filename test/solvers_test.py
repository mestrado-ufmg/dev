import sys
sys.path.append('./')

from pybird.solver.utils.solvers import linear_system, non_linear_system

from time import time
import numpy as np

if __name__ == '__main__':

    n = 100

    for _ in range(3):

        a = np.random.random((n, n))
        b = np.random.random(n)
        
        t1 = time()
        sol1 = linear_system(a, b)
        t2 = time()

        def func(x):
            out = np.zeros_like(x)
            for i in range(n):
                out[i] = np.sum(a[i, :] * x) - b[i]
            return np.sum(out * out)
        
        def jac(x):
            out = np.zeros((n, n))
            for i in range(n):
                f = np.sum(a[i, :] * x) - b[i]
                out[i, :] = a[i, :] * f
            return out
        
        x0 = np.random.random(n)

        t3 = time()
        sol2 = non_linear_system(func, None, x0)
        t4 = time()

        is_equal = np.allclose(sol1, sol2)
        dt1 = t2 - t1
        dt2 = t4 - t3
        
        print('dt: [{:.8f}, {:.8f}] | bigger: {} | equal: {}'.format(dt1, dt2, dt2 > dt1, is_equal))