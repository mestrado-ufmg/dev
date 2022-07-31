import numpy as np
from time import time

from lib_wrapper import func, solve

def pass_data():
    input = np.random.random((5000, 5000))
    
    t1 = time()
    # out1 = func(input, 1) # npct.as_ctypes e converte 1D para 2D em c
    t2 = time()

    print(t2 - t1)

    # out2 = func(input, 2) # convertArrayToPointer
    t3 = time()

    print(t3 - t2)

    out3 = func(input, 3) # ctypeslib.ndpointer e passa 1D
    t4 = time()

    print(t4 - t3)

    # print('out1: {}, vals: {}'.format(out1.shape, out1[0, :3]))
    # print('out2: {}, vals: {}'.format(out2.shape, out2[0, :3]))
    print('out3: {}, vals: {}'.format(out3.shape, out3[0, :3]))

    return

def solve_system():
    n = 1500
    a = np.random.random((n, n))
    b = np.random.random(n)

    # t1 = time()
    # sol1 = solve(a, b, type=1) # convertArrayToPointer
    # t2 = time()

    # print('Time: {} s'.format(t2 - t1))
    aux = a.reshape(-1)
    t3 = time()
    sol2 = solve(aux, b, type=2) # ctypeslib.ndpointer
    t4 = time()

    print('Time: {} s'.format(t4 - t3))

    t5 = time()
    sol3 = np.linalg.solve(a, b) # linalog
    t6 = time()

    print('Time: {} s'.format(t6 - t5))

    # print(np.allclose(sol1, sol2), np.max(sol1 - sol2))
    print(np.allclose(sol2, sol3), np.max(sol2 - sol3))

if __name__ == '__main__':
    pass_data()