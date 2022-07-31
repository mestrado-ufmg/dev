import numpy as np
from time import time

import sol

if __name__ == '__main__':

    a = np.random.random(1000)
    b = np.random.random(1000)

    t1 = time()
    res1 = sol.func(a, b)
    t2 = time()
    res2 = a / b
    t3 = time()

    print(t2 - t1)
    print(t3 - t2)