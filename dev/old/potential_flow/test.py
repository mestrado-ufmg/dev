import numpy as np
import linear_system_wrapper

if __name__ == '__main__':
    
    n = 100
    a = np.random.random((n, n))
    b = np.random.random(n)
    x0 = np.random.random(n)
    
    for i in range(n):
        a[i, i] = a[i, i] + 1

    sol1 = linear_system_wrapper.solver(a, b, x0)
    sol2 = np.linalg.solve(a, b)

    # print(sol1)
    # print(sol2)

    # print('')

    print(np.dot(a, sol1) - b)
    # print(np.dot(a, sol2) - b)