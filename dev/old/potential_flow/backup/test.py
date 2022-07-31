from curses import wrapper
import numpy as np
import lib_wrapper

def func():
    n = 3
    b = np.random.random(n)
    
    sol = lib_wrapper.func(n, b)

    print(b)
    print(sol)

def func2():
    n = 3
    b = np.random.random((n, n))
    
    sol = lib_wrapper.func3(n, b)

    print('')
    print(b)
    print(sol)


if __name__ == '__main__':
    func2()