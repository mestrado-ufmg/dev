import numpy as np

if __name__ == '__main__':
    a = np.array([1, 0, 2, 3, 5, 4, 0])
    print(a)
    a[(a >= 2) & (a < 4)] = -1
    print(a)