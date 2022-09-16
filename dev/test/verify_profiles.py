import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt('./src/file2.txt')

    plt.figure()
    plt.plot(data[:, 1], data[:, 0], label='U')
    plt.plot(data[:, 2], data[:, 0], label='W')
    plt.plot(data[:, 3], data[:, 0], label='S')
    plt.plot(data[:, 4], data[:, 0], label='T')
    plt.plot(data[:, 5], data[:, 0], label='R')
    plt.plot(data[:, 6], data[:, 0], label='dUdeta')
    plt.plot(data[:, 7], data[:, 0], label='dWdeta')
    plt.grid()
    plt.legend()

    plt.show()