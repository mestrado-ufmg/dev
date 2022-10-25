import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    data = np.loadtxt('conv.txt')

    plt.figure()
    plt.plot(data[:, 0], data[:, 1])
    plt.xlabel('Interactions')
    plt.ylabel('Error')
    plt.yscale('log')
    plt.grid()

    plt.figure()
    plt.plot(data[:, 0], data[:, 2])
    plt.xlabel('Interactions')
    plt.ylabel('Drag')
    plt.grid()

    plt.figure()
    plt.plot(data[:, 1], data[:, 2])
    plt.xlabel('Error')
    plt.ylabel('Drag')
    plt.grid()
    plt.xscale('log')

    plt.show()