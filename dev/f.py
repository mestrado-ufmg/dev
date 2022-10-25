from re import I
import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':

    n = 100
    Re_theta = np.linspace(0, 1e6, num=n)
    Hk = np.linspace(1.8, 5, num=n)

    x1, x2 = np.meshgrid(Re_theta, Hk)

    f1 = 0.01 * np.sqrt(np.power(2.4 * x2 - 3.7 + 2.5 * np.tanh(1.5 * x2 - 4.65), 2) + 0.25)
    f2 = np.power(10, (1.415 / (x2 - 1) - 0.489) * np.tanh(20 / (x2 - 1) - 12.9) + 3.295 / (x2 - 1) + 0.44)
    y = f1 * (x1 - f2)

    plt.figure()
    plt.contourf(x1, x2, y)
    plt.xlabel('Re_theta')
    plt.ylabel('Hk')
    plt.colorbar()
    plt.grid()
    plt.show()