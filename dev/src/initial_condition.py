import matplotlib.pyplot as plt
import numpy as np

def delta_calc(U, rho, mu, x):
    return 5.0 * np.sqrt(mu * x / (rho * U))

if __name__ == '__main__':
    
    U = 10.0
    rho = 1.225
    mu = 1e-5
    Cf_max = 0.02

    # delta
    x = np.linspace(0.1, 0.5, num=200)
    delta = delta_calc(U, rho, mu, x)

    # Shear stress
    tau_w = 0.332 * mu * U * np.sqrt(rho * U / (mu * x))

    # A
    A = 0.332 * 5.0 * U * np.ones(x.size)

    # Plot
    # plt.figure()
    # plt.plot(x, tau_w)
    # plt.title('tau_w')
    # plt.grid()

    plt.figure()
    plt.plot(x, A)
    plt.title('A')
    plt.grid()

    # plt.figure()
    # plt.plot(x, delta)
    # plt.title('delta')
    # plt.grid()

    plt.show()