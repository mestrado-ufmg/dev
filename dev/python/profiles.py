import numpy as np

def get_profiles(n_layers: int, delta: float, A: float, B: float, Psi: float, velocity, density, viscosity, mach):

    eta = np.linspace(0, 1, num=n_layers)

    f0 = 6 * np.power(eta, 2) - 8 * np.power(eta, 3) + 3 * np.power(eta, 4)
    f1 = eta - 3 * np.power(eta, 2) + 3 * np.power(eta, 3) - np.power(eta, 4)
    f2 = (eta - 4 * np.power(eta, 2) + 6 * np.power(eta, 3) - 4 * np.power(eta, 4) + np.power(eta, 5)) * np.power(1 - eta, 2)
    f3 = (np.power(eta, 2) - 3 * np.power(eta, 3) + 3 * np.power(eta, 4) - np.power(eta, 5)) * np.power(1 - eta, 2)

    df0_deta = 12 * eta - 24 * np.power(eta, 2) + 12 * np.power(eta, 3)
    df1_deta = 1 - 6 * eta + 9 * np.power(eta, 2) - 4 * np.power(eta, 3)
    df2_deta = (1 - 6 * eta + 9 * np.power(eta, 2) - 4 * np.power(eta, 3)) * np.power(1 - eta, 2) + 2 * (eta - 4 * np.power(eta, 2) + 6 * np.power(eta, 3) - 4 * np.power(eta, 4) + np.power(eta, 5)) * (1 - eta)
    df3_deta = (2 * eta - 9 * np.power(eta, 2) + 12 * np.power(eta, 3) - 5 * np.power(eta, 4)) * np.power(1 - eta, 2) + 2 * (np.power(eta, 2) - 3 * np.power(eta, 3) + 3 * np.power(eta, 4) - np.power(eta, 5)) * (1 - eta)

    # Velocities
    U = f0 + A * (1 - 0.6 * (A - 3) * eta * eta * eta) * f1
    W = B * f2 + Psi * f3

    # Velocity derivatives
    dU_deta = df0_deta + A * (1 - 0.6 * (A - 3) * eta * eta * eta) * df1_deta - A * 1.8 * (A - 3) * np.power(eta, 2) * f1
    dW_deta = B * df2_deta + Psi * df3_deta

    # Density
    R = 1 / (1 + 0.2 * mach * mach * (1 - U * U - W * W))

    # Angle
    psi = np.arctan(W / U)
    dpsidy = (1 / delta) * ((dW_deta * U - W * dU_deta) / (U * U)) * (1 / (1 + (W / U) * (W / U)))

    # Real values
    y = eta * delta
    u = U * velocity
    w = W * velocity
    dudy = dU_deta * velocity / delta
    dwdy = dW_deta * velocity / delta
    rho = R * density
    taux = viscosity * dudy
    tauy = viscosity * dwdy

    return [y, u, w, dudy, dwdy, rho, psi, dpsidy, taux, tauy]

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt

    delta = 0.1
    A = 1
    B = 1
    Psi = 1

    velocity = 1.0
    density = 1.0
    viscosity = 1e-5
    mach = 0.1
    
    data = get_profiles(300, delta, A, B, Psi, velocity, density, viscosity, mach)

    # plt.figure()
    # plt.plot(data[1], data[0], label='u')
    # plt.plot(data[2], data[0], label='w')
    # plt.grid()
    # plt.legend()

    # plt.figure()
    # plt.plot(data[3], data[0], label='dudy')
    # plt.plot(data[4], data[0], label='dwdy')
    # plt.grid()
    # plt.legend()

    # plt.figure()
    # plt.plot(data[8], data[0], label='taux')
    # plt.plot(data[9], data[0], label='tauy')
    # plt.grid()
    # plt.legend()

    # plt.figure()
    # plt.plot(data[5], data[0], label='rho')
    # plt.grid()
    # plt.legend()

    plt.figure()
    plt.plot(data[6], data[0], label='psi')
    plt.grid()
    plt.legend()

    plt.figure()
    plt.plot(data[7], data[0], label='dpsidy')
    plt.grid()
    plt.legend()

    plt.show()