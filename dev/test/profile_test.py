import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def sublayer(y_plus: float) -> float:
    return y_plus

def log_law_layer(y_plus: float) -> float:
    k = 0.41
    C = 5.0
    return (1 / k) * np.log(y_plus) + C

def buffer_layer(y_plus: float, y_min: float, y_max: float) -> float:

    u_plus_1 = sublayer(y_min)
    u_plus_2 = log_law_layer(y_max)

    der_u_plus_1 = 1
    der_u_plus_2 = 1 / (0.41 * y_max)

    y = (u_plus_1 - u_plus_2 - y_min * der_u_plus_1 + y_max * der_u_plus_2) / (der_u_plus_2 - der_u_plus_1)
    u = y * der_u_plus_2 + u_plus_2 - y_max * der_u_plus_2

    t = (y_plus - y_min) / (y_max - y_min)

    return (1 - t) * (1 - t) * u_plus_1 + 2 * (1 - t) * t * u + t * t * u_plus_2

def find_bezier_curve():
    
    # Interpolation limites
    y_min = 5
    y_max = 200
    u_min = sublayer(y_min)
    u_max = log_law_layer(y_max)

    # Bezier curves
    def func_bezier_order_4(u_min, y_min, u_max, y_max, x, y_plus):

        dudy_min = 10 * np.log10(y_min) * np.log(10)
        dudy_max = 1 / (0.41 * np.log10(np.e))

        c_1 = u_min + dudy_min * x[0]
        c_2 = x[1]
        c_3 = u_max - dudy_max * x[2]

        t = (np.log10(y_plus) - np.log10(y_min)) / (np.log10(y_max) - np.log10(y_min))

        return ((1 - t) ** 4) * u_min + 4 * ((1 - t) ** 3) * t * c_1 + 6 * ((1 - t) ** 2) * (t ** 2) * c_2 + 4 * (1 - t) * (t ** 3) * c_3 + (t ** 4) * u_max

    # Fit curves
    u_fit = np.linspace(u_min, u_max, num=100)
    y_fit = u_fit + np.exp(-0.41 * 5) * (np.exp(0.41 * u_fit) - 1 - 0.41 * u_fit - 0.5 * ((0.41 * u_fit) ** 2) - (1 / 6) * ((0.41 * u_fit) ** 3))

    # Objective function
    def obj_func(x, *args):
        u, y, func, u_min, y_min, u_max, y_max = args
        # f = func([u_1] + x.tolist() + [u_2], y)
        f = func(u_min, y_min, u_max, y_max, x, y)
        return np.sum(np.power(u - f, 2))
    
    # Find points
    sol = minimize(obj_func, [0.1, 0.5 * (u_min + u_max), 0.1], args=(u_fit, y_fit, func_bezier_order_4, u_min, y_min, u_max, y_max), bounds=((0, 1), (u_min, u_max), (0, 1)))

    # Create curves
    n = 1000
    y_plus = np.geomspace(0.1, 1000, num=n)
    u_plus = np.zeros(n)

    for i in range(n):
        if y_plus[i] <= y_min:
            u_plus[i] = sublayer(y_plus[i])
        elif y_min < y_plus[i] < y_max:
            u_plus[i] = func_bezier_order_4(u_min, y_min, u_max, y_max, sol.x, y_plus[i])
        else:
            u_plus[i] = log_law_layer(y_plus[i])
    
    u_plus_eq = np.linspace(0.1, 21.7, num=1000)
    y_plus_eq = u_plus_eq + np.exp(-0.41 * 5) * (np.exp(0.41 * u_plus_eq) - 1 - 0.41 * u_plus_eq - 0.5 * ((0.41 * u_plus_eq) ** 2) - (1 / 6) * ((0.41 * u_plus_eq) ** 3))

    # Show curves
    plt.figure()
    plt.plot(y_plus, u_plus)
    plt.plot(y_plus_eq, u_plus_eq, '--k')
    plt.grid()
    plt.xscale('log')
    plt.show()

    return

def test_profile():

    k = 0.41
    C = 5.0
    u_min, y_min = 5.0, 5.0
    u_max, y_max = 17.922725284263503, 200

    log_y_min, log_y_max = np.log10(y_min), np.log10(y_max)
    a = u_min + 10 * log_y_min * np.log(10) * 0.26957378
    b = 14.2135593
    c = u_max - (1 / (k * np.log10(np.e))) * 0.51958278

    print(a, b, c)

    n = 1000
    y_plus = np.geomspace(0.1, 1000, num=n)
    u_plus = np.zeros(n)

    for i in range(n):
        if y_plus[i] <= y_min:
            u_plus[i] = y_plus[i]
        elif y_min < y_plus[i] <= y_max:
            t = (np.log10(y_plus[i]) - log_y_min) / (log_y_max - log_y_min)
            f1 = ( (1 - t) ** 4 )
            f2 = 4 * ( (1 - t) ** 3 ) * t
            f3 = 6 * ( (1 - t) ** 2 ) * (t ** 2)
            f4 = 4 * (1 - t) * (t ** 3)
            f5 = t ** 4
            u_plus[i] = f1 * u_min + f2 * a + f3 * b + f4 * c + f5 * u_max
        else:
            u_plus[i] = (1 / k) * np.log(y_plus[i]) + C

    u_plus_eq = np.linspace(0.1, 21.7, num=1000)
    y_plus_eq = u_plus_eq + np.exp(-0.41 * 5) * (np.exp(0.41 * u_plus_eq) - 1 - 0.41 * u_plus_eq - 0.5 * ((0.41 * u_plus_eq) ** 2) - (1 / 6) * ((0.41 * u_plus_eq) ** 3))

    # Show curves
    plt.figure()
    plt.plot(y_plus, u_plus)
    plt.plot(y_plus_eq, u_plus_eq, '--k')
    plt.grid()
    # plt.xscale('log')
    plt.show()

    return

if __name__ == '__main__':
    # find_bezier_curve()
    test_profile()