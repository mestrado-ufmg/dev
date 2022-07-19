import numpy as np
from math import log, atan
import matplotlib.pyplot as plt

"""Fonte em um painel"""
ZERO_ERROR = 1e-8
def division(a: float, b: float) -> float:
    if -ZERO_ERROR < b < ZERO_ERROR:
        sign = 1 if b > 0 else -1
        return sign * a / ZERO_ERROR
    else:
        return a / b

def source_velocity(sigma, p, p1, p2, p3):

    # Params
    r1 = ((p[0] - p1[0]) ** 2 + (p[1] - p1[1]) ** 2 + p[2] ** 2) ** 0.5
    r2 = ((p[0] - p2[0]) ** 2 + (p[1] - p2[1]) ** 2 + p[2] ** 2) ** 0.5
    r3 = ((p[0] - p3[0]) ** 2 + (p[1] - p3[1]) ** 2 + p[2] ** 2) ** 0.5
    e1 = (p[0] - p1[0]) ** 2 + p[2] ** 2
    e2 = (p[0] - p2[0]) ** 2 + p[2] ** 2
    e3 = (p[0] - p3[0]) ** 2 + p[2] ** 2
    h1 = (p[0] - p1[0]) * (p[1] - p1[1])
    h2 = (p[0] - p2[0]) * (p[1] - p2[1])
    h3 = (p[0] - p3[0]) * (p[1] - p3[1])

    # Point 1 and 2
    d12 = ((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2) ** 0.5
    m12 = division((p2[1] - p1[1]), (p2[0] - p1[0]))

    # Point 2 and 3
    d23 = ((p3[0] - p2[0]) ** 2 + (p3[1] - p2[1]) ** 2) ** 0.5
    m23 = division((p3[1] - p2[1]), (p3[0] - p2[0]))

    # Point 3 and 1
    d31 = ((p1[0] - p3[0]) ** 2 + (p1[1] - p3[1]) ** 2) ** 0.5
    m31 = division((p1[1] - p3[1]), (p1[0] - p3[0]))

    # Velocities
    ln12 = log( division((r1 + r2 - d12), (r1 + r2 + d12)) )
    ln23 = log( division((r2 + r3 - d23), (r2 + r3 + d23)) )
    ln31 = log( division((r3 + r1 - d31), (r3 + r1 + d31)) )

    mult = (sigma / (4 * np.pi))

    u = -mult * ( division((p2[1] - p1[1]), d12) * ln12 + division((p3[1] - p2[1]), d23) * ln23 + division((p1[1] - p3[1]), d31) * ln31 )
    v = mult * ( division((p2[0] - p1[0]), d12) * ln12 + division((p3[0] - p2[0]), d23) * ln23 + division((p1[0] - p3[0]), d31) * ln31 )
    w = -mult * ( atan(division(m12 * e1 - h1, z * r1)) - atan(division(m12 * e2 - h2, z * r2)) + atan(division(m23 * e2 - h2, z * r2)) - atan(division(m23 * e3 - h3, z * r3)) + atan(division(m31 * e3 - h3, z * r3)) - atan(division(m31 * e1 - h1, z * r1)) )
    
    return [u, v, w]

if __name__ == '__main__':
    n = 100
    x = np.linspace(-2, 2, num=n)
    y = np.linspace(-2, 2, num=n)

    z = -0.1
    X, Y = np.meshgrid(x, y)
    u = np.zeros((n, n))
    v = np.zeros((n, n))
    w = np.zeros((n, n))

    p1 = [0, 0]
    p2 = [1, 0]
    p3 = [1, 1]

    for i in range(n):
        for j in range(n):
            u[i, j], v[i, j], w[i, j] = source_velocity(1, [X[i, j], Y[i, j], z], p1, p2, p3)

    plt.figure()
    plt.quiver(X, Y, u, v)
    plt.axis('equal')
    plt.grid()
    plt.title('uv')

    plt.figure()
    plt.contourf(X, Y, v)
    plt.axis('equal')
    plt.grid()
    plt.colorbar()
    plt.title('v')

    plt.figure()
    plt.contourf(X, Y, w)
    plt.axis('equal')
    plt.grid()
    plt.colorbar()
    plt.title('w')

    plt.show()