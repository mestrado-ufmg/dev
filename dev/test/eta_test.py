import numpy as np
import matplotlib.pyplot as plt

def solve(S: float, a0: float, n: int, r0: float) -> float:

    # Convergence
    tol = 1e-8

    f = lambda x: a0 * (1 - x ** n) / (1 - x)  - S 

    # Check convergence
    if 1 - 1e-8 < r0 < 1 + 1e-8: r0 = 1 + 1e-8

    aux = f(r0)

    if - tol < aux < tol:
        return r0

    # Find edges
    step = 1e-2

    if aux < 0:
        fa = aux
        a = r0
        while aux <= 0:
            r0 = r0 + step
            if 1 - 1e-8 < r0 < 1 + 1e-8: r0 = 1 + 1e-8
            aux = f(r0)
            if aux < 0:
                fa = aux
                a = r0
        b = r0
        fb = aux
    else:
        fb = aux
        b = r0
        while aux >= 0:
            r0 = r0 - step
            if 1 - 1e-8 < r0 < 1 + 1e-8: r0 = 1 + 1e-8
            aux = f(r0)
            if aux > 0:
                fb = aux
                b = r0
        a = r0
        fa = aux

    # Find root
    for i in range(100):

        dfdr = (fb - fa) / (b - a)
        dx = - fa / dfdr
        aux = a + dx
        faux = f(aux)

        print(i, faux)

        if abs(faux) < tol:
            out = aux
            break
        
        if faux < 0:
            a = aux
            fa = faux
        else:
            b = aux
            fb = faux

    return out

if __name__ == '__main__':

    S = 100
    a0 = 1
    n = 100
    r0 = 1

    y_plus = []
    eta = []

    if S / n > a0:

        res = solve(S, a0, n, r0)

        for i in range(n):
            if i == 0:
                y_plus.append(0.0)
                eta.append(0.0)
            elif i == n - 1:
                y_plus.append(S)
                eta.append(1)
            else:
                y_plus.append(y_plus[i - 1] + a0 * res ** (i - 1))
                eta.append(y_plus[i] / S)
    
    else:

        res = -1

        eta = [i * 1 / (n - 1) for i in range(n)]
        y_plus = [i * S / (n - 1) for i in range(n)]

    print([res, eta[-1], y_plus[-1]])

    plt.figure()
    plt.plot(eta, y_plus)
    plt.scatter(eta, y_plus)
    plt.grid()
    plt.show()