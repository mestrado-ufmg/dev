from scipy.integrate import simpson
import numpy as np

def get_integrals(y: np.ndarray, u: np.ndarray, w: np.ndarray, dudy: np.ndarray, dwdy: np.ndarray, rho: np.ndarray, psi: np.ndarray, taux: np.ndarray, tauy: np.ndarray, velocity: float, density: float):
    
    M = np.array([
        simpson(y, density * velocity - rho * u),
        simpson(y, - rho * w),
    ])

    J = np.array([
        [simpson(y, density * velocity * velocity - rho * u * u), simpson(y, - rho * u * w)],
        [simpson(y, - rho * w * u), simpson(y, - rho * w * w)],
    ])

    E = np.array([
        simpson(y, density * velocity * velocity * velocity - rho * u * (u * u + w * w)),
        simpson(y, - rho * w * (u * u + w * w)),
    ])

    Q = np.array([
        simpson(y, velocity - u),
        simpson(y, - w),
    ])

    Ko = np.array([
        simpson(y, - psi * (u * u + w * w) * rho * u),
        simpson(y, - psi * (u * u + w * w) * rho * w),
    ])

    Qo = np.array([
        simpson(y, - psi * u),
        simpson(y, - psi * w),
    ])

    D = simpson(y, taux * dudy + tauy * dwdy)

    Dx = simpson(y, taux * dwdy - tauy * dudy)
    
    return