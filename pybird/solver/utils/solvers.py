from typing import Callable
from numpy import ndarray
from numpy.linalg import solve
from scipy.optimize import fsolve, minimize

def linear_system(a: ndarray, b: ndarray) -> ndarray:
    """
    Solves the linear system a * x = b

    Parameters
    ----------
    - a: matrix
    - b: array

    Output
    ------
    - x
    """

    return solve(a, b)

def non_linear_system(func: Callable[[ndarray], ndarray], jac: Callable[[ndarray], ndarray], x0: ndarray) -> ndarray:
    """
    Solves the non linear system func(x) = 0

    Parameters
    ----------
    - fun: equations
    - jac: jacobian
    - x0: initial condition

    Output
    ------
    - x
    """

    sol = fsolve(func, x0, fprime=jac)

    return sol