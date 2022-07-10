from math import fabs
from tabnanny import check
from typing import Any, List
from matplotlib import bezier
from numpy import dot, linspace, zeros
from scipy.optimize import minimize, minimize_scalar

from pybird.models.types import Curve, Point
from pybird.geo.utils import vector

def quadratic(points: List[Point], t: Any) -> Any:
    """Returns a Curve or Point depending on t"""

    def func(i: float) -> Point:
        return (1 - i) * (1 - i) * points[0] + 2 * i * (1 - i) * points[1] + i * i * points[2]

    if isinstance(t, Curve):

        n = len(t)
        out = zeros((n, 3))

        for i in range(n):
            out[i, :] = func(t[i])
        
        return out

    else:

        return func(t)

def quadratic_derivative(points: List[Point], t: Any) -> Any:
    """Returns a Curve or Point depending on t"""

    def func(i: float) -> Point:
        return 2 * (1 - i) * (points[1] - points[0]) + 2 * i * (points[2] - points[1])

    if isinstance(t, Curve):

        n = len(t)
        out = zeros((n, 3))

        for i in range(n):
            out[i, :] = func(t[i])
        
        return out

    else:

        return func(t)

def cubic_derivative(points: List[Point], t: Any) -> Any:
    """Returns a Curve or Point depending on t"""

    def func(i: float) -> Point:
        return 3 * (1 - i) * (1 - i) * (points[1] - points[0]) + 6 * (1 - i) * i * (points[2] - points[1]) + 3 * i * i * (points[3] - points[2])

    if isinstance(t, Curve):

        n = len(t)
        out = zeros((n, 3))

        for i in range(n):
            out[i, :] = func(t[i])
        
        return out

    else:

        return func(t)

def cubic(points: List[Point], t: Any) -> Any:
    """Returns a Curve or Point depending on t"""

    def func(i: float) -> Point:
        return (1 - i) * (1 - i) * (1 - i) * points[0] + 3 * i * (1 - i) * (1 - i) * points[1] + 3 * i * i * (1 - i) * points[2] + i * i * i * points[3]

    if isinstance(t, Curve):

        n = len(t)
        out = zeros((n, 3))

        for i in range(n):
            out[i, :] = func(t[i])
        
        return out

    else:
        
        return func(t)

def fit_curbic_curve(curve: List[Point], t1: float, t2: float) -> List[Point]:
    """Finds the control points that best fit a cubic curve into part of a quadratic curve"""
    
    def func(x, *args):
        t, p1, p2, dp1dt, dp2dt, curve, f = args

        c1, c2 = p1 + x[0] * dp1dt, p2 + x[1] * dp2dt
        points = [p1, c1, c2, p2]
        value = [0, 0, 0]

        for i in range(len(t)):
            f2 = f(points, t[i])
            value[0] = value[0] + (curve[i][0] - f2[0]) ** 2
            value[1] = value[1] + (curve[i][1] - f2[1]) ** 2
            value[2] = value[2] + (curve[i][2] - f2[2]) ** 2

        return value[0] + value[1] + value[2]

    if len(curve) == 3:
        p1 = quadratic(curve, t1)
        p2 = quadratic(curve, t2)
        dp1dt = quadratic_derivative(curve, t1)
        dp2dt = quadratic_derivative(curve, t2)
    else:
        p1 = cubic(curve, t1)
        p2 = cubic(curve, t2)
        dp1dt = cubic_derivative(curve, t1)
        dp2dt = cubic_derivative(curve, t2)

    n = 20
    tquad = linspace(t1, t2, num=n)
    tcub = linspace(0, 1, num=n)

    if len(curve) == 3:
        ref_curve = [quadratic(curve, tquad[i]) for i in range(n)]
    else:
        ref_curve = [cubic(curve, tquad[i]) for i in range(n)]
    
    sol = minimize(func, [0, 0], args=(tcub, p1, p2, dp1dt, dp2dt, ref_curve, cubic,))

    return [p1 + sol.x[0] * dp1dt, p2 + sol.x[1] * dp2dt]

def fit_quadratic_curve(curve: List[Point], t1: float, t2: float, useLast: bool) -> Point:
    """Finds the control point that best fit a quadratic curve into part of a quadratic curve"""
    
    def func(x, *args):
        t, p1, p2, dpdt, curve, f, check = args

        c = p2 + x[0] * dpdt if check else p1 + x[0] * dpdt
        points = [p1, c, p2]
        value = [0, 0, 0]

        for i in range(len(t)):
            f2 = f(points, t[i])
            value[0] = value[0] + (curve[i][0] - f2[0]) ** 2
            value[1] = value[1] + (curve[i][1] - f2[1]) ** 2
            value[2] = value[2] + (curve[i][2] - f2[2]) ** 2

        return value[0] + value[1] + value[2]

    p1 = quadratic(curve, t1)
    p2 = quadratic(curve, t2)
    dpdt = quadratic_derivative(curve, t2) if useLast else quadratic_derivative(curve, t1)

    n = 20
    tquad = linspace(t1, t2, num=n)

    ref_curve = []
    for i in range(n):
        ref_curve.append(quadratic(curve, tquad[i]))
    

    sol = minimize(func, [0], args=(linspace(0, 1, num=n), p1, p2, dpdt, ref_curve, quadratic, useLast))

    if useLast:
        return p2 + sol.x[0] * dpdt
    else:
        return p1 + sol.x[0] * dpdt

def find_plane_intersection(curve: List[Point], p: Point) -> List[Point]:
    """Finds the intersection between the Bezier curve and a plane with p beeing a point in the plane"""
    
    def func(x, *args):
        p, curve, f, dfdt = args

        aux = f(curve, x)
        der = dfdt(curve, x)

        v = p - aux

        val = fabs(dot(v / vector.norm(v), der / vector.norm(der)))

        return val
    
    f = cubic if len(curve) == 4 else quadratic
    dfdt = cubic_derivative if len(curve) == 4 else quadratic_derivative
    sol = minimize_scalar(func, method='bounded', bounds=(0, 1), args=(p, curve, f, dfdt))

    return sol.x

def point_is_contained(curve: List[Point], p: Point) -> float:

    def func(x, *args):
        f, p, curve = args
        pCurve = f(curve, x)
        return (pCurve[0] - p[0]) ** 2 + (pCurve[1] - p[1]) ** 2 + (pCurve[2] - p[2]) ** 2
    
    f = cubic if len(curve) == 4 else quadratic
    sol = minimize_scalar(func, method='bounded', bounds=(0, 1), args=(f, p, curve))

    pCurve = f(curve, sol.x)

    dist = ((pCurve[0] - p[0]) ** 2 + (pCurve[1] - p[1]) ** 2 + (pCurve[2] - p[2]) ** 2) ** 0.5
    
    return dist