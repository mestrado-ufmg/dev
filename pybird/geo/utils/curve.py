from numpy import asarray, gradient, linspace, zeros
from scipy.interpolate import splprep, splev
from scipy.integrate import simpson

from pybird.models.types import Curve

def interpolate_2D(points: Curve, n: int) -> Curve:
    tck, _ = splprep([points[:, 0], points[:, 1]], s=0)
    u = linspace(0, 1, num=n)
    tulple_points = splev(u, tck)
    new_points = zeros((n, 2))
    new_points[:, 0], new_points[:, 1] = asarray(list(tulple_points[0])), asarray(list(tulple_points[1]))
    return new_points

def interpolate_3D(points: Curve, n: int) -> Curve:
    tck, _ = splprep([points[:, 0], points[:, 1], points[:, 2]], s=0)
    u = linspace(0, 1, num=n)
    tulple_points = splev(u, tck)
    new_points = zeros((n, 3))
    new_points[:, 0], new_points[:, 1], new_points[:, 2] = asarray(list(tulple_points[0])), asarray(list(tulple_points[1])), asarray(list(tulple_points[2]))
    return new_points

def lenght(points: Curve) -> float:
    
    n = 200
    interp = interpolate_3D(points, n)

    u = linspace(0, 1, num=n)
    df = gradient(interp[:, 0], u)
    dg = gradient(interp[:, 1], u)
    dh = gradient(interp[:, 2], u)

    df[0] = df[1] + 0.5 * (df[0] - df[2])
    dg[0] = dg[1] + 0.5 * (dg[0] - dg[2])
    dh[0] = dh[1] + 0.5 * (dh[0] - dh[2])

    df[n - 1] = df[n - 2] + 0.5 * (df[n - 2] - df[n - 4])
    dg[n - 1] = dg[n - 2] + 0.5 * (dg[n - 2] - dg[n - 4])
    dh[n - 1] = dh[n - 2] + 0.5 * (dh[n - 2] - dh[n - 4])

    sol = simpson((df ** 2 + dg ** 2 + dh ** 2) ** 0.5, u)

    return sol

def interpolate_fixed_ds(points: Curve, dl: float, removeEdges: bool = False) -> Curve:

    l = lenght(points)
    n = round(l / dl)
    new_points = interpolate_3D(points, n + 1)

    if removeEdges:
        return new_points[1:n, :]
    else:
        return new_points