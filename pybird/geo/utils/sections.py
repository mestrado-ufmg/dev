from typing import List
from numpy import array, cross, dot, flip, loadtxt, zeros
from numpy.linalg import norm
from pytools import argmax

from pybird.models.types import Point
from pybird.geo.utils import vector
from pybird.geo.utils.curve import interpolate_2D, interpolate_3D

def process_section(file: str, v1: Point, v2: Point, x: Point, z: Point, n: int) -> List:

    nFoil = 300

    # Load airfoil
    foil = interpolate_2D(loadtxt(file), nFoil)

    # Scale and center foil
    scale = 1 / (max(foil[:, 0]) - min(foil[:, 0]))
    foil[:, 0], foil[:, 1] = -foil[:, 0] * scale, foil[:, 1] * scale
    foil[:, 0], foil[:, 1] = foil[:, 0] - foil[0, 0],  foil[:, 1] - foil[0, 1]

    # Align foil
    chord = norm(v1 - v2)
    xAux = vector.unary(v1 - v2)
    rotData = vector.angle(x, xAux)
    x = xAux * chord
    z = vector.unary(vector.rot(z, rotData[0], rotData[1])) * chord

    # Find leading edge index
    index2 = argmax(foil[:, 0])
    index1 = round(index2 / 2)
    index4 = 3 * index1

    # Create curve
    curve = zeros((nFoil, 3))
    xCurve = vector.unary(foil[index2, :] - foil[0, :])
    zCurve = vector.unary(cross(array([0, 0, 1]), array([xCurve[0], xCurve[1], 0])))[:2]

    for i in range(nFoil):
        vec = foil[i, :]
        curve[i, :] = x * dot(vec, xCurve) + z * dot(vec, zCurve) + v2

    # Separate curves and points
    pUpper = curve[index1, :]
    pLower = curve[index4, :]

    upperSurface2 = interpolate_3D(curve[:index1 + 1, :], n + 2)[1:n + 1]
    upperSurface1 = interpolate_3D(curve[index1:index2 + 1, :], n + 2)[1:n + 1]
    lowerSurface1 = flip(interpolate_3D(curve[index2:index4 + 1, :], n + 2)[1:n + 1], axis=0)
    lowerSurface2 = flip(interpolate_3D(curve[index4:, :], n + 2)[1:n + 1], axis=0)

    return pUpper, pLower, upperSurface1, lowerSurface1, upperSurface2, lowerSurface2