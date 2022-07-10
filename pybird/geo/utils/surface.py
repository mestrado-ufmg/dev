from matplotlib.pyplot import axis
from numpy import cross, dot, linspace, zeros

from pybird.models.types import Curve, Point
from pybird.geo.utils import vector

def interpolate_curve(p1: Point, p2: Point, left: Curve, right: Curve) -> Curve:
    """Interpolate a curve at a given index from bottom to top, using left and right as bounds"""

    # Bottom and top base vectors
    nSide = len(left)
    start1, start2 = left[0, :], left[nSide - 1, :]

    x1 = vector.unary(right[0, :] - left[0, :])
    y1 = vector.unary(cross(x1, p1 - start1))
    z1 = vector.unary(cross(x1, y1))

    x2 = vector.unary(right[nSide - 1, :] - left[nSide - 1, :])
    y2 = vector.unary(cross(x2, p2 - start2))
    z2 = vector.unary(cross(x2, y2))

    # Index values
    p1x, p1z = dot(p1 - start1, x1), dot(p1 - start1, z1)
    p2x, p2z = dot(p2 - start2, x2), dot(p2 - start2, z2)

    # Find the angle to rotate z
    rotData = vector.angle(x1, x2)
    zAux = vector.rot(z1, rotData[0], rotData[1])
    zAngle, axis = vector.angle(zAux, z2)
    rotSign = 1 if dot(axis, x2) > 0 else -1

    # Output
    nPoints = nSide - 2
    points = zeros((nPoints, 3))

    for i in range(nPoints):

        i1 = ((i + 1) / (nPoints))
        i2 = (1 - i1)

        # Section base vectors
        x = vector.unary(right[i + 1, :] - left[i + 1, :])
        rotData = vector.angle(x1, x)
        zAux = vector.rot(z1, rotData[0], rotData[1])
        z = vector.rot(zAux, zAngle * i1, rotSign * x)

        # Save
        points[i, :] = i2 * (p1x * x + p1z * z) + i1 * (p2x * x + p2z * z) + left[i + 1, :]
    
    return points

def interpolate_tip_curve(p1: Point, p2: Point) -> Curve:
    """Interpolate a curve at a given index from bottom to top, using left and right as bounds"""

    # Output
    nPoints = 20
    t = linspace(0, 1, num=nPoints + 2)
    points = zeros((nPoints, 3))

    for i in range(nPoints):
        points[i, :] = p2 * t[i + 1] + p1 * (1 - t[i + 1])
    
    return points