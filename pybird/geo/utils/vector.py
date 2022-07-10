from math import acos, fabs, pi
from typing import List
from numpy import cross, deg2rad, dot, zeros
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R

from pybird.models.types import Point

def unary(a: Point) -> Point:
    """Returns a unary vector"""
    lenght = norm(a)
    if -1e-8 < lenght < 1e-8:
        return zeros(len(a))
    
    return a / lenght

def angle(a: Point, b: Point) -> List:
    """Returns the angle between two Points and the vector responsible to rotate a to b"""
    a = unary(a)
    b = unary(b)
    val = dot(a, b)
    if val > 1: val = 1
    if val < -1: val = -1
    theta = acos(val) * 180 / pi
    axis = unary(cross(a, b))
    return [theta, axis]

def rot(a: Point, theta: float, axis: Point) -> Point:
    r = R.from_rotvec(deg2rad(theta) * axis)
    return r.apply(a)