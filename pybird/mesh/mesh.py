from time import time
from dataclasses import dataclass
from math import pi
from numpy import cross, ndarray
from numpy.linalg import norm

from pybird.geo.geo import Geo
from pybird.mesh.utils.surface_mesh import SurfaceMesh
from pybird.mesh.utils.wake_mesh import build_wake
from pybird.models.wake_model import WakeModel

@dataclass
class _Wake:
    leftWing: WakeModel
    rightWing: WakeModel
    tail: WakeModel

def unary(a: ndarray):

    if a.ndim == 2:
        n = a.shape[0]
        for i in range(n):
            a[i, :] = a[i, :] / norm(a[i, :])
        return a

    return a / norm(a)

class Mesh:

    def __init__(self, geo: Geo) -> None:
        self._geo = geo

        self.vertices: ndarray = None
        self.edges: ndarray = None
        self.faces: ndarray = None

        self.wake: _Wake = None

        self.e1: ndarray = None
        self.e2: ndarray = None
        self.e3: ndarray = None
        self.facesCenters: ndarray = None
        self.controlPoints: ndarray = None
        self.p1Local: ndarray = None
        self.p2Local: ndarray = None
        self.p3Local: ndarray = None
        self.controlPointsDistance: ndarray = None
        return
    
    def build(self, size: float = None,
                    n_wing_le: float = None,
                    n_wing_te: float = None,
                    n_head: float = None,
                    n_tail_le: float = None,
                    n_tail_te: float = None,
                    n_body: float = None,
                    wake_dist: float = None,
                    accom_dist: float = None,
                    alpha: float = None,
                    beta: float = None,) -> None:
                    
        print('- Building mesh')

        self._surface = SurfaceMesh(self._geo, size, n_wing_le, n_wing_te, n_head, n_tail_le, n_tail_te, n_body)
        self.vertices, self.edges, self.faces, self.leftWingFirstSectionFacesTags, self.leftWingSecondSectionFacesTags, self.leftWingThirdSectionFacesTags, self.rightWingFirstSectionFacesTags, self.rightWingSecondSectionFacesTags, self.rightWingThirdSectionFacesTags, self.bodyFacesTags, self.headFacesTags, self.tailFacesTags = self._surface.build()

        leftWingWake, rightWingWake, tailWake = build_wake(self.vertices, self.edges, self.faces, self._geo, wake_dist, accom_dist, alpha, beta)
        self.wake = _Wake(leftWingWake, rightWingWake, tailWake)

        self._posproc()
        
        return
    
    def _posproc(self) -> None:

        eps = 1e-5

        # Faces center
        self.facesCenter = (1 / 3) * (self.vertices[self.faces[:, 0], :] + self.vertices[self.faces[:, 1], :] + self.vertices[self.faces[:, 2], :])

        # Base vectors
        self.e3 = unary(cross(self.vertices[self.faces[:, 1], :] - self.vertices[self.faces[:, 0], :], self.vertices[self.faces[:, 2], :] - self.vertices[self.faces[:, 0], :]))
        self.e1 = unary(self.vertices[self.faces[:, 1], :] - self.facesCenter)
        self.e2 = unary(cross(self.e3, self.e1))
        
        # Control points
        self.controlPoints = self.facesCenter + eps * self.e3

        # Faces areas
        auxVec = cross(self.vertices[self.faces[:, 1]] - self.vertices[self.faces[:, 0]], self.vertices[self.faces[:, 2]] - self.vertices[self.faces[:, 0]])
        auxNorm = (auxVec[:, 0] ** 2 + auxVec[:, 1] ** 2 + auxVec[:, 2] ** 2) ** 0.5
        self.facesAreas = 0.5 * auxNorm

        # Faces max distance
        self.facesMaxDistance = 10 * (4 * self.facesAreas / pi) ** 0.5

        return
