from dataclasses import dataclass
from numpy import concatenate, cross, zeros
from sympy import E1
from pybird.geo.geo import Geo
from pybird.mesh.utils.surface_mesh import SurfaceMesh
from pybird.mesh.utils.wake_mesh import build_wake
from pybird.models.wake_model import WakeModel

@dataclass
class _Wake:
    leftWing: WakeModel
    rightWing: WakeModel
    tail: WakeModel

class Mesh:

    def __init__(self, geo: Geo) -> None:
        self._geo = geo

        self.vertices = None
        self.edges = None
        self.faces = None

        self.wake = None

        self.e1 = None
        self.e2 = None
        self.e3 = None
        self.facesCenters = None
        self.controlPoints = None
        self.p1Local = None
        self.p2Local = None
        self.p3Local = None
        self.controlPointsDistance = None
        return
    
    def build(self, size: float = None,
                    n_wing_le: float = None,
                    n_wing_te: float = None,
                    n_head: float = None,
                    n_tail_le: float = None,
                    n_tail_te: float = None,
                    wake_dist: float = None,
                    accom_dist: float = None,
                    alpha: float = None,
                    beta: float = None,) -> None:
        print('- Building mesh')
        print('  > Surface')
        self._surface = SurfaceMesh(self._geo, size, n_wing_le, n_wing_te, n_head, n_tail_le, n_tail_te)
        self.vertices, self.edges, self.faces, self.leftWingFirstSectionFacesTags, self.leftWingSecondSectionFacesTags, self.leftWingThirdSectionFacesTags, self.rightWingFirstSectionFacesTags, self.rightWingSecondSectionFacesTags, self.rightWingThirdSectionFacesTags, self.bodyFacesTags, self.headFacesTags, self.tailFacesTags = self._surface.build()
        print('  > Wake')
        leftWingWake, rightWingWake, tailWake = build_wake(self.vertices, self.edges, self.faces, self._geo, wake_dist, accom_dist, alpha, beta)
        self.wake = _Wake(leftWingWake, rightWingWake, tailWake)
        print('  > Posproc')
        self._posproc()
        return
    
    def _posproc(self) -> None:

        eps = 1e-5
        n = len(self.faces[:, 0])

        # Base vectors
        # e1 = self.vertices[self.faces[:, 1], :] - self.vertices[self.faces[:, 0]]
        # norms = (e1[:, 0] ** 2 + e1[:, 1] ** 2 + e1[:, 2] ** 2) ** 0.5
        # e1[:, 0], e1[:, 1], e1[:, 2] = e1[:, 0] / norms, e1[:, 1] / norms, e1[:, 2] / norms

        # e3 = cross(e1, self.vertices[self.faces[:, 2], :] - self.vertices[self.faces[:, 0], :])
        # norms = (e3[:, 0] ** 2 + e3[:, 1] ** 2 + e3[:, 2] ** 2) ** 0.5
        # e3[:, 0], e3[:, 1], e3[:, 2] = e3[:, 0] / norms, e3[:, 1] / norms, e3[:, 2] / norms

        # e2 = cross(e3, e1)

        # self.e1 = e1
        # self.e2 = e2
        # self.e3 = e3

        self.e1 = zeros((n, 3))
        self.e2 = zeros((n, 3))
        self.e3 = zeros((n, 3))
        self.facesCenters = zeros((n, 3))
        self.controlPoints = zeros((n, 3))

        for i in range(n):

            # Face
            face = self.faces[i, :]

            p0 = self.vertices[face[0], :]
            p1 = self.vertices[face[1], :]
            p2 = self.vertices[face[2], :]

            # Base system
            vec1 = p1 - p0
            vec2 = p2 - p0
            vec3 = cross(vec1, vec2)

            vec1Norm = (vec1[0] ** 2 + vec1[1] ** 2 + vec1[2] ** 2) ** 0.5
            vec3Norm = (vec3[0] ** 2 + vec3[1] ** 2 + vec3[2] ** 2) ** 0.5

            e1 = vec1 / vec1Norm
            e3 = vec3 / vec3Norm
            e2 = cross(e3, e1)

            self.e1[i, :] = e1
            self.e2[i, :] = e2
            self.e3[i, :] = e3

            # Face center
            self.facesCenters[i, :] = (p0 + p1 + p2) / 3
            self.controlPoints[i, :] = self.facesCenters[i, :] + eps * self.e3[i, :]
            

        # # Face centers
        # self.facesCenters = (self.vertices[self.faces[:, 0], :] + self.vertices[self.faces[:, 1], :] + self.vertices[self.faces[:, 2], :]) / 3

        # # Control points
        # self.controlPoints = self.facesCenters + eps * self.e3

        # Vertices in local system
        p1 = self.vertices[self.faces[:, 0], :] - self.facesCenters
        p2 = self.vertices[self.faces[:, 1], :] - self.facesCenters
        p3 = self.vertices[self.faces[:, 2], :] - self.facesCenters

        p1LocalX = self.e1[:, 0] * p1[:, 0] + self.e1[:, 1] * p1[:, 1] + self.e1[:, 2] * p1[:, 2]
        p1LocalY = self.e2[:, 0] * p1[:, 0] + self.e2[:, 1] * p1[:, 1] + self.e2[:, 2] * p1[:, 2]
        p1LocalZ = self.e3[:, 0] * p1[:, 0] + self.e3[:, 1] * p1[:, 1] + self.e3[:, 2] * p1[:, 2]
        self.p1Local = concatenate([p1LocalX.reshape((n, 1)), p1LocalY.reshape((n, 1)), p1LocalZ.reshape((n, 1))], axis=1)

        p2LocalX = self.e1[:, 0] * p2[:, 0] + self.e1[:, 1] * p2[:, 1] + self.e1[:, 2] * p2[:, 2]
        p2LocalY = self.e2[:, 0] * p2[:, 0] + self.e2[:, 1] * p2[:, 1] + self.e2[:, 2] * p2[:, 2]
        p2LocalZ = self.e3[:, 0] * p2[:, 0] + self.e3[:, 1] * p2[:, 1] + self.e3[:, 2] * p2[:, 2]
        self.p2Local = concatenate([p2LocalX.reshape((n, 1)), p2LocalY.reshape((n, 1)), p2LocalZ.reshape((n, 1))], axis=1)

        p3LocalX = self.e1[:, 0] * p3[:, 0] + self.e1[:, 1] * p3[:, 1] + self.e1[:, 2] * p3[:, 2]
        p3LocalY = self.e2[:, 0] * p3[:, 0] + self.e2[:, 1] * p3[:, 1] + self.e2[:, 2] * p3[:, 2]
        p3LocalZ = self.e3[:, 0] * p3[:, 0] + self.e3[:, 1] * p3[:, 1] + self.e3[:, 2] * p3[:, 2]
        self.p3Local = concatenate([p3LocalX.reshape((n, 1)), p3LocalY.reshape((n, 1)), p3LocalZ.reshape((n, 1))], axis=1)

        # Distance between one control point to another
        self.controlPointsDistance = zeros((n, n, 3))
        for i in range(n):

            x = self.controlPoints[:, 0] - self.facesCenters[i, 0]
            y = self.controlPoints[:, 1] - self.facesCenters[i, 1]
            z = self.controlPoints[:, 2] - self.facesCenters[i, 2]

            xLocal = self.e1[:, 0] * x + self.e1[:, 1] * y + self.e1[:, 2] * z
            yLocal = self.e2[:, 0] * x + self.e2[:, 1] * y + self.e2[:, 2] * z
            zLocal = self.e3[:, 0] * x + self.e3[:, 1] * y + self.e3[:, 2] * z

            self.controlPointsDistance[i, :, 0] = xLocal
            self.controlPointsDistance[i, :, 1] = yLocal
            self.controlPointsDistance[i, :, 2] = zLocal
        
        return