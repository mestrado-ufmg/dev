from dataclasses import dataclass
from numpy import cross
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

        self.facesNormals = None
        self.facesCenters = None
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

        # Faces normals
        normals = cross(self.vertices[self.faces[:, 1], :] - self.vertices[self.faces[:, 0], :], self.vertices[self.faces[:, 2], :] - self.vertices[self.faces[:, 0], :])
        norms = (normals[:, 0] ** 2 + normals[:, 1] ** 2 + normals[:, 2] ** 2) ** 0.5
        normals[:, 0], normals[:, 1], normals[:, 2] = normals[:, 0] / norms, normals[:, 1] / norms, normals[:, 2] / norms
        self.facesNormals = normals

        # Face centers
        self.facesCenters = (self.vertices[self.faces[:, 0], :] + self.vertices[self.faces[:, 1], :] + self.vertices[self.faces[:, 2], :]) / 3
        
        return