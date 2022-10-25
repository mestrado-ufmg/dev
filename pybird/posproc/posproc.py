from typing import Iterable
from pybird.mesh.mesh import Mesh
from pybird.solver.solver import Solver
from numpy import asarray, cross, dot, double, ndarray, zeros

class Posproc:

    def __init__(self, mesh: Mesh, solver: Solver) -> None:
        self._mesh = mesh
        self._solver = solver
        
        self.CF: ndarray = None     # force coefficient
        self.CM: float = None       # moment coefficient

        return
        
    def forcesMoments(self, center: Iterable[float] = (.0, .0, .0),
                            area: float = 1.0,
                            chord: float = 1.0) -> None:

        self.center = asarray(center)

        forceVector = zeros(3, dtype=double)
        momentVector = zeros(3, dtype=double)

        for face in range(self._mesh.faces.shape[0]):
            forceVector = forceVector - self._mesh.e3[face, :] * self._solver.facesParams.cp[face] * self._mesh.facesAreas[face]
            momentVector = momentVector + cross(self._mesh.facesCenter[face, :] - self.center, - self._mesh.e3[face, :] * self._solver.facesParams.cp[face] * self._mesh.facesAreas[face])
        
        self.CF = forceVector / area
        self.CM = momentVector / (area * chord)
        
        return