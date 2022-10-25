from numpy import array, deg2rad
from scipy.spatial.transform import Rotation as R

from pybird.mesh.mesh import Mesh
from pybird.solver.utils.lib_wrapper import FacesParams, VerticesParams, solve

class Solver:

    def __init__(self, mesh: Mesh) -> None:
        self._mesh = mesh
        self.freestream = 1
        
        self.facesParams: FacesParams = None
        self.verticesParams: VerticesParams = None

        return
    
    @property
    def done(self) -> bool:
        return self.facesParams is not None and self.verticesParams is not None
    
    def setEnvironment(self, freestream: float) -> None:
        self.freestream = freestream
        return
    
    def solve(self, alpha: float = None, beta: float = None) -> None:

        print('- Solving external flow')
        
        # Freestream vector
        x = array([-1.0, 0.0, 0.0])
        y = array([0.0, -1.0, 0.0])
        z = array([0.0, 0.0, 1.0])

        if alpha is not None:
            r = R.from_rotvec(-deg2rad(alpha) * y)
            x = r.apply(x)
            z = r.apply(z)

        if beta is not None:
            r = R.from_rotvec(-deg2rad(beta) * z)
            x = r.apply(x)
            y = r.apply(y)
        
        self.U = x * self.freestream

        # Solve
        data = solve(0,
                     self._mesh.vertices,
                     self._mesh.faces,
                     self._mesh.facesAreas,
                     self._mesh.facesMaxDistance,
                     self._mesh.facesCenter,
                     self._mesh.controlPoints,
                     self._mesh.p1Local, self._mesh.p2Local, self._mesh.p3Local,
                     self._mesh.e1, self._mesh.e2, self._mesh.e3,
                     self.U,
                     self._mesh.wake.leftWing.grid, self._mesh.wake.leftWing.vertices, self._mesh.wake.leftWing.faces,
                     self._mesh.wake.rightWing.grid, self._mesh.wake.rightWing.vertices, self._mesh.wake.rightWing.faces,
                     self._mesh.wake.tail.grid, self._mesh.wake.tail.vertices, self._mesh.wake.tail.faces)
        
        self.facesParams, self.verticesParams = data

        return