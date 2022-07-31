from time import time
from numpy import array, deg2rad, ndarray, zeros
from numpy.linalg import solve
from scipy.spatial.transform import Rotation as R

from pybird.mesh.mesh import Mesh
from pybird.solver.utils.system_wrapper import create

class Solver:

    def __init__(self, mesh: Mesh) -> None:
        self._mesh = mesh
        self.freestream = 1
        
        self.sigma: ndarray = None
        self.doublet: ndarray = None
        self.velNorm: ndarray = None
        self.velField: ndarray = None
        self.cp: ndarray = None

        return
    
    def setEnvironment(self, freestream: float) -> None:
        self.freestream = freestream
        return
    
    def solve(self, alpha: float = None, beta: float = None) -> None:

        print('- Solving potential flow')

        t1 = time()

        self._n = len(self._mesh.faces[:, 0])

        self.sigma = zeros(self._n)
        self.doublet = zeros(self._n)
        self.velNorm = zeros(self._n)
        self.velField = zeros((self._n, 3))
        self.cp = zeros(self._n)
        
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
        
        U = x * self.freestream
        self._UNorm2 = (U[0] ** 2 + U[1] ** 2 + U[2] ** 2) ** 0.5

        self.sigma = self._mesh.e3[:, 0] * U[0] + self._mesh.e3[:, 1] * U[1] + self._mesh.e3[:, 2] * U[2]

        self.transpiration = zeros(self._n)

        t2 = time()

        print('    -> {} - data created'.format(t2 - t1))

        self.matrix, self.array, self.velXMatrix, self.velYMatrix, self.velZMatrix, self.velArray = create(
            vertices=self._mesh.vertices,
            faces=self._mesh.faces,
            facesAreas=self._mesh.facesAreas,
            facesMaxDistance=self._mesh.facesMaxDistance,
            facesCenter=self._mesh.facesCenter,
            controlPoints=self._mesh.controlPoints,
            e1=self._mesh.e1,
            e2=self._mesh.e2,
            e3=self._mesh.e3,
            freestream=U,
            sigma=self.sigma,
        )

        t3 = time()

        print('    -> {} - linear system'.format(t3 - t2))

        self._potentialFlow(self.transpiration)

        return
    
    def _potentialFlow(self, transpiration: ndarray) -> None:

        t4 = time()

        self.doublet = solve(self.matrix, self.array + transpiration)

        t5 = time()

        print('    -> {} - solve'.format(t5 - t4))

        for i in range(self._n):
            self.velField[i, 0] = sum(self.velXMatrix[i, :] * self.doublet) + self.velArray[i, 0]
            self.velField[i, 1] = sum(self.velYMatrix[i, :] * self.doublet) + self.velArray[i, 1]
            self.velField[i, 2] = sum(self.velZMatrix[i, :] * self.doublet) + self.velArray[i, 2]

            self.velNorm[i] = (self.velField[i, 0] ** 2 + self.velField[i, 1] ** 2 + self.velField[i, 2] ** 2) ** 0.5

            self.cp[i] = 1 - (self.velNorm[i] ** 2) / self._UNorm2

        t6 = time()

        print('    -> {} - data calculated'.format(t6 - t5))

        return
    