from typing import List
from numpy import arctan, array, deg2rad, ndarray, pi, zeros, zeros_like, log
from scipy.spatial.transform import Rotation as R
from tqdm import tqdm

from pybird.mesh.mesh import Mesh
from pybird.solver.utils.singularities import source_velocity
from pybird.solver.utils.solvers import linear_system
from pybird.solver.utils import velocities, solution

ZERO_ERROR = 1e-8
def division(a: ndarray, b: ndarray) -> float:
    if isinstance(a, float or int):
        if (-ZERO_ERROR < b) and (b < 0):
            b = -ZERO_ERROR
        elif (0 < b) and (b < ZERO_ERROR):
            b = ZERO_ERROR
        return a / b
    b[(-ZERO_ERROR < b) & (b < 0)] = -ZERO_ERROR
    b[(0 < b) & (b < ZERO_ERROR)] = ZERO_ERROR
    return a / b

class Solver:

    def __init__(self, mesh: Mesh) -> None:
        self._mesh = mesh
        return
    
    def solve(self, freestream: float, alpha: float = None, beta: float = None) -> None:
        
        print('-Solver')

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
        
        U = x * freestream

        # self.sigma = solution.source_distribution(
        #     U,
        #     zeros(len(self._mesh.faces[:, 0])),
        #     self._mesh.controlPointsDistance,
        #     self._mesh.p1Local,
        #     self._mesh.p2Local,
        #     self._mesh.p3Local,
        #     self._mesh.e1,
        #     self._mesh.e2,
        #     self._mesh.e3
        # )

        self.sigma = self._source(U, zeros(len(self._mesh.faces[:, 0])))
        self.velField, self.velNorm = self._velocity(self.sigma, U)

        return
    
    def _velocity(self, sigma: ndarray, U: ndarray) -> List[ndarray]:

        print('  > Calculating velocity')

        n = len(self._mesh.faces[:, 0])

        a = zeros((n, 3))
        a_norm = zeros(n)

        for i in tqdm(range(n)):

            a[i, 0] = sum(self.velField[i, :, 0] * self.sigma[i]) + U[0]
            a[i, 1] = sum(self.velField[i, :, 1] * self.sigma[i]) + U[1]
            a[i, 2] = sum(self.velField[i, :, 2] * self.sigma[i]) + U[2]

            a_norm[i] = (a[i, 0] ** 2 + a[i, 1] ** 2 + a[i, 2] ** 2) ** 0.5

        return a, a_norm
    
    def _source(self, U: ndarray, transpiration: ndarray) -> ndarray:
        # sigma = -(self._mesh.e3[:, 0] * U[0] + self._mesh.e3[:, 1] * U[1] + self._mesh.e3[:, 2] * U[2])
        # return sigma

        print('  > Creating linear system [source]')
        n = len(self._mesh.faces[:, 0])

        a = zeros((n, n))
        b = zeros_like(transpiration)

        velField = zeros((n, n, 3))

        mult = (-1 / (4 * pi))
        vel = zeros((n, 3))
        p = zeros((n, 3))

        for i in tqdm(range(n)):
            p1 = self._mesh.p1Local[i, :]
            p2 = self._mesh.p2Local[i, :]
            p3 = self._mesh.p3Local[i, :]
            pxAux = self._mesh.controlPoints[:, 0] - self._mesh.facesCenters[i, 0]
            pyAux = self._mesh.controlPoints[:, 1] - self._mesh.facesCenters[i, 1]
            pzAux = self._mesh.controlPoints[:, 2] - self._mesh.facesCenters[i, 2]
            e1 = self._mesh.e1[i, :]
            e2 = self._mesh.e2[i, :]
            e3 = self._mesh.e3[i, :]
            p[:, 0] = pxAux * e1[0] + pyAux * e1[1] + pzAux * e1[2]
            p[:, 1] = pxAux * e2[0] + pyAux * e2[1] + pzAux * e2[2]
            p[:, 2] = pxAux * e3[0] + pyAux * e3[1] + pzAux * e3[2]

            b[i] = transpiration[i] - (U[0] * e3[0] + U[1] * e3[1] + U[2] * e3[2])

            pZSquare = p[:, 2] ** 2
            dpXSquare = (p[:, 0] - p1[0]) ** 2
            dpYSquare = (p[:, 0] - p2[0]) ** 2
            dpZSquare = (p[:, 0] - p3[0]) ** 2

            r1 = (dpXSquare + (p[:, 1] - p1[1]) ** 2 + pZSquare) ** 0.5
            r2 = (dpYSquare + (p[:, 1] - p2[1]) ** 2 + pZSquare) ** 0.5
            r3 = (dpZSquare + (p[:, 1] - p3[1]) ** 2 + pZSquare) ** 0.5
            l1 = dpXSquare + pZSquare
            l2 = dpYSquare + pZSquare
            l3 = dpZSquare + pZSquare
            h1 = (p[:, 0] - p1[0]) * (p[:, 1] - p1[1])
            h2 = (p[:, 0] - p2[0]) * (p[:, 1] - p2[1])
            h3 = (p[:, 0] - p3[0]) * (p[:, 1] - p3[1])

            d12 = ((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2) ** 0.5
            m12 = division((p2[1] - p1[1]), (p2[0] - p1[0]))
            d23 = ((p3[0] - p2[0]) ** 2 + (p3[1] - p2[1]) ** 2) ** 0.5
            m23 = division((p3[1] - p2[1]), (p3[0] - p2[0]))
            d31 = ((p1[0] - p3[0]) ** 2 + (p1[1] - p3[1]) ** 2) ** 0.5
            m31 = division((p1[1] - p3[1]), (p1[0] - p3[0]))

            ln12 = log( division((r1 + r2 - d12), (r1 + r2 + d12)) )
            ln23 = log( division((r2 + r3 - d23), (r2 + r3 + d23)) )
            ln31 = log( division((r3 + r1 - d31), (r3 + r1 + d31)) )

            u = -mult * ( division((p2[1] - p1[1]), d12) * ln12 + division((p3[1] - p2[1]), d23) * ln23 + division((p1[1] - p3[1]), d31) * ln31 )
            v = mult * ( division((p2[0] - p1[0]), d12) * ln12 + division((p3[0] - p2[0]), d23) * ln23 + division((p1[0] - p3[0]), d31) * ln31 )
            w = -mult * ( arctan(division(m12 * l1 - h1, p[:, 2] * r1)) - arctan(division(m12 * l2 - h2, p[:, 2] * r2)) + arctan(division(m23 * l2 - h2, p[:, 2] * r2)) - arctan(division(m23 * l3 - h3, p[:, 2] * r3)) + arctan(division(m31 * l3 - h3, p[:, 2] * r3)) - arctan(division(m31 * l1 - h1, p[:, 2] * r1)) )

            vel[:, 0] = u * e1[0] + v * e2[0] + w * e3[0]
            vel[:, 1] = u * e1[1] + v * e2[1] + w * e3[1]
            vel[:, 2] = u * e1[2] + v * e2[2] + w * e3[2]

            velField[i, :, 0] = vel[:, 0]
            velField[i, :, 1] = vel[:, 1]
            velField[i, :, 2] = vel[:, 2]

            a[i, :] = vel[:, 0] * e3[0] + vel[:, 1] * e3[1] + vel[:, 2] * e3[2]
        
        print('  > Solving linear system [source]')
        sol = linear_system(a, b)
        self.velField = velField

        return -(self._mesh.e3[:, 0] * U[0] + self._mesh.e3[:, 1] * U[1] + self._mesh.e3[:, 2] * U[2])
        # return sol