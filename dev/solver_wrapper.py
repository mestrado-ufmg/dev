from ctypes import CDLL, c_int
from typing import List
from numpy import dot, empty, ndarray, copy, ctypeslib, float64
from numpy.linalg import solve as solvesys 

ND_POINTER_1 = ctypeslib.ndpointer(dtype=float64, ndim=1, flags="C")

def solve(facesAreas: ndarray,
          facesMaxDistance: ndarray,
          facesCenter: ndarray,
          controlPoints: ndarray,
          p1: ndarray, p2: ndarray, p3: ndarray,
          e1: ndarray, e2: ndarray, e3: ndarray,
          freestream: ndarray,
          sigma: ndarray) -> List[ndarray]:

    # Load library
    lib = CDLL('./solver_lib.so')

    # Set input and output
    lib.solve.argtypes = [
        c_int,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
        ND_POINTER_1,
    ]

    lib.solve.restype = None
    
    # Change input shape
    n = facesAreas.shape[0]
    
    facesAreas = copy(facesAreas.reshape(facesAreas.size).astype(float64))
    facesMaxDistance = copy(facesMaxDistance.reshape(facesMaxDistance.size).astype(float64))
    facesCenter = copy(facesCenter.reshape(facesCenter.size).astype(float64))
    controlPoints = copy(controlPoints.reshape(controlPoints.size).astype(float64))
    p1 = copy(p1.reshape(p1.size).astype(float64))
    p2 = copy(p2.reshape(p2.size).astype(float64))
    p3 = copy(p3.reshape(p3.size).astype(float64))
    e1Vec = copy(e1.reshape(e1.size).astype(float64))
    e2Vec = copy(e2.reshape(e2.size).astype(float64))
    e3Vec = copy(e3.reshape(e3.size).astype(float64))
    freestream = copy(freestream)
    sigma = copy(sigma)

    # Create output
    matrix = empty(n * n, dtype=float64)
    array = empty(n, dtype=float64)
    matrixVelx = empty(n * n, dtype=float64)
    matrixVely = empty(n * n, dtype=float64)
    matrixVelz = empty(n * n, dtype=float64)
    arrayVel = empty(n * 3, dtype=float64)

    # Call library function
    lib.solve(n,
              facesAreas,
              facesMaxDistance,
              facesCenter,
              controlPoints,
              p1, p2, p3,
              e1Vec, e2Vec, e3Vec,
              freestream,
              sigma,
              matrix,
              array,
              matrixVelx,
              matrixVely,
              matrixVelz,
              arrayVel)
    
    matrix = matrix.reshape((n, n))

    # Doublet
    doublet = solvesys(matrix, array)

    # Velocity
    matrixVelx = matrixVelx.reshape((n, n))
    matrixVely = matrixVely.reshape((n, n))
    matrixVelz = matrixVelz.reshape((n, n))
    arrayVel = arrayVel.reshape((n, 3))

    velxDoublet = dot(matrixVelx, doublet)
    velyDoublet = dot(matrixVely, doublet)
    velzDoublet = dot(matrixVelz, doublet)

    velx = velxDoublet + arrayVel[:, 0]
    vely = velyDoublet + arrayVel[:, 1]
    velz = velzDoublet + arrayVel[:, 2]

    velAux = velx * velx + vely * vely + velz * velz
    velNorm = (velAux) ** 0.5

    # Cp
    freestreamSquare = freestream[0] ** 2 + freestream[1] ** 2 + freestream[2] ** 2
    cp = 1 - velAux / freestreamSquare

    # Transpiration
    transpiration = velx * e3[:, 0] + vely * e3[:, 1] + velz * e3[:, 2]

    return [doublet, cp, velNorm, velx, vely, velz, transpiration]