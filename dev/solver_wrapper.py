from ctypes import CDLL, c_int
from typing import List
from numpy import dot, empty, int32, ndarray, copy, ctypeslib, float64
from numpy.linalg import solve as solvesys 

ND_POINTER_1 = ctypeslib.ndpointer(dtype=float64, ndim=1, flags="C")
ND_POINTER_2 = ctypeslib.ndpointer(dtype=int32, ndim=1, flags="C")

def solve(facesAreas: ndarray,
          facesMaxDistance: ndarray,
          facesCenter: ndarray,
          controlPoints: ndarray,
          p1: ndarray, p2: ndarray, p3: ndarray,
          e1: ndarray, e2: ndarray, e3: ndarray,
          freestream: ndarray,
          sigma: ndarray,
          leftWingGrid: ndarray, leftWingVertices: ndarray, leftWingFaces: ndarray,
          rightWingGrid: ndarray, rightWingVertices: ndarray, rightWingFaces: ndarray,
          tailGrid: ndarray, tailVertices: ndarray, tailFaces: ndarray) -> List[ndarray]:

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
        c_int,
        c_int,
        ND_POINTER_2,
        ND_POINTER_1,
        ND_POINTER_2,
        c_int,
        c_int,
        ND_POINTER_2,
        ND_POINTER_1,
        ND_POINTER_2,
        c_int,
        c_int,
        ND_POINTER_2,
        ND_POINTER_1,
        ND_POINTER_2,
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

    nSpanLeftWing = leftWingGrid.shape[0]
    nWakeLeftWing = leftWingGrid.shape[1]

    nSpanRightWing = rightWingGrid.shape[0]
    nWakeRightWing = rightWingGrid.shape[1]

    nSpanTail = tailGrid.shape[0]
    nWakeTail = tailGrid.shape[1]

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
    freestream = copy(freestream.astype(float64))
    sigma = copy(sigma.astype(float64))
    leftWingGrid = copy(leftWingGrid.reshape(leftWingGrid.size).astype(int32))
    leftWingVertices = copy(leftWingVertices.reshape(leftWingVertices.size).astype(float64))
    leftWingFaces = copy(leftWingFaces.reshape(leftWingFaces.size).astype(int32))
    rightWingGrid = copy(rightWingGrid.reshape(rightWingGrid.size).astype(int32))
    rightWingVertices = copy(rightWingVertices.reshape(rightWingVertices.size).astype(float64))
    rightWingFaces = copy(rightWingFaces.reshape(rightWingFaces.size).astype(int32))
    tailGrid = copy(tailGrid.reshape(tailGrid.size).astype(int32))
    tailVertices = copy(tailVertices.reshape(tailVertices.size).astype(float64))
    tailFaces = copy(tailFaces.reshape(tailFaces.size).astype(int32))

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
              nSpanLeftWing,
              nWakeLeftWing,
              leftWingGrid,
              leftWingVertices,
              leftWingFaces,
              nSpanRightWing,
              nWakeRightWing,
              rightWingGrid,
              rightWingVertices,
              rightWingFaces,
              nSpanTail,
              nWakeTail,
              tailGrid,
              tailVertices,
              tailFaces,
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