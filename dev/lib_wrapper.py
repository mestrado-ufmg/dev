from typing import List, NamedTuple
from ctypes import CDLL, c_double, c_int
from numpy import ndarray, ctypeslib, empty, double, int32, zeros
from numpy.linalg import norm

ND_POINTER_DOUBLE = ctypeslib.ndpointer(dtype=double, ndim=1, flags="C")
ND_POINTER_INT = ctypeslib.ndpointer(dtype=int32, ndim=1, flags="C")

class FacesParams(NamedTuple):
    doublet: ndarray
    sigma: ndarray
    cp: ndarray
    vel_field: ndarray
    vel_norm: ndarray
    transpiration: ndarray
    delta: ndarray
    A: ndarray
    B: ndarray
    Psi: ndarray
    Ctau1: ndarray
    Ctau2: ndarray
    tau_wall: ndarray

class VerticesParams(NamedTuple):
    doublet: ndarray
    sigma: ndarray
    cp: ndarray
    vel_field: ndarray
    vel_norm: ndarray
    transpiration: ndarray
    delta: ndarray
    A: ndarray
    B: ndarray
    Psi: ndarray
    Ctau1: ndarray
    Ctau2: ndarray
    tau_wall: ndarray

def unary(a: ndarray):

    if a.ndim == 2:
        n = a.shape[0]
        for i in range(n):
            a[i, :] = a[i, :] / norm(a[i, :])
        return a

    return a / norm(a)

def solve(type: int,
          vertices: ndarray,
          faces: ndarray,
          facesAreas: ndarray,
          facesMaxDistance: ndarray,
          facesCenter: ndarray,
          controlPoints: ndarray,
          p1: ndarray, p2: ndarray, p3: ndarray,
          e1: ndarray, e2: ndarray, e3: ndarray,
          freestream: ndarray,
          leftWingGrid: ndarray,
          leftWingVertices: ndarray,
          leftWingFaces: ndarray,
          rightWingGrid: ndarray,
          rightWingVertices: ndarray,
          rightWingFaces: ndarray,
          tailGrid: ndarray,
          tailVertices: ndarray,
          tailFaces: ndarray,
          delta: ndarray, A: ndarray, B: ndarray, Psi: ndarray, Ctau1: ndarray, Ctau2: ndarray) -> List[ndarray]:
    
    # Calculated parameters
    nf = facesAreas.size
    nv = vertices.shape[0]
    sigma = -(e3[:, 0] * freestream[0] + e3[:, 1] * freestream[1] + e3[:, 2] * freestream[2])
    nSpanLeftWing = leftWingGrid.shape[0]; nWakeLeftWing = leftWingGrid.shape[1]
    nSpanRightWing = rightWingGrid.shape[0]; nWakeRightWing = rightWingGrid.shape[1] 
    nSpanTail = tailGrid.shape[0]; nWakeTail = tailGrid.shape[1]

    # Reshape
    vertices = vertices.reshape(vertices.size)
    faces = faces.reshape(faces.size)
    facesCenter = facesCenter.reshape(facesCenter.size)
    controlPoints = controlPoints.reshape(controlPoints.size)
    p1 = p1.reshape(p1.size)
    p2 = p2.reshape(p2.size)
    p3 = p3.reshape(p3.size)
    e1 = e1.reshape(e1.size)
    e2 = e2.reshape(e2.size)
    e3 = e3.reshape(e3.size)
    leftWingGrid = leftWingGrid.reshape(leftWingGrid.size)
    leftWingVertices = leftWingVertices.reshape(leftWingVertices.size)
    leftWingFaces = leftWingFaces.reshape(leftWingFaces.size)
    rightWingGrid = rightWingGrid.reshape(rightWingGrid.size)
    rightWingVertices = rightWingVertices.reshape(rightWingVertices.size)
    rightWingFaces = rightWingFaces.reshape(rightWingFaces.size)
    tailGrid = tailGrid.reshape(tailGrid.size)
    tailVertices = tailVertices.reshape(tailVertices.size)
    tailFaces = tailFaces.reshape(tailFaces.size)

    # Output

    # Faces values
    doublet = empty(nf, dtype=double)
    cp = empty(nf, dtype=double)
    velx = empty(nf, dtype=double)
    vely = empty(nf, dtype=double)
    velz = empty(nf, dtype=double)
    velNorm = empty(nf, dtype=double)
    mach = empty(nf, dtype=double)
    transpiration = zeros(nf, dtype=double)

    # delta = empty(nf, dtype=double)
    # A = empty(nf, dtype=double)
    # B = empty(nf, dtype=double)
    # Psi = empty(nf, dtype=double)
    # Ctau1 = empty(nf, dtype=double)
    # Ctau2 = empty(nf, dtype=double)
    tau_x = empty(nf, dtype=double)
    tau_y = empty(nf, dtype=double)
    tau_z = empty(nf, dtype=double)

    # Vertices values
    sigma_v = empty(nv, dtype=double)
    doublet_v = empty(nv, dtype=double)
    cp_v = empty(nv, dtype=double)
    velx_v = empty(nv, dtype=double)
    vely_v = empty(nv, dtype=double)
    velz_v = empty(nv, dtype=double)
    velNorm_v = empty(nv, dtype=double)
    transpiration_v = empty(nv, dtype=double)

    delta_v = empty(nv, dtype=double)
    A_v = empty(nv, dtype=double)
    B_v = empty(nv, dtype=double)
    Psi_v = empty(nv, dtype=double)
    Ctau1_v = empty(nv, dtype=double)
    Ctau2_v = empty(nv, dtype=double)
    tau_x_v = empty(nv, dtype=double)
    tau_y_v = empty(nv, dtype=double)
    tau_z_v = empty(nv, dtype=double)

    # Load library
    lib = CDLL('./bin/libsolver.so')

    # Set input and output
    lib.solve.argtypes = [
        c_int,
        c_int,
        c_int,
        ND_POINTER_DOUBLE,
        ND_POINTER_INT,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        c_int,
        c_int,
        ND_POINTER_INT,
        ND_POINTER_DOUBLE,
        ND_POINTER_INT,
        c_int,
        c_int,
        ND_POINTER_INT,
        ND_POINTER_DOUBLE,
        ND_POINTER_INT,
        c_int,
        c_int,
        ND_POINTER_INT,
        ND_POINTER_DOUBLE,
        ND_POINTER_INT,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        c_double,
        c_double,
        c_double,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE
    ]

    lib.solve.restype = None

    # Solve linear system
    lib.solve(type,
              nv,
              nf,
              vertices,
              faces,
              facesAreas,
              facesMaxDistance,
              facesCenter,
              controlPoints,
              p1, p2, p3,
              e1, e2, e3,
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
              doublet,
              velx, vely, velz,
              velNorm,
              cp,
              mach,
              delta,
              A, B,
              Psi,
              Ctau1, Ctau2,
              tau_x, tau_y, tau_z,
              1.225, 1e-5, 340,
              transpiration,
              sigma_v,
              doublet_v,
              cp_v,
              velx_v, vely_v, velz_v, velNorm_v,
              transpiration_v,
              delta_v,
              A_v, B_v,
              Psi_v,
              Ctau1_v, Ctau2_v,
              tau_x_v, tau_y_v, tau_z_v)
    
    # Create output

    velField = empty((nf, 3), dtype=double)
    velField[:, 0], velField[:, 1], velField[:, 2] = velx, vely, velz

    tauWall = empty((nf, 3), dtype=double)
    tauWall[:, 0], tauWall[:, 1], tauWall[:, 2] = tau_x, tau_y, tau_z

    facesParams = FacesParams(doublet, sigma, cp, velField, velNorm, transpiration, delta, A, B, Psi, Ctau1, Ctau2, tauWall)

    velField_v = empty((nv, 3), dtype=double)
    velField_v[:, 0], velField_v[:, 1], velField_v[:, 2] = velx_v, vely_v, velz_v

    tauWall_v = empty((nv, 3), dtype=double)
    tauWall_v[:, 0], tauWall_v[:, 1], tauWall_v[:, 2] = tau_x_v, tau_y_v, tau_z_v

    verticesParams = VerticesParams(doublet_v, sigma_v, cp_v, velField_v, velNorm_v, transpiration_v, delta_v, A_v, B_v, Psi_v, Ctau1_v, Ctau2_v, tauWall_v)
    
    return [facesParams, verticesParams]