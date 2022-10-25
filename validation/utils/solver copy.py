import numpy as np
from typing import NamedTuple
import ctypes
from scipy.spatial.transform import Rotation as R

from utils.mesh import MeshData

ND_POINTER_DOUBLE = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags="C")
ND_POINTER_INT = np.ctypeslib.ndpointer(dtype=np.int32, ndim=1, flags="C")

class SolverData(NamedTuple):
    cp_f: np.ndarray
    velx_f: np.ndarray
    vely_f: np.ndarray
    velz_f: np.ndarray
    transpiration_f: np.ndarray
    cp_v: np.ndarray
    velx_v: np.ndarray
    vely_v: np.ndarray
    velz_v: np.ndarray
    transpiration_v: np.ndarray

def solve(mesh: MeshData, alpha: float, freestream: float) -> SolverData:

    # Freestream vector
    x = np.array([-1.0, 0.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])

    if alpha is not None:
        r = R.from_rotvec(-np.deg2rad(alpha) * z)
        x = r.apply(x)
    
    freestream = freestream * x

    # Calculated parameters
    nf = mesh.facesAreas.size
    nv = mesh.vertices.shape[0]
    sigma = -(mesh.e3[:, 0] * freestream[0] + mesh.e3[:, 1] * freestream[1] + mesh.e3[:, 2] * freestream[2])
    nSpanWake = mesh.wake_grid.shape[0]
    nWake = mesh.wake_grid.shape[1] 

    # Type
    vertices = mesh.vertices.astype(np.double)
    faces = mesh.faces.astype(np.int32)
    facesCenter = mesh.facesCenter.astype(np.double)
    facesAreas = mesh.facesAreas.astype(np.double)
    facesMaxDistance = mesh.facesMaxDistance.astype(np.double)
    controlPoints = mesh.controlPoints.astype(np.double)
    p1 = mesh.p1Local.astype(np.double)
    p2 = mesh.p2Local.astype(np.double)
    p3 = mesh.p3Local.astype(np.double)
    e1 = mesh.e1.astype(np.double)
    e2 = mesh.e2.astype(np.double)
    e3 = mesh.e3.astype(np.double)
    wakeGrid = mesh.wake_grid.astype(np.int32)
    wakeVertices = mesh.wake_vertices.astype(np.double)
    wakeFaces = mesh.wake_faces.astype(np.int32)

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
    wakeGrid = wakeGrid.reshape(wakeGrid.size)
    wakeVertices = wakeVertices.reshape(wakeVertices.size)
    wakeFaces = wakeFaces.reshape(wakeFaces.size)

    # Faces values
    doublet = np.empty(nf, dtype=np.double)
    cp = np.empty(nf, dtype=np.double)
    velx = np.empty(nf, dtype=np.double)
    vely = np.empty(nf, dtype=np.double)
    velz = np.empty(nf, dtype=np.double)
    velNorm = np.empty(nf, dtype=np.double)
    mach = np.empty(nf, dtype=np.double)
    transpiration = np.zeros(nf, dtype=np.double)

    delta = np.empty(nf, dtype=np.double)
    A = np.empty(nf, dtype=np.double)
    B = np.empty(nf, dtype=np.double)
    Psi = np.empty(nf, dtype=np.double)
    Ctau1 = np.empty(nf, dtype=np.double)
    Ctau2 = np.empty(nf, dtype=np.double)
    tau_x = np.empty(nf, dtype=np.double)
    tau_y = np.empty(nf, dtype=np.double)
    tau_z = np.empty(nf, dtype=np.double)

    # Vertices values
    sigma_v = np.empty(nv, dtype=np.double)
    doublet_v = np.empty(nv, dtype=np.double)
    cp_v = np.empty(nv, dtype=np.double)
    velx_v = np.empty(nv, dtype=np.double)
    vely_v = np.empty(nv, dtype=np.double)
    velz_v = np.empty(nv, dtype=np.double)
    velNorm_v = np.empty(nv, dtype=np.double)
    transpiration_v = np.empty(nv, dtype=np.double)

    delta_v = np.empty(nv, dtype=np.double)
    A_v = np.empty(nv, dtype=np.double)
    B_v = np.empty(nv, dtype=np.double)
    Psi_v = np.empty(nv, dtype=np.double)
    Ctau1_v = np.empty(nv, dtype=np.double)
    Ctau2_v = np.empty(nv, dtype=np.double)
    tau_x_v = np.empty(nv, dtype=np.double)
    tau_y_v = np.empty(nv, dtype=np.double)
    tau_z_v = np.empty(nv, dtype=np.double)

    # Load library
    lib = ctypes.CDLL('./utils/bin/libsolver.so')

    # Set input and output
    lib.solve.argtypes = [
        ctypes.c_int,
        ctypes.c_int,
        ctypes.c_int,
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
        ctypes.c_int,
        ctypes.c_int,
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
        ctypes.c_double,
        ctypes.c_double,
        ctypes.c_double,
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
    lib.solve(1,
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
              nSpanWake,
              nWake,
              wakeGrid,
              wakeVertices,
              wakeFaces,
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
    
    data = SolverData(
        cp_f=cp,
        velx_f=velx,
        vely_f=vely,
        velz_f=velz,
        transpiration_f=transpiration,
        cp_v=cp_v,
        velx_v=velx_v,
        vely_v=vely_v,
        velz_v=velz_v,
        transpiration_v=transpiration_v,
    )

    return data