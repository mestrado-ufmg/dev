from typing import NamedTuple
import numpy as np
from scipy.spatial.transform import Rotation as R
from utils.bin.wrapper import wrapper
from utils.mesh import MeshData

class SolverData(NamedTuple):
    cp_v: np.ndarray
    velx_v: np.ndarray
    vely_v: np.ndarray
    velz_v: np.ndarray
    sigma_v: np.ndarray
    doublet_v: np.ndarray
    transpiration_v: np.ndarray

def solve(mesh: MeshData,
          freestream: float,
          alpha: float,
          density: float,
          viscosity: float,
          soundSpeed: float) -> None:

    # Freestream vector
    x = np.array([-1.0, 0.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])

    if alpha is not None:
        r = R.from_rotvec(-np.deg2rad(alpha) * z)
        x = r.apply(x)
    
    freestream = freestream * x

    sol = wrapper(
        mesh.vertices,
        mesh.faces,
        mesh.facesAreas,
        mesh.facesMaxDistance,
        mesh.facesCenter,
        mesh.controlPoints,
        np.zeros((2, 2), dtype=np.int32),
        np.zeros((2, 2), dtype=np.double),
        np.zeros((2, 2), dtype=np.int32),
        np.zeros((3, 3), dtype=np.int32),
        np.zeros((3, 3), dtype=np.double),
        np.zeros((3, 3), dtype=np.int32),
        mesh.wake_grid,
        mesh.wake_vertices,
        mesh.wake_faces,
        mesh.p1Local, mesh.p2Local, mesh.p3Local,
        mesh.e1, mesh.e2, mesh.e3,
        freestream, density, viscosity, soundSpeed)
    
    out = SolverData(
        cp_v=sol[0],
        velx_v=sol[1],
        vely_v=sol[2],
        velz_v=sol[3],
        transpiration_v=sol[4],
        sigma_v=sol[5],
        doublet_v=sol[6]
    )

    return out