from math import sqrt
import numpy as np
import ctypes

ND_POINTER_DOUBLE = np.ctypeslib.ndpointer(dtype=np.double, ndim=1, flags="C")

#---------------------------------------------#
#                    Input                    #
#---------------------------------------------#
class SURFACE_MESH(ctypes.Structure):
    _fields_ = [
        ("nv", ctypes.c_int),
        ("nf", ctypes.c_int),
        ("vertices", ctypes.POINTER(ctypes.c_double)),
        ("faces", ctypes.POINTER(ctypes.c_int)),
        ("facesAreas", ctypes.POINTER(ctypes.c_double)),
        ("facesMaxDistance", ctypes.POINTER(ctypes.c_double)),
        ("facesCenter", ctypes.POINTER(ctypes.c_double)),
        ("controlPoints", ctypes.POINTER(ctypes.c_double)),
        ("p1", ctypes.POINTER(ctypes.c_double)),
        ("p2", ctypes.POINTER(ctypes.c_double)),
        ("p3", ctypes.POINTER(ctypes.c_double)),
        ("e1", ctypes.POINTER(ctypes.c_double)),
        ("e2", ctypes.POINTER(ctypes.c_double)),
        ("e3", ctypes.POINTER(ctypes.c_double)),
    ]

class WAKE_MESH_PART(ctypes.Structure):
    _fields_ = [
        ("nSpan", ctypes.c_int),
        ("nWake", ctypes.c_int),
        ("grid", ctypes.POINTER(ctypes.c_int)),
        ("vertices", ctypes.POINTER(ctypes.c_double)),
        ("faces", ctypes.POINTER(ctypes.c_int)),
    ]

class WAKE_MESH(ctypes.Structure):
    _fields_ = [
        ("left", WAKE_MESH_PART),
        ("right", WAKE_MESH_PART),
        ("tail", WAKE_MESH_PART),
    ]

class MESH(ctypes.Structure):
    _fields_ = [
        ("surface", SURFACE_MESH),
        ("wake", WAKE_MESH),
    ]

class ENVIRONMENT(ctypes.Structure):
    _fields_ = [
        ("vel_x", ctypes.c_double),
        ("vel_y", ctypes.c_double),
        ("vel_z", ctypes.c_double),
        ("velNorm", ctypes.c_double),
        ("density", ctypes.c_double),
        ("viscosity", ctypes.c_double),
        ("soundSpeed", ctypes.c_double),
    ]

class INPUT(ctypes.Structure):
    _fields_ = [
        ("type", ctypes.c_int),
        ("mesh", MESH),
        ("environment", ENVIRONMENT),
    ]

#---------------------------------------------#
#                   WRAPPER                   #
#---------------------------------------------#
def wrapper(vertices: np.ndarray,
            faces: np.ndarray,
            facesAreas: np.ndarray,
            facesMaxDistance: np.ndarray,
            facesCenter: np.ndarray,
            controlPoints: np.ndarray,
            gridWakeLeft: np.ndarray,
            verticesWakeLeft: np.ndarray,
            facesWakeLeft: np.ndarray,
            gridWakeRight: np.ndarray,
            verticesWakeRight: np.ndarray,
            facesWakeRight: np.ndarray,
            gridWakeTail: np.ndarray,
            verticesWakeTail: np.ndarray,
            facesWakeTail: np.ndarray,
            p1: np.ndarray, p2: np.ndarray, p3: np.ndarray,
            e1: np.ndarray, e2: np.ndarray, e3: np.ndarray,
            freestream: np.ndarray,
            density: float,
            viscosity: float,
            soundSpeed: float):

    nv = vertices.shape[0]
    nf = faces.shape[0]

    # Input/output
    surfaceMesh = SURFACE_MESH(
        nv,
        nf,
        np.ctypeslib.as_ctypes(vertices.astype(np.double).reshape(vertices.size)),
        np.ctypeslib.as_ctypes(faces.astype(np.int32).reshape(faces.size)),
        np.ctypeslib.as_ctypes(facesAreas.astype(np.double).reshape(facesAreas.size)),
        np.ctypeslib.as_ctypes(facesMaxDistance.astype(np.double).reshape(facesMaxDistance.size)),
        np.ctypeslib.as_ctypes(facesCenter.astype(np.double).reshape(facesCenter.size)),
        np.ctypeslib.as_ctypes(controlPoints.astype(np.double).reshape(controlPoints.size)),
        np.ctypeslib.as_ctypes(p1.astype(np.double).reshape(p1.size)),
        np.ctypeslib.as_ctypes(p2.astype(np.double).reshape(p2.size)),
        np.ctypeslib.as_ctypes(p3.astype(np.double).reshape(p3.size)),
        np.ctypeslib.as_ctypes(e1.astype(np.double).reshape(e1.size)),
        np.ctypeslib.as_ctypes(e2.astype(np.double).reshape(e2.size)),
        np.ctypeslib.as_ctypes(e3.astype(np.double).reshape(e3.size))
    )

    left = WAKE_MESH_PART(
        gridWakeLeft.shape[0],
        gridWakeLeft.shape[1],
        np.ctypeslib.as_ctypes(gridWakeLeft.astype(np.int32).reshape(gridWakeLeft.size)),
        np.ctypeslib.as_ctypes(verticesWakeLeft.astype(np.double).reshape(verticesWakeLeft.size)),
        np.ctypeslib.as_ctypes(facesWakeLeft.astype(np.int32).reshape(facesWakeLeft.size))
    )

    right = WAKE_MESH_PART(
        gridWakeRight.shape[0],
        gridWakeRight.shape[1],
        np.ctypeslib.as_ctypes(gridWakeRight.astype(np.int32).reshape(gridWakeRight.size)),
        np.ctypeslib.as_ctypes(verticesWakeRight.astype(np.double).reshape(verticesWakeRight.size)),
        np.ctypeslib.as_ctypes(facesWakeRight.astype(np.int32).reshape(facesWakeRight.size))
    )

    tail = WAKE_MESH_PART(
        gridWakeTail.shape[0],
        gridWakeTail.shape[1],
        np.ctypeslib.as_ctypes(gridWakeTail.astype(np.int32).reshape(gridWakeTail.size)),
        np.ctypeslib.as_ctypes(verticesWakeTail.astype(np.double).reshape(verticesWakeTail.size)),
        np.ctypeslib.as_ctypes(facesWakeTail.astype(np.int32).reshape(facesWakeTail.size))
    )

    wakeMesh = WAKE_MESH(
        left,
        right,
        tail
    )

    mesh = MESH(
        surfaceMesh,
        wakeMesh
    )

    environment = ENVIRONMENT(
        freestream[0],
        freestream[1],
        freestream[2],
        sqrt(freestream[0] * freestream[0] + freestream[1] * freestream[1] + freestream[2] * freestream[2]),
        density,
        viscosity,
        soundSpeed
    )

    input = INPUT(
        1,
        mesh,
        environment
    )

    cp_v = np.empty(nv, dtype=np.double)
    vel_x_v = np.empty(nv, dtype=np.double)
    vel_y_v = np.empty(nv, dtype=np.double)
    vel_z_v = np.empty(nv, dtype=np.double)
    transpiration_v = np.empty(nv, dtype=np.double)
    sigma_v = np.empty(nv, dtype=np.double)
    doublet_v = np.empty(nv, dtype=np.double)
    forces = np.empty(3, dtype=np.double)

    # Load library
    lib = ctypes.CDLL('./utils/bin/libsolver.so')

    # Set input and output
    lib.solve.argtypes = [
        INPUT,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
        ND_POINTER_DOUBLE,
    ]

    lib.solve.restype = None

    lib.solve(input, vel_x_v, vel_y_v, vel_z_v, transpiration_v, sigma_v, doublet_v, forces)

    print('Forces: {}'.format(forces))

    cp_v = 1 - (vel_x_v * vel_x_v + vel_y_v * vel_y_v + vel_z_v * vel_z_v) / (freestream[0] ** 2 + freestream[1] ** 2 + freestream[2] ** 2)

    return [
        cp_v,
        vel_x_v,
        vel_y_v,
        vel_z_v,
        transpiration_v,
        sigma_v,
        doublet_v,
    ]