from typing import List
import numpy as np

from singularities import source_velocity

def unary(a: np.ndarray):

    if a.ndim == 2:
        n = a.shape[0]
        for i in range(n):
            a[i, :] = a[i, :] / np.linalg.norm(a[i, :])
        return a

    return a / np.linalg.norm(a)

def solver(vertices: np.ndarray, faces: np.ndarray, freestream: np.ndarray) -> List[np.ndarray]:

    nv = vertices.shape[0]
    nf = faces.shape[0]

    # Faces center
    facesCenter = (1 / 3) * (vertices[faces[:, 0], :] + vertices[faces[:, 1], :] + vertices[faces[:, 2], :])

    # Base vectors
    e3 = unary(np.cross(vertices[faces[:, 1], :] - vertices[faces[:, 0], :], vertices[faces[:, 2], :] - vertices[faces[:, 0], :]))
    e1 = unary(vertices[faces[:, 1], :] - facesCenter)
    e2 = unary(np.cross(e3, e1))

    # Control points
    controlPoints = facesCenter + 1e-2 * e3

    # Linear system
    m = np.zeros((nf, nf))
    a = np.zeros(nf)

    for i in range(nf):
        for j in range(nf):
            
            p1_inertial = vertices[faces[j, 0], :] - facesCenter[j, :]
            p2_inertial = vertices[faces[j, 1], :] - facesCenter[j, :]
            p3_inertial = vertices[faces[j, 2], :] - facesCenter[j, :]

            p1 = np.array([np.dot(p1_inertial, e1[j, :]), np.dot(p1_inertial, e2[j, :]), np.dot(p1_inertial, e3[j, :])])
            p2 = np.array([np.dot(p2_inertial, e1[j, :]), np.dot(p2_inertial, e2[j, :]), np.dot(p2_inertial, e3[j, :])])
            p3 = np.array([np.dot(p3_inertial, e1[j, :]), np.dot(p3_inertial, e2[j, :]), np.dot(p3_inertial, e3[j, :])])

            p_inertial = controlPoints[i, :] - facesCenter[j, :]
            p = np.array([np.dot(p_inertial, e1[j, :]), np.dot(p_inertial, e2[j, :]), np.dot(p_inertial, e3[j, :])])

            v = source_velocity(
                sigma=1,
                p1=p1,
                p2=p2,
                p3=p3,
                p=p,
                e1=e1[j, :],
                e2=e2[j, :],
                e3=e3[j, :],
            )

            m[i, j] = np.dot(v, e3[i, :])
        
        a[i] = -np.dot(freestream, e3[i, :])

    # Solve
    sol = np.linalg.solve(m, a)

    # Calculate the velocity
    vel = np.zeros((nf, 3))

    for i in range(nf):
        for j in range(nf):
            
            p1_inertial = vertices[faces[j, 0], :] - facesCenter[j, :]
            p2_inertial = vertices[faces[j, 1], :] - facesCenter[j, :]
            p3_inertial = vertices[faces[j, 2], :] - facesCenter[j, :]

            p1 = np.array([np.dot(p1_inertial, e1[j, :]), np.dot(p1_inertial, e2[j, :]), np.dot(p1_inertial, e3[j, :])])
            p2 = np.array([np.dot(p2_inertial, e1[j, :]), np.dot(p2_inertial, e2[j, :]), np.dot(p2_inertial, e3[j, :])])
            p3 = np.array([np.dot(p3_inertial, e1[j, :]), np.dot(p3_inertial, e2[j, :]), np.dot(p3_inertial, e3[j, :])])

            p_inertial = controlPoints[i, :] - facesCenter[j, :]
            p = np.array([np.dot(p_inertial, e1[j, :]), np.dot(p_inertial, e2[j, :]), np.dot(p_inertial, e3[j, :])])

            v = source_velocity(
                sigma=sol[j],
                p1=p1,
                p2=p2,
                p3=p3,
                p=p,
                e1=e1[j, :],
                e2=e2[j, :],
                e3=e3[j, :],
            )

            vel[i, :] = vel[i, :] + v
        
        vel[i, :] = vel[i, :] + freestream

    return facesCenter, controlPoints, e1, e2, e3, vel