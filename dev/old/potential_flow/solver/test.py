import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt
from matplotlib import cm

from system_wrapper_2 import create

if __name__ == '__main__':

    vertices = np.array([
        [0, -1.5, 0],
        [0, -0.5, 0],
        [0, 0.5, 0],
        [0, 1.5, 0],
        [-1, -1, 0.1],
        [-1, 0, 0.1],
        [-1, 1, 0.1],
        [-2, -0.5, 0.3],
        [-2, 0.5, 0.3],
        [-3, 0, 0.8],
    ]).astype(np.double)

    vertices[:, 0] = vertices[:, 0] - 0.5 * (max(vertices[:, 0]) + min(vertices[:, 0]))
    vertices[:, 1] = vertices[:, 1] - 0.5 * (max(vertices[:, 1]) + min(vertices[:, 1]))
    vertices[:, 2] = vertices[:, 2] - 0.5 * (max(vertices[:, 2]) + min(vertices[:, 2]))

    faces = np.array([
        [0, 1, 4],
        [1, 2, 5],
        [2, 3, 6],
        [1, 5, 4],
        [2, 6, 5],
        [4, 5, 7],
        [5, 6, 8],
        [5, 8, 7],
        [7, 8, 9],
    ]).astype(np.int32)

    n = len(faces[:, 0])

    auxVec = np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]], vertices[faces[:, 2]] - vertices[faces[:, 0]])
    auxNorm = (auxVec[:, 0] ** 2 + auxVec[:, 1] ** 2 + auxVec[:, 2] ** 2) ** 0.5
    facesAreas = 0.5 * auxNorm
    
    facesMaxDistance = 5 * (4 * facesAreas / np.pi) ** 0.5

    facesCenter = (1 / 3) * (vertices[faces[:, 0]] + vertices[faces[:, 1]] + vertices[faces[:, 2]])

    e3 = np.zeros((n, 3))
    e3[:, 0] = auxVec[:, 0] / auxNorm
    e3[:, 1] = auxVec[:, 1] / auxNorm
    e3[:, 2] = auxVec[:, 2] / auxNorm

    auxVec = vertices[faces[:, 1]] - vertices[faces[:, 0]]
    auxNorm = (auxVec[:, 0] ** 2 + auxVec[:, 1] ** 2 + auxVec[:, 2] ** 2) ** 0.5
    e1 = np.zeros((n, 3))
    e1[:, 0] = auxVec[:, 0] / auxNorm
    e1[:, 1] = auxVec[:, 1] / auxNorm
    e1[:, 2] = auxVec[:, 2] / auxNorm
    
    e2 = np.cross(e3, e1)

    controlPoints = facesCenter + 1e-5 * e3


    freestream = np.array([1, 0, 0])

    sigma = e3[:, 0] * freestream[0] + e3[:, 1] * freestream[1] + e3[:, 2] * freestream[2]
    
    transpiration = np.ones(len(faces[:, 0])) * (-0.5)
    
    matrix, array, velMatrix, velArray = create(
        vertices=vertices,
        faces=faces,
        facesAreas=facesAreas,
        facesMaxDistance=facesMaxDistance,
        facesCenter=facesCenter,
        controlPoints=controlPoints,
        e1=e1,
        e2=e2,
        e3=e3,
        freestream=freestream,
        sigma=sigma,
    )

    doublet = np.linalg.solve(matrix, array + transpiration)

    freestreamNorm = np.linalg.norm(freestream) ** 2

    vel = np.zeros((n, 3))
    velNorm = np.zeros(n)
    deltaPressure = np.zeros(n)
    for i in range(n):
        vel[i, 0] = np.sum(velMatrix[i, :, 0] * doublet) + velArray[i, 0]
        vel[i, 1] = np.sum(velMatrix[i, :, 1] * doublet) + velArray[i, 1]
        vel[i, 2] = np.sum(velMatrix[i, :, 2] * doublet) + velArray[i, 2]

        velNorm[i] = np.linalg.norm(vel[i, :])

        deltaPressure[i] = 0.5 * 1.225 * (freestreamNorm - velNorm[i] ** 2)
    
    print(vel[:, 0] * e3[:, 0] + vel[:, 1] * e3[:, 1] + vel[:, 2] * e3[:, 2])
    
    mapper = cm.ScalarMappable(cmap='Wistia')
    colors = [(r, g, b) for r, g, b, a in mapper.to_rgba(velNorm)]

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    poly3d = []

    for i in range(n):
        poly3d.append([vertices[faces[i, 0], :], vertices[faces[i, 1], :], vertices[faces[i, 2], :]])

    ax.add_collection3d(Poly3DCollection(poly3d, linewidths=1, facecolors=colors))

    if True:
        ax.quiver(controlPoints[:, 0], controlPoints[:, 1], controlPoints[:, 2], e1[:, 0], e1[:, 1], e1[:, 2], color='g')
        ax.quiver(controlPoints[:, 0], controlPoints[:, 1], controlPoints[:, 2], e2[:, 0], e2[:, 1], e2[:, 2], color='g')
        ax.quiver(controlPoints[:, 0], controlPoints[:, 1], controlPoints[:, 2], e3[:, 0], e3[:, 1], e3[:, 2], color='g')

        ax.quiver(controlPoints[:, 0], controlPoints[:, 1], controlPoints[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], color='r')

    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_zlim(-1.5, 1.5)

    plt.show()