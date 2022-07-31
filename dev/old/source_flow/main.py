import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import matplotlib.pyplot as plt

from singularities import source_velocity
from solver import solver

def velFieldTest():

    freestream = np.array([1, 0, 0])

    n = 5
    x = np.linspace(0, n, num=n)

    faces = np.zeros((n * n, 3))
    vertices = np.zeros((n * n * 3, 3))

    countF = 0
    countV = 0
    for i in range(n):
        for j in range(n):
            
            faces[countF, :] = np.array([countV, countV + 1, countV + 2])
            countF += 1

            vertices[countV, :] = np.random.random(3) + np.array([x[i], x[j], 0])
            countV += 1
            vertices[countV, :] = np.random.random(3) + np.array([x[i], x[j], 0])
            countV += 1
            vertices[countV, :] = np.random.random(3) + np.array([x[i], x[j], 0])
            countV += 1
    
    faces = faces.astype(np.int16)
    vertices = vertices.astype(np.float64)

    facesCenter, controlPoints, e1, e2, e3, vel = solver(vertices, faces, freestream)

    for i in range(len(faces[:, 0])):
        print('Normal velocity: {}'.format(np.dot(vel[i, :], e3[i, :])))

    # View
    plotBase = True
    plotVel = True

    facesPlot = []

    for i in range(n * n):
        facesPlot.append([vertices[faces[i, 0], :], vertices[faces[i, 1], :], vertices[faces[i, 2], :]])
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    ax.add_collection3d(Poly3DCollection(facesPlot, linewidths=1))

    if plotBase:
        ax.quiver(facesCenter[:, 0], facesCenter[:, 1], facesCenter[:, 2], e1[:, 0], e1[:, 1], e1[:, 2], color='k')
        ax.quiver(facesCenter[:, 0], facesCenter[:, 1], facesCenter[:, 2], e2[:, 0], e2[:, 1], e2[:, 2], color='k')
        ax.quiver(facesCenter[:, 0], facesCenter[:, 1], facesCenter[:, 2], e3[:, 0], e3[:, 1], e3[:, 2], color='k')
    
    if plotVel:
        ax.quiver(facesCenter[:, 0], facesCenter[:, 1], facesCenter[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], color='r')

    ax.set_xlim(-0.5, n + 0.5)
    ax.set_ylim(-0.5, n + 0.5)
    ax.set_zlim(-(n + 1) / 2, (n + 1) / 2)

    plt.show()

    return

def sourceFlowField():

    n = 100
    x = np.linspace(-2, 2, num=n)
    y = np.linspace(-2, 2, num=n)
    z = -0.1

    X, Y = np.meshgrid(x, y)

    vel = np.zeros((n, n, 3))

    p1 = np.array([1, 0, 0])
    p2 = np.array([0, 1, 0])
    p3 = np.array([0, 0, 0])
    pCenter = (1 / 3) * (p1 + p2 + p3)

    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    for i in range(n):
        for j in range(n):

            v = source_velocity(
                sigma=1,
                p1=p1,
                p2=p2,
                p3=p3,
                p=np.array([X[i, j], Y[i, j], z]) - pCenter,
                e1=e1,
                e2=e2,
                e3=e3,
            )

            vel[i, j, :] = v
    
    plt.figure()
    plt.contourf(X, Y, vel[:, :, 2])
    plt.colorbar()
    plt.axis('equal')

    plt.figure()
    plt.quiver(X, Y, vel[:, :, 0], vel[:, :, 1])
    plt.axis('equal')

    plt.show()

    return

if __name__ == '__main__':
    velFieldTest()