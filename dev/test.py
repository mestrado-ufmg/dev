import numpy as np
from time import time

from lib_wrapper import solve
from gen_vtk2 import genVTK

if __name__ == '__main__':

    t1 = time()

    path = './data/mesh2/'

    # Load
    vertices = np.loadtxt(path + 'vertices.txt', dtype=np.double)
    faces = np.loadtxt(path + 'faces.txt', dtype=np.int32)
    facesAreas = np.loadtxt(path + 'facesAreas.txt', dtype=np.double)
    facesMaxDistance = np.loadtxt(path + 'facesMaxDistance.txt', dtype=np.double)
    facesCenter = np.loadtxt(path + 'facesCenter.txt', dtype=np.double)
    controlPoints = np.loadtxt(path + 'controlPoints.txt', dtype=np.double)
    p1 = np.loadtxt(path + 'p1.txt', dtype=np.double)
    p2 = np.loadtxt(path + 'p2.txt', dtype=np.double)
    p3 = np.loadtxt(path + 'p3.txt', dtype=np.double)
    e1 = np.loadtxt(path + 'e1.txt', dtype=np.double)
    e2 = np.loadtxt(path + 'e2.txt', dtype=np.double)
    e3 = np.loadtxt(path + 'e3.txt', dtype=np.double)
    freestream = np.loadtxt(path + 'freestream.txt', dtype=np.double)
    leftWingGrid = np.loadtxt(path + 'leftWingGrid.txt', dtype=np.int32)
    leftWingVertices = np.loadtxt(path + 'leftWingVertices.txt', dtype=np.double)
    leftWingFaces = np.loadtxt(path + 'leftWingFaces.txt', dtype=np.int32)
    rightWingGrid = np.loadtxt(path + 'rightWingGrid.txt', dtype=np.int32)
    rightWingVertices = np.loadtxt(path + 'rightWingVertices.txt', dtype=np.double)
    rightWingFaces = np.loadtxt(path + 'rightWingFaces.txt', dtype=np.int32)
    tailGrid = np.loadtxt(path + 'tailGrid.txt', dtype=np.int32)
    tailVertices = np.loadtxt(path + 'tailVertices.txt', dtype=np.double)
    tailFaces = np.loadtxt(path + 'tailFaces.txt', dtype=np.int32)
    delta = np.loadtxt(path + 'delta.txt', dtype=np.double)
    A = np.loadtxt(path + 'A.txt', dtype=np.double)
    B = np.loadtxt(path + 'B.txt', dtype=np.double)
    Psi = np.loadtxt(path + 'Psi.txt', dtype=np.double)
    Ctau1 = np.loadtxt(path + 'Ctau1.txt', dtype=np.double)
    Ctau2 = np.loadtxt(path + 'Ctau2.txt', dtype=np.double)

    data = solve(1,
                 vertices,
                 faces,
                 facesAreas,
                 facesMaxDistance,
                 facesCenter,
                 controlPoints,
                 p1, p2, p3,
                 e1, e2, e3,
                 freestream,
                 leftWingGrid, leftWingVertices, leftWingFaces,
                 rightWingGrid, rightWingVertices, rightWingFaces,
                 tailGrid, tailVertices, tailFaces,
                 delta, A, B, Psi, Ctau1, Ctau2)
    
    facesParams, verticesParams = data
    
    np.savetxt(path + 'delta.txt', facesParams.delta)
    np.savetxt(path + 'A.txt', facesParams.A)
    np.savetxt(path + 'B.txt', facesParams.B)
    np.savetxt(path + 'Psi.txt', facesParams.Psi)
    np.savetxt(path + 'Ctau1.txt', facesParams.Ctau1)
    np.savetxt(path + 'Ctau2.txt', facesParams.Ctau2)
    
    genVTK('new5', './data/vtk/', faces, vertices, verticesParams)

    t2 = time()

    print("Execution time: {} s".format(t2 - t1))