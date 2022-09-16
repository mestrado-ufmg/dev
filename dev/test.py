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
                 tailGrid, tailVertices, tailFaces)
    
    facesParams, verticesParams = data
    
    genVTK('new', './data/vtk/', faces, vertices, verticesParams)

    t2 = time()

    print("Execution time: {} s".format(t2 - t1))