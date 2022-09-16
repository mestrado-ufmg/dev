import sys

import numpy as np
sys.path.append('../')

import pybird

if __name__ == '__main__':
    
    # Create model
    model = pybird.model('test')

    # Load data
    model.load('./data/data.case')

    # Create geometry
    model.geo.build()

    # Create mesh
    model.mesh.build(
        size=0.01,
        accom_dist=10,
        n_head=2,
        n_wing_le=2,
        n_wing_te=2,
        n_tail_le=2,
        n_tail_te=2,
        n_body=2,
    )

    model.view.paraview()

    path = './data/mesh3/'

    # Save
    np.savetxt(path + 'faces.txt', model.mesh.faces)
    np.savetxt(path + 'vertices.txt', model.mesh.vertices)

    np.savetxt(path + 'facesAreas.txt', model.mesh.facesAreas)
    np.savetxt(path + 'facesMaxDistance.txt', model.mesh.facesMaxDistance)
    np.savetxt(path + 'facesCenter.txt', model.mesh.facesCenter)
    np.savetxt(path + 'controlPoints.txt', model.mesh.controlPoints)

    np.savetxt(path + 'leftWingGrid.txt', model.mesh.wake.leftWing.grid)
    np.savetxt(path + 'leftWingVertices.txt', model.mesh.wake.leftWing.vertices)
    np.savetxt(path + 'leftWingFaces.txt', model.mesh.wake.leftWing.faces)

    np.savetxt(path + 'rightWingGrid.txt', model.mesh.wake.rightWing.grid)
    np.savetxt(path + 'rightWingVertices.txt', model.mesh.wake.rightWing.vertices)
    np.savetxt(path + 'rightWingFaces.txt', model.mesh.wake.rightWing.faces)

    np.savetxt(path + 'tailGrid.txt', model.mesh.wake.tail.grid)
    np.savetxt(path + 'tailVertices.txt', model.mesh.wake.tail.vertices)
    np.savetxt(path + 'tailFaces.txt', model.mesh.wake.tail.faces)

    np.savetxt(path + 'p1.txt', model.mesh.p1Local)
    np.savetxt(path + 'p2.txt', model.mesh.p2Local)
    np.savetxt(path + 'p3.txt', model.mesh.p3Local)

    np.savetxt(path + 'e1.txt', model.mesh.e1)
    np.savetxt(path + 'e2.txt', model.mesh.e2)
    np.savetxt(path + 'e3.txt', model.mesh.e3)

    freestream = np.array([-10, 0, 0])

    np.savetxt(path + 'freestream.txt', freestream)