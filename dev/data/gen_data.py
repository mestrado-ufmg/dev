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
        n_head=3,
        n_wing_le=3,
        n_wing_te=2,
        n_tail_le=3,
        n_tail_te=2,
    )

    path = './data/mesh2/'

    # Save
    np.savetxt(path + 'faces.txt', model.mesh.faces)
    np.savetxt(path + 'vertices.txt', model.mesh.vertices)

    np.savetxt(path + 'facesAreas.txt', model.mesh.facesAreas)
    np.savetxt(path + 'facesMaxDistance.txt', model.mesh.facesMaxDistance)
    np.savetxt(path + 'facesCenter.txt', model.mesh.facesCenter)
    np.savetxt(path + 'controlPoints.txt', model.mesh.controlPoints)

    p1 = model.mesh.vertices[model.mesh.faces[:, 0], :] - model.mesh.facesCenter
    p2 = model.mesh.vertices[model.mesh.faces[:, 1], :] - model.mesh.facesCenter
    p3 = model.mesh.vertices[model.mesh.faces[:, 2], :] - model.mesh.facesCenter

    n = model.mesh.faces.shape[0]

    p1Local = np.empty((n, 2))
    p2Local = np.empty((n, 2))
    p3Local = np.empty((n, 2))

    p1Local[:, 0], p1Local[:, 1] = (p1 * model.mesh.e1).sum(axis=1), (p1 * model.mesh.e2).sum(axis=1)
    p2Local[:, 0], p2Local[:, 1] = (p2 * model.mesh.e1).sum(axis=1), (p2 * model.mesh.e2).sum(axis=1)
    p3Local[:, 0], p3Local[:, 1] = (p3 * model.mesh.e1).sum(axis=1), (p3 * model.mesh.e2).sum(axis=1)

    np.savetxt(path + 'p1.txt', p1Local)
    np.savetxt(path + 'p2.txt', p2Local)
    np.savetxt(path + 'p3.txt', p3Local)

    np.savetxt(path + 'e1.txt', model.mesh.e1)
    np.savetxt(path + 'e2.txt', model.mesh.e2)
    np.savetxt(path + 'e3.txt', model.mesh.e3)

    freestream = np.array([-10, 0, 0])
    sigma = -(model.mesh.e3[:, 0] * freestream[0] + model.mesh.e3[:, 1] * freestream[1] + model.mesh.e3[:, 2] * freestream[2])

    np.savetxt(path + 'freestream.txt', freestream)
    np.savetxt(path + 'sigma.txt', sigma)