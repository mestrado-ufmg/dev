import sys

import numpy as np
sys.path.append('./')

import pybird

import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    # Create model
    model = pybird.model('test')

    # Load data
    model.load('./data/data.case')

    # Create geometry
    model.geo.build()

    # Create mesh
    model.mesh.build(
        size=0.15,
        accom_dist=10,
        n_head=3,
        n_wing_le=3,
        n_wing_te=3,
        n_tail_le=3,
        n_tail_te=3,
    )

    d = np.zeros_like(model.mesh.e1)
    d[:, 0] = model.mesh.e1[:, 0] * model.mesh.e2[:, 0] + model.mesh.e1[:, 1] * model.mesh.e2[:, 1] + model.mesh.e1[:, 2] * model.mesh.e2[:, 2]
    d[:, 1] = model.mesh.e3[:, 0] * model.mesh.e2[:, 0] + model.mesh.e3[:, 1] * model.mesh.e2[:, 1] + model.mesh.e3[:, 2] * model.mesh.e2[:, 2]
    d[:, 1] = model.mesh.e3[:, 0] * model.mesh.e1[:, 0] + model.mesh.e3[:, 1] * model.mesh.e1[:, 1] + model.mesh.e3[:, 2] * model.mesh.e1[:, 2]

    print(d)

    # View
    model.view.faceVectors()