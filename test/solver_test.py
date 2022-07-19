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
        size=0.16,
        accom_dist=10,
        n_head=3,
        n_wing_le=3,
        n_wing_te=3,
        n_tail_le=3,
        n_tail_te=3,
    )

    # Solve
    model.solver.solve(1)

    # plt.figure()
    # plt.plot(model.solver.sigma)
    # plt.grid()
    # plt.show()
    
    # View
    # model.view.surfaceParameter(model.solver.velNorm)
    model.view.faceVectors([model.solver.velField])