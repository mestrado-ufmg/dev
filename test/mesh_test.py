import sys
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
        size=0.06,
        accom_dist=10,
        n_head=3,
        n_wing_le=3,
        n_wing_te=3,
        n_tail_le=3,
        n_tail_te=3,
    )

    # View
    model.view.mesh()