import sys


sys.path.append('./')

import pybird

if __name__ == '__main__':
    
    # Create model
    model = pybird.model('test')

    # Load data
    model.load('./data/data.case')

    # Create geometry
    model.geo.build()

    # Create mesh
    model.mesh.build(size=0.04, alpha=20, accom_dist=10)

    # View
    model.view.mesh()