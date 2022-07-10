import sys

sys.path.append('./')

import pybird

import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    # Create model
    model = pybird.model('test')

    # Load data
    model.load('./data/data.case')

    # Clean shoulder angles
    model.data.geo.wing.theta0 = 0
    model.data.geo.wing.theta1_e = 0
    model.data.geo.wing.theta2_e = 0
    model.data.geo.wing.theta3_e = 0

    # Build geometry
    model.geo.build()

    # Countor
    fig = plt.figure()

    plt.scatter(
        [model.geo.p1e[0], model.geo.p13e[0], model.geo.p14e[0], model.geo.p15e[0]],
        [model.geo.p1e[2], model.geo.p13e[2], model.geo.p14e[2], model.geo.p15e[2]],
        color='k'
    )

    plt.scatter(model.geo.curve13e[:, 0], model.geo.curve13e[:, 2], color='r')
    plt.scatter(model.geo.curve15e[:, 0], model.geo.curve15e[:, 2], color='b')
    plt.scatter(model.geo.curve14e[:, 0], model.geo.curve14e[:, 2], color='g')
    plt.scatter(model.geo.curve16e[:, 0], model.geo.curve16e[:, 2], color='y')
    
    plt.axis('equal')
    plt.grid()

    plt.show()