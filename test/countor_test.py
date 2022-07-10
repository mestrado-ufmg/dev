import sys

import numpy as np

sys.path.append('./')

import pybird
from pybird.geo.utils import bezier

import matplotlib.pyplot as plt

if __name__ == '__main__':
    
    # Create model
    model = pybird.model('test')

    # Load data
    model.load('./data/data.case')

    # Create geometry and mesh
    model.geo.build()
    model.mesh.build()

    # Countor
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    # Surface points
    ax.scatter(
        [model.geo.p1e[0], model.geo.p2e[0], model.geo.p3e[0], model.geo.p4e[0], model.geo.p5e[0], model.geo.p6e[0], model.geo.p7e[0], model.geo.p8e[0], model.geo.p9e[0], model.geo.p10e[0], model.geo.p11e[0], model.geo.p12e[0], model.geo.p13e[0]],
        [model.geo.p1e[1], model.geo.p2e[1], model.geo.p3e[1], model.geo.p4e[1], model.geo.p5e[1], model.geo.p6e[1], model.geo.p7e[1], model.geo.p8e[1], model.geo.p9e[1], model.geo.p10e[1], model.geo.p11e[1], model.geo.p12e[1], model.geo.p13e[1]],
        [model.geo.p1e[2], model.geo.p2e[2], model.geo.p3e[2], model.geo.p4e[2], model.geo.p5e[2], model.geo.p6e[2], model.geo.p7e[2], model.geo.p8e[2], model.geo.p9e[2], model.geo.p10e[2], model.geo.p11e[2], model.geo.p12e[2], model.geo.p13e[2]],
        color='g',
    )

    ax.scatter(
        [model.geo.aux1e[0], model.geo.aux2e[0], model.geo.aux3e[0], model.geo.aux4e[0]],
        [model.geo.aux1e[1], model.geo.aux2e[1], model.geo.aux3e[1], model.geo.aux4e[1]],
        [model.geo.aux1e[2], model.geo.aux2e[2], model.geo.aux3e[2], model.geo.aux4e[2]],
        color='b',
    )

    ax.scatter(
        [model.geo.c1e[0], model.geo.c2e[0], model.geo.c3e[0], model.geo.c4e[0], model.geo.c5e[0], model.geo.c6e[0], model.geo.c7e[0], model.geo.c8e[0], model.geo.c9e[0], model.geo.c10e[0], model.geo.c11e[0], model.geo.c12e[0], model.geo.c13e[0], model.geo.c14e[0], model.geo.c15e[0], model.geo.c16e[0], model.geo.c17e[0], model.geo.c18e[0], model.geo.c19e[0], model.geo.c20e[0], model.geo.c21e[0]],
        [model.geo.c1e[1], model.geo.c2e[1], model.geo.c3e[1], model.geo.c4e[1], model.geo.c5e[1], model.geo.c6e[1], model.geo.c7e[1], model.geo.c8e[1], model.geo.c9e[1], model.geo.c10e[1], model.geo.c11e[1], model.geo.c12e[1], model.geo.c13e[1], model.geo.c14e[1], model.geo.c15e[1], model.geo.c16e[1], model.geo.c17e[1], model.geo.c18e[1], model.geo.c19e[1], model.geo.c20e[1], model.geo.c21e[1]],
        [model.geo.c1e[2], model.geo.c2e[2], model.geo.c3e[2], model.geo.c4e[2], model.geo.c5e[2], model.geo.c6e[2], model.geo.c7e[2], model.geo.c8e[2], model.geo.c9e[2], model.geo.c10e[2], model.geo.c11e[2], model.geo.c12e[2], model.geo.c13e[2], model.geo.c14e[2], model.geo.c15e[2], model.geo.c16e[2], model.geo.c17e[2], model.geo.c18e[2], model.geo.c19e[2], model.geo.c20e[2], model.geo.c21e[2]],
        color='r',
    )

    t = np.linspace(0, 1, 50)
    b1 = bezier.quadratic([model.geo.p1e, model.geo.c1e, model.geo.p2e], t)
    b2 = bezier.cubic([model.geo.p2e, model.geo.c2e, model.geo.c3e, model.geo.p3e], t)
    b3 = bezier.cubic([model.geo.p3e, model.geo.c4e, model.geo.c5e, model.geo.p4e], t)
    b4 = bezier.quadratic([model.geo.p4e, model.geo.p4e, model.geo.p5e], t)
    b5 = bezier.quadratic([model.geo.p5e, model.geo.aux2e, model.geo.p7e], t)
    b6 = bezier.quadratic([model.geo.p7e, model.geo.aux3e, model.geo.p10e], t)
    b7 = bezier.cubic([model.geo.p10e, model.geo.c16e, model.geo.c17e, model.geo.p11e], t)
    b8 = bezier.quadratic([model.geo.p11e, model.geo.aux4e, model.geo.p13e], t)
    b9 = bezier.cubic([model.geo.p5e, model.geo.c6e, model.geo.c7e, model.geo.p6e], t)
    b10 = bezier.cubic([model.geo.p6e, model.geo.c8e, model.geo.c9e, model.geo.p7e], t)
    b11 = bezier.cubic([model.geo.p7e, model.geo.c10e, model.geo.c11e, model.geo.p8e], t)
    b12 = bezier.cubic([model.geo.p8e, model.geo.c12e, model.geo.c13e, model.geo.p9e], t)
    b13 = bezier.cubic([model.geo.p9e, model.geo.c14e, model.geo.c15e, model.geo.p10e], t)
    b14 = bezier.cubic([model.geo.p11e, model.geo.c18e, model.geo.c19e, model.geo.p12e], t)
    b15 = bezier.cubic([model.geo.p12e, model.geo.c20e, model.geo.c21e, model.geo.p13e], t)

    ax.plot(b1[:, 0], b1[:, 1], b1[:, 2], 'k')
    ax.plot(b2[:, 0], b2[:, 1], b2[:, 2], 'k')
    ax.plot(b3[:, 0], b3[:, 1], b3[:, 2], 'k')
    ax.plot(b4[:, 0], b4[:, 1], b4[:, 2], 'k')
    ax.plot(b5[:, 0], b5[:, 1], b5[:, 2], '--k')
    ax.plot(b6[:, 0], b6[:, 1], b6[:, 2], '--k')
    ax.plot(b7[:, 0], b7[:, 1], b7[:, 2], 'k')
    ax.plot(b8[:, 0], b8[:, 1], b8[:, 2], '--k')
    ax.plot(b9[:, 0], b9[:, 1], b9[:, 2], 'y')
    ax.plot(b10[:, 0], b10[:, 1], b10[:, 2], 'y')
    ax.plot(b11[:, 0], b11[:, 1], b11[:, 2], 'y')
    ax.plot(b12[:, 0], b12[:, 1], b12[:, 2], 'y')
    ax.plot(b13[:, 0], b13[:, 1], b13[:, 2], 'y')
    ax.plot(b14[:, 0], b14[:, 1], b14[:, 2], 'y')
    ax.plot(b15[:, 0], b15[:, 1], b15[:, 2], 'y')

    ax.plot(model.geo.curve13e[:, 0], model.geo.curve13e[:, 1], model.geo.curve13e[:, 2], 'k')
    ax.plot(model.geo.curve14e[:, 0], model.geo.curve14e[:, 1], model.geo.curve14e[:, 2], 'k')
    ax.plot(model.geo.curve15e[:, 0], model.geo.curve15e[:, 1], model.geo.curve15e[:, 2], 'k')
    ax.plot(model.geo.curve16e[:, 0], model.geo.curve16e[:, 1], model.geo.curve16e[:, 2], 'k')
    ax.scatter([model.geo.p14e[0]], [model.geo.p14e[1]], [model.geo.p14e[2]], color='r')
    ax.scatter([model.geo.p15e[0]], [model.geo.p15e[1]], [model.geo.p15e[2]], color='r')

    ax.plot(model.geo.curve25e[:, 0], model.geo.curve25e[:, 1], model.geo.curve25e[:, 2], 'k')
    ax.plot(model.geo.curve26e[:, 0], model.geo.curve26e[:, 1], model.geo.curve26e[:, 2], 'k')
    ax.plot(model.geo.curve27e[:, 0], model.geo.curve27e[:, 1], model.geo.curve27e[:, 2], 'k')
    ax.plot(model.geo.curve28e[:, 0], model.geo.curve28e[:, 1], model.geo.curve28e[:, 2], 'k')
    ax.scatter([model.geo.p18e[0]], [model.geo.p18e[1]], [model.geo.p18e[2]], color='r')
    ax.scatter([model.geo.p19e[0]], [model.geo.p19e[1]], [model.geo.p19e[2]], color='r')

    ax.plot(model.geo.curve31e[:, 0], model.geo.curve31e[:, 1], model.geo.curve31e[:, 2], 'k')
    ax.plot(model.geo.curve32e[:, 0], model.geo.curve32e[:, 1], model.geo.curve32e[:, 2], 'k')
    ax.plot(model.geo.curve33e[:, 0], model.geo.curve33e[:, 1], model.geo.curve33e[:, 2], 'k')
    ax.plot(model.geo.curve34e[:, 0], model.geo.curve34e[:, 1], model.geo.curve34e[:, 2], 'k')
    ax.scatter([model.geo.p20e[0]], [model.geo.p20e[1]], [model.geo.p20e[2]], color='r')
    ax.scatter([model.geo.p21e[0]], [model.geo.p21e[1]], [model.geo.p21e[2]], color='r')

    ax.plot(model.geo.curve43e[:, 0], model.geo.curve43e[:, 1], model.geo.curve43e[:, 2], 'k')
    ax.plot(model.geo.curve44e[:, 0], model.geo.curve44e[:, 1], model.geo.curve44e[:, 2], 'k')
    ax.plot(model.geo.curve45e[:, 0], model.geo.curve45e[:, 1], model.geo.curve45e[:, 2], 'k')
    ax.plot(model.geo.curve46e[:, 0], model.geo.curve46e[:, 1], model.geo.curve46e[:, 2], 'k')
    ax.scatter([model.geo.p24e[0]], [model.geo.p24e[1]], [model.geo.p24e[2]], color='r')
    ax.scatter([model.geo.p25e[0]], [model.geo.p25e[1]], [model.geo.p25e[2]], color='r')

    ax.plot(model.geo.curve19e[:, 0], model.geo.curve19e[:, 1], model.geo.curve19e[:, 2], 'k')
    ax.plot(model.geo.curve20e[:, 0], model.geo.curve20e[:, 1], model.geo.curve20e[:, 2], 'k')
    ax.plot(model.geo.curve21e[:, 0], model.geo.curve21e[:, 1], model.geo.curve21e[:, 2], 'k')
    ax.plot(model.geo.curve22e[:, 0], model.geo.curve22e[:, 1], model.geo.curve22e[:, 2], 'k')
    ax.scatter([model.geo.p16e[0]], [model.geo.p16e[1]], [model.geo.p16e[2]], color='r')
    ax.scatter([model.geo.p17e[0]], [model.geo.p17e[1]], [model.geo.p17e[2]], color='r')

    ax.plot(model.geo.curve37e[:, 0], model.geo.curve37e[:, 1], model.geo.curve37e[:, 2], 'k')
    ax.plot(model.geo.curve38e[:, 0], model.geo.curve38e[:, 1], model.geo.curve38e[:, 2], 'k')
    ax.plot(model.geo.curve39e[:, 0], model.geo.curve39e[:, 1], model.geo.curve39e[:, 2], 'k')
    ax.plot(model.geo.curve40e[:, 0], model.geo.curve40e[:, 1], model.geo.curve40e[:, 2], 'k')
    ax.scatter([model.geo.p22e[0]], [model.geo.p22e[1]], [model.geo.p22e[2]], color='r')
    ax.scatter([model.geo.p23e[0]], [model.geo.p23e[1]], [model.geo.p23e[2]], color='r')

    ax.plot(model.geo.curve17e[:, 0], model.geo.curve17e[:, 1], model.geo.curve17e[:, 2], 'k')
    ax.plot(model.geo.curve18e[:, 0], model.geo.curve18e[:, 1], model.geo.curve18e[:, 2], 'k')
    ax.plot(model.geo.curve23e[:, 0], model.geo.curve23e[:, 1], model.geo.curve23e[:, 2], 'k')
    ax.plot(model.geo.curve24e[:, 0], model.geo.curve24e[:, 1], model.geo.curve24e[:, 2], 'k')
    ax.plot(model.geo.curve29e[:, 0], model.geo.curve29e[:, 1], model.geo.curve29e[:, 2], 'k')
    ax.plot(model.geo.curve30e[:, 0], model.geo.curve30e[:, 1], model.geo.curve30e[:, 2], 'k')
    ax.plot(model.geo.curve35e[:, 0], model.geo.curve35e[:, 1], model.geo.curve35e[:, 2], 'k')
    ax.plot(model.geo.curve36e[:, 0], model.geo.curve36e[:, 1], model.geo.curve36e[:, 2], 'k')
    ax.plot(model.geo.curve41e[:, 0], model.geo.curve41e[:, 1], model.geo.curve41e[:, 2], 'k')
    ax.plot(model.geo.curve42e[:, 0], model.geo.curve42e[:, 1], model.geo.curve42e[:, 2], 'k')
    ax.plot(model.geo.curve47e[:, 0], model.geo.curve47e[:, 1], model.geo.curve47e[:, 2], 'k')
    ax.plot(model.geo.curve48e[:, 0], model.geo.curve48e[:, 1], model.geo.curve48e[:, 2], 'k')
    
    ax.set_xlim(-0.300, 0.250)
    ax.set_ylim( 0.050, 0.500)
    ax.set_zlim(-0.275, 0.275)
    
    plt.show()