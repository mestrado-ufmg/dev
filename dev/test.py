import numpy as np

from solver_wrapper import solve
from gen_vtk import gen

if __name__ == '__main__':

    path = './data/mesh/'

    # Load
    facesAreas = np.loadtxt(path + 'facesAreas.txt')
    facesMaxDistance = np.loadtxt(path + 'facesMaxDistance.txt')
    facesCenter = np.loadtxt(path + 'facesCenter.txt')
    controlPoints = np.loadtxt(path + 'controlPoints.txt')
    p1 = np.loadtxt(path + 'p1.txt')
    p2 = np.loadtxt(path + 'p2.txt')
    p3 = np.loadtxt(path + 'p3.txt')
    e1 = np.loadtxt(path + 'e1.txt')
    e2 = np.loadtxt(path + 'e2.txt')
    e3 = np.loadtxt(path + 'e3.txt')
    freestream = np.loadtxt(path + 'freestream.txt')
    sigma = np.loadtxt(path + 'sigma.txt')

    # Call function
    out = solve(facesAreas,
                facesMaxDistance,
                facesCenter,
                controlPoints,
                p1,
                p2,
                p3,
                e1,
                e2,
                e3,
                freestream,
                sigma)

    n = facesAreas.shape[0]

    transpiration = out.velx * e3[:, 0] + out.vely * e3[:, 1] + out.velz * e3[:, 2]

    gen('first',
        './data/mesh/',
        './data/vtk/',
        sigma,
        out.doublet,
        out.velNorm,
        out.velx,
        out.vely,
        out.velz,
        out.cp,
        transpiration)