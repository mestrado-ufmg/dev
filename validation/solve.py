import numpy as np

from utils.mesh import MeshData
from utils.solver import solve
from utils.vtk import gen_vtk_file

#---------------------------------------------#
#                   POSPROC                   #
#---------------------------------------------#
def getForces():

    from scipy.integrate import dblquad

    def dJ( u, v, p1, p2, p3 ):
        x1, y1 = p1
        x2, y2 = p2
        x3, y3 = p3
        dxdu = ( (1-v)*x2 + v*x3 - x1 )
        dxdv = ( u*x3 - u*x2 )
        dydu = ( (1-v)*y2 + v*y3 - y1 )
        dydv = ( u*y3 - u*y2 )
        return np.abs( dxdu*dydv - dxdv*dydu )

    def tridblquad( integrand, p1, p2, p3 ):
        '''
        Perform double quadtrature integration on triangular domain.
        Input: function to integrate, points of triangle as tuples.
        Output: integral and estimated absolute error as a tuple.
        '''
        x1, y1 = p1 ; x2, y2 = p2 ; x3, y3 = p3

        # transformation to the unit square
        g = lambda u, v, c1, c2, c3: (1-u)*c1 + u*( (1-v)*c2 + v*c3 )

        # transformation for the integrand, 
        # including the Jacobian scaling factor
        def h( u, v ):
            x = g( u, v, x1, x2, x3 )
            y = g( u, v, y1, y2, y3 )
            I = integrand( x, y )
            I *= dJ( u, v, p1, p2, p3 )
            return I
        
        # perfrom the double integration using quadrature in the transformed space
        integral, error = dblquad( h, 0, 1, lambda x: 0, lambda x: 1, epsrel=1e-6, epsabs=0 )
        return integral, error

    v1 = 0.0
    v2 = 0.1
    v3 = 0.15

    p1 = (0., 0.)
    p2 = (1., 0.)
    p3 = (0., 1.)

    area, error = tridblquad(lambda x, y: 2., p1, p2, p3)
    print('Area: {},   Error: {:.3e}'.format( area, error ))


    return

if __name__ == '__main__':

    # getForces()

    alpha = 0.0
    freestream = 10.0
    density = 1.225
    viscosity = 1e-5
    soundSpeed = 340

    path = './data/mesh/NACA0012-AoA-0/'

    mesh = MeshData(
        vertices=np.loadtxt(path + 'vertices.txt', dtype=np.double),
        edges=np.loadtxt(path + 'edges.txt', dtype=np.int32),
        faces=np.loadtxt(path + 'faces.txt', dtype=np.int32),
        facesCenter=np.loadtxt(path + 'facesCenter.txt', dtype=np.double),
        e1=np.loadtxt(path + 'e1.txt', dtype=np.double),
        e2=np.loadtxt(path + 'e2.txt', dtype=np.double),
        e3=np.loadtxt(path + 'e3.txt', dtype=np.double),
        controlPoints=np.loadtxt(path + 'controlPoints.txt', dtype=np.double),
        facesAreas=np.loadtxt(path + 'facesAreas.txt', dtype=np.double),
        facesMaxDistance=np.loadtxt(path + 'facesMaxDistance.txt', dtype=np.double),
        p1Local=np.loadtxt(path + 'p1Local.txt', dtype=np.double),
        p2Local=np.loadtxt(path + 'p2Local.txt', dtype=np.double),
        p3Local=np.loadtxt(path + 'p3Local.txt', dtype=np.double),
        wake_faces=np.loadtxt(path + 'wake_faces.txt', dtype=np.int32),
        wake_grid=np.loadtxt(path + 'wake_grid.txt', dtype=np.int32),
        wake_vertices=np.loadtxt(path + 'wake_vertices.txt', dtype=np.double),
    )

    sol = solve(mesh, freestream, alpha, density, viscosity, soundSpeed)

    # vtk
    gen_vtk_file('./data/vtk/solution-NACA0012-AoA-0', mesh, sol)

    # Save
    path = './data/solution/NACA0012-AoA-0/'
    
    np.savetxt(path + 'cp_v.txt', sol.cp_v)
    np.savetxt(path + 'velx_v.txt', sol.velx_v)
    np.savetxt(path + 'vely_v.txt', sol.vely_v)
    np.savetxt(path + 'velz_v.txt', sol.velz_v)
    np.savetxt(path + 'transpiration_v.txt', sol.transpiration_v)
    np.savetxt(path + 'sigma_v.txt', sol.sigma_v)
    np.savetxt(path + 'doublet_v.txt', sol.doublet_v)