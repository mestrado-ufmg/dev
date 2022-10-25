import numpy as np

from utils.mesh import gen_mesh
from utils.vtk import gen_vtk_file



if __name__ == '__main__':

    # Parameters
    span = 3.0
    chord = 1.0
    foil = np.loadtxt('./foils/NACA0012.dat')
    cell_size = chord / 10
    le_ratio = 0.1
    te_ratio = 0.2
    alpha = 0.0
    wake_lenght = 20
    wake_accom_dist = 3

    # Mesh
    mesh = gen_mesh(foil, span, chord, cell_size, le_ratio, te_ratio, alpha, wake_lenght, wake_accom_dist)

    path = './data/mesh/NACA0012-AoA-0/'

    # Save
    np.savetxt(path + 'controlPoints.txt', mesh.controlPoints)
    np.savetxt(path + 'e1.txt', mesh.e1)
    np.savetxt(path + 'e2.txt', mesh.e2)
    np.savetxt(path + 'e3.txt', mesh.e3)
    np.savetxt(path + 'edges.txt', mesh.edges)
    np.savetxt(path + 'faces.txt', mesh.faces)
    np.savetxt(path + 'facesAreas.txt', mesh.facesAreas)
    np.savetxt(path + 'facesCenter.txt', mesh.facesCenter)
    np.savetxt(path + 'facesMaxDistance.txt', mesh.facesMaxDistance)
    np.savetxt(path + 'p1Local.txt', mesh.p1Local)
    np.savetxt(path + 'p2Local.txt', mesh.p2Local)
    np.savetxt(path + 'p3Local.txt', mesh.p3Local)
    np.savetxt(path + 'vertices.txt', mesh.vertices)
    np.savetxt(path + 'wake_faces.txt', mesh.wake_faces)
    np.savetxt(path + 'wake_grid.txt', mesh.wake_grid)
    np.savetxt(path + 'wake_vertices.txt', mesh.wake_vertices)

    # Generate vtk file
    gen_vtk_file('./data/vtk/geometry-NACA0012-AoA-0', mesh, wake=False)