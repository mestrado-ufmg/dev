import numpy as np

if __name__ == '__main__':

    area = 10.0
    force = np.array([0.0, 0.0, 0.0])
    x = np.array([1.0, 0.0, 0.0])
    y = np.array([0.0, 1.0, 0.0])
    z = np.array([0.0, 0.0, 1.0])

    facesAreas = np.loadtxt('./data/mesh/NACA0012-AoA-0/facesAreas.txt', dtype=np.double)
    e3 = np.loadtxt('./data/mesh/NACA0012-AoA-0/e3.txt', dtype=np.double)
    cp = np.loadtxt('./data/solution/NACA0012-AoA-0/cp_f.txt', dtype=np.double)

    for i in range(cp.size):
        force = force - cp[i] * e3[i, :] * facesAreas[i]
    
    print(force / area)