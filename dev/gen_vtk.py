import numpy as np
import vtk

def genVTK(name: str,
           path: str,
           save_in: str,
           sigma: np.ndarray,
           doublet: np.ndarray,
           vel: np.ndarray,
           velx: np.ndarray,
           vely: np.ndarray,
           velz: np.ndarray,
           cp: np.ndarray,
           transpiration: np.ndarray):

    vertices = np.loadtxt(path + 'vertices.txt')
    faces = np.loadtxt(path + 'faces.txt').astype(np.int64)

    pd = vtk.vtkPolyData()

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()

    # Surface
    for vertice in vertices:
        points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
        
    for face in faces:
        cell = vtk.vtkTriangle()
        cell.GetPointIds().SetId(0, face[0])
        cell.GetPointIds().SetId(1, face[1])
        cell.GetPointIds().SetId(2, face[2])
        cells.InsertNextCell(cell)
    
    # Parameters
    sigmaVTK = vtk.vtkFloatArray()
    sigmaVTK.SetName('sigma')

    doubletVTK = vtk.vtkFloatArray()
    doubletVTK.SetName('doublet')

    velNormVTK = vtk.vtkFloatArray()
    velNormVTK.SetName('velNorm')

    velFieldVTK = vtk.vtkFloatArray()
    velFieldVTK.SetNumberOfComponents(3)
    velFieldVTK.SetName('velField')

    cpVTK = vtk.vtkFloatArray()
    cpVTK.SetName('cp')

    transpirationVTK = vtk.vtkFloatArray()
    transpirationVTK.SetName('transpiration')
    
    n = len(vertices[:, 0])

    sigmaVerticesValues = np.zeros(n)
    doubletVerticesValues = np.zeros(n)
    velNormVerticesValues = np.zeros(n)
    velFieldVerticesValues = np.zeros((n, 3))
    cpVerticesValues = np.zeros(n)
    transpirationVerticesValues = np.zeros(n)
    nVerticesValues = np.zeros(n)

    for i in range(len(faces[:, 0])):
        sigmaVerticesValues[faces[i, 0]] = sigmaVerticesValues[faces[i, 0]] + sigma[i]
        sigmaVerticesValues[faces[i, 1]] = sigmaVerticesValues[faces[i, 1]] + sigma[i]
        sigmaVerticesValues[faces[i, 2]] = sigmaVerticesValues[faces[i, 2]] + sigma[i]
                
        doubletVerticesValues[faces[i, 0]] = doubletVerticesValues[faces[i, 0]] + doublet[i]
        doubletVerticesValues[faces[i, 1]] = doubletVerticesValues[faces[i, 1]] + doublet[i]
        doubletVerticesValues[faces[i, 2]] = doubletVerticesValues[faces[i, 2]] + doublet[i]
                
        velNormVerticesValues[faces[i, 0]] = velNormVerticesValues[faces[i, 0]] + vel[i]
        velNormVerticesValues[faces[i, 1]] = velNormVerticesValues[faces[i, 1]] + vel[i]
        velNormVerticesValues[faces[i, 2]] = velNormVerticesValues[faces[i, 2]] + vel[i]

        velFieldVerticesValues[faces[i, 0], 0] = velFieldVerticesValues[faces[i, 0], 0] + velx[i]
        velFieldVerticesValues[faces[i, 1], 0] = velFieldVerticesValues[faces[i, 1], 0] + velx[i]
        velFieldVerticesValues[faces[i, 2], 0] = velFieldVerticesValues[faces[i, 2], 0] + velx[i]
        velFieldVerticesValues[faces[i, 0], 1] = velFieldVerticesValues[faces[i, 0], 1] + vely[i]
        velFieldVerticesValues[faces[i, 1], 1] = velFieldVerticesValues[faces[i, 1], 1] + vely[i]
        velFieldVerticesValues[faces[i, 2], 1] = velFieldVerticesValues[faces[i, 2], 1] + vely[i]
        velFieldVerticesValues[faces[i, 0], 2] = velFieldVerticesValues[faces[i, 0], 2] + velz[i]
        velFieldVerticesValues[faces[i, 1], 2] = velFieldVerticesValues[faces[i, 1], 2] + velz[i]
        velFieldVerticesValues[faces[i, 2], 2] = velFieldVerticesValues[faces[i, 2], 2] + velz[i]

        cpVerticesValues[faces[i, 0]] = cpVerticesValues[faces[i, 0]] + cp[i]
        cpVerticesValues[faces[i, 1]] = cpVerticesValues[faces[i, 1]] + cp[i]
        cpVerticesValues[faces[i, 2]] = cpVerticesValues[faces[i, 2]] + cp[i]

        transpirationVerticesValues[faces[i, 0]] = transpirationVerticesValues[faces[i, 0]] + transpiration[i]
        transpirationVerticesValues[faces[i, 1]] = transpirationVerticesValues[faces[i, 1]] + transpiration[i]
        transpirationVerticesValues[faces[i, 2]] = transpirationVerticesValues[faces[i, 2]] + transpiration[i]

        nVerticesValues[faces[i, 0]] = nVerticesValues[faces[i, 0]] + 1
        nVerticesValues[faces[i, 1]] = nVerticesValues[faces[i, 1]] + 1
        nVerticesValues[faces[i, 2]] = nVerticesValues[faces[i, 2]] + 1
            
    sigmaValues = sigmaVerticesValues / (nVerticesValues + 1e-8)
    doubletValues = doubletVerticesValues / (nVerticesValues + 1e-8)
    velNormValues = velNormVerticesValues / (nVerticesValues + 1e-8)
    velxFieldValues = velFieldVerticesValues[:, 0] / (nVerticesValues + 1e-8)
    velyFieldValues = velFieldVerticesValues[:, 1] / (nVerticesValues + 1e-8)
    velzFieldValues = velFieldVerticesValues[:, 2] / (nVerticesValues + 1e-8)
    cpValues = cpVerticesValues / (nVerticesValues + 1e-8)
    transpirationValues = transpirationVerticesValues / (nVerticesValues + 1e-8)

    for i in range(n):
        sigmaVTK.InsertNextTuple1(sigmaValues[i])
        doubletVTK.InsertNextTuple1(doubletValues[i])
        velNormVTK.InsertNextTuple1(velNormValues[i])
        velFieldVTK.InsertNextTuple3(velxFieldValues[i], velyFieldValues[i], velzFieldValues[i])
        cpVTK.InsertNextTuple1(cpValues[i])
        transpirationVTK.InsertNextTuple1(transpirationValues[i])

    pd.SetPoints(points)
    pd.SetPolys(cells)
    pd.GetPointData().AddArray(sigmaVTK)
    pd.GetPointData().AddArray(doubletVTK)
    pd.GetPointData().AddArray(velNormVTK)
    pd.GetPointData().AddArray(velFieldVTK)
    pd.GetPointData().AddArray(cpVTK)
    pd.GetPointData().AddArray(transpirationVTK)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(save_in + '{}.vtp'.format(name))
    writer.SetInputData(pd)
    writer.Write()
    
    return