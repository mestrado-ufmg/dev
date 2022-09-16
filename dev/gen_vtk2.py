import numpy as np
import vtk

from lib_wrapper import VerticesParams

def genVTK(name: str,
           save_in: str,
           faces: np.ndarray,
           vertices: np.ndarray,
           params: VerticesParams):
    
    # Data structure
    pd = vtk.vtkPolyData()
    
    #----------------------------------#
    #               MESH               #
    #----------------------------------#

    # Parameters
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
    
    #----------------------------------#
    #             VERTICES             #
    #----------------------------------#

    # Parameters
    sigma = vtk.vtkFloatArray(); sigma.SetName('sigma')
    doublet = vtk.vtkFloatArray(); doublet.SetName('doublet')
    velNorm = vtk.vtkFloatArray(); velNorm.SetName('velNorm')
    velField = vtk.vtkFloatArray(); velField.SetNumberOfComponents(3); velField.SetName('velField')
    cp = vtk.vtkFloatArray(); cp.SetName('cp')
    transpiration = vtk.vtkFloatArray(); transpiration.SetName('transpiration')
    delta = vtk.vtkFloatArray(); delta.SetName('delta')
    A = vtk.vtkFloatArray(); A.SetName('A')
    B = vtk.vtkFloatArray(); B.SetName('B')
    Psi = vtk.vtkFloatArray(); Psi.SetName('Psi')
    Ctau1 = vtk.vtkFloatArray(); Ctau1.SetName('Ctau1')
    Ctau2 = vtk.vtkFloatArray(); Ctau2.SetName('Ctau2')
    tauWall = vtk.vtkFloatArray(); tauWall.SetNumberOfComponents(3); tauWall.SetName('tauWall')
    amplification = vtk.vtkFloatArray(); amplification.SetName('amplification')

    for i in range(params.sigma.shape[0]):
        sigma.InsertNextTuple1(params.sigma[i])
        doublet.InsertNextTuple1(params.doublet[i])
        velNorm.InsertNextTuple1(params.vel_norm[i])
        velField.InsertNextTuple3(params.vel_field[i, 0], params.vel_field[i, 1], params.vel_field[i, 2])
        cp.InsertNextTuple1(params.cp[i])
        transpiration.InsertNextTuple1(params.transpiration[i])
        delta.InsertNextTuple1(params.delta[i])
        A.InsertNextTuple1(params.A[i])
        B.InsertNextTuple1(params.B[i])
        Psi.InsertNextTuple1(params.Psi[i])
        Ctau1.InsertNextTuple1(params.Ctau1[i])
        Ctau2.InsertNextTuple1(params.Ctau2[i])
        tauWall.InsertNextTuple3(params.tau_wall[i, 0], params.tau_wall[i, 1], params.tau_wall[i, 2])
        amplification.InsertNextTuple1((params.Ctau1[i] ** 2 + params.Ctau2[i] ** 2) ** 0.5)

    pd.SetPoints(points)
    pd.SetPolys(cells)
    pd.GetPointData().AddArray(sigma)
    pd.GetPointData().AddArray(doublet)
    pd.GetPointData().AddArray(velNorm)
    pd.GetPointData().AddArray(velField)
    pd.GetPointData().AddArray(cp)
    pd.GetPointData().AddArray(transpiration)
    pd.GetPointData().AddArray(delta)
    pd.GetPointData().AddArray(A)
    pd.GetPointData().AddArray(B)
    pd.GetPointData().AddArray(Psi)
    pd.GetPointData().AddArray(Ctau1)
    pd.GetPointData().AddArray(Ctau2)
    pd.GetPointData().AddArray(tauWall)
    pd.GetPointData().AddArray(amplification)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(save_in + '{}.vtp'.format(name))
    writer.SetInputData(pd)
    writer.Write()
    
    return