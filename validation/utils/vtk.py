import numpy as np
import vtk
from utils.mesh import MeshData

from utils.solver import SolverData

def gen_vtk_file(file: str, mesh: MeshData, sol: SolverData = None, wake: bool = False) -> None:

    pd = vtk.vtkPolyData()

    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()

    for vertice in mesh.vertices:
        points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
        
    for face in mesh.faces:
        cell = vtk.vtkTriangle()
        cell.GetPointIds().SetId(0, face[0])
        cell.GetPointIds().SetId(1, face[1])
        cell.GetPointIds().SetId(2, face[2])
        cells.InsertNextCell(cell)
    
    pd.SetPoints(points)
    pd.SetPolys(cells)

    if wake:

        m = len(mesh.vertices[:, 0])
        lines = vtk.vtkCellArray()

        nlw = len(mesh.wake_vertices[:, 0])
        for i in range(nlw):
            points.InsertNextPoint(mesh.wake_vertices[i, 0], mesh.wake_vertices[i, 1], mesh.wake_vertices[i, 2])

        nGridWake = len(mesh.wake_grid[0, :])
        nGridSpan = len(mesh.wake_grid[:, 0])
        for i in range(nGridSpan):
            polyLine = vtk.vtkPolyLine()
            polyLine.GetPointIds().SetNumberOfIds(nGridWake)
            for j in range(nGridWake):
                polyLine.GetPointIds().SetId(j, m + mesh.wake_grid[i, j])
            lines.InsertNextCell(polyLine)
        
        pd.SetLines(lines)

    if sol is not None:
        
        cp = vtk.vtkFloatArray()
        cp.SetName('Cp')

        transpiration = vtk.vtkFloatArray()
        transpiration.SetName('Transpiration')

        vel = vtk.vtkFloatArray()
        vel.SetNumberOfComponents(3)
        vel.SetName('Velocity')

        sigma = vtk.vtkFloatArray()
        sigma.SetName('Sigma')

        doublet = vtk.vtkFloatArray()
        doublet.SetName('Doublet')

        for i in range(sol.cp_v.shape[0]):
            cp.InsertNextTuple1(sol.cp_v[i])
            vel.InsertNextTuple3(sol.velx_v[i], sol.vely_v[i], sol.velz_v[i])
            transpiration.InsertNextTuple1(sol.transpiration_v[i])
            sigma.InsertNextTuple1(sol.sigma_v[i])
            doublet.InsertNextTuple1(sol.doublet_v[i])
        
        pd.GetPointData().AddArray(cp)
        pd.GetPointData().AddArray(vel)
        pd.GetPointData().AddArray(transpiration)
        pd.GetPointData().AddArray(sigma)
        pd.GetPointData().AddArray(doublet)
    
    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(file + '.vtp')
    writer.SetInputData(pd)
    writer.Write()

    return