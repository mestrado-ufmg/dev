import os
import vtk

from pybird.mesh.mesh import Mesh
from pybird.solver.solver import Solver

class View:

    def __init__(self, mesh: Mesh, solver: Solver, name: str) -> None:
        self._name = name
        self._mesh = mesh
        self._solver = solver
        return
    
    def paraview(self, path: str = '.', wake: bool = False) -> None:

        print('- Preparing view')
        
        pd = vtk.vtkPolyData()

        #----------------------------------#
        #               MESH               #
        #----------------------------------#
        points = vtk.vtkPoints()
        cells = vtk.vtkCellArray()

        # Surface
        for vertice in self._mesh.vertices:
            points.InsertNextPoint(vertice[0], vertice[1], vertice[2])
        
        for face in self._mesh.faces:
            cell = vtk.vtkTriangle()
            cell.GetPointIds().SetId(0, face[0])
            cell.GetPointIds().SetId(1, face[1])
            cell.GetPointIds().SetId(2, face[2])
            cells.InsertNextCell(cell)
        
        # Wake
        if wake:

            m = len(self._mesh.vertices[:, 0])
            lines = vtk.vtkCellArray()

            # Left wing
            nlw = len(self._mesh.wake.leftWing.vertices[:, 0])
            for i in range(nlw):
                points.InsertNextPoint(self._mesh.wake.leftWing.vertices[i, 0], self._mesh.wake.leftWing.vertices[i, 1], self._mesh.wake.leftWing.vertices[i, 2])

            nGridWake = len(self._mesh.wake.leftWing.grid[0, :])
            nGridSpan = len(self._mesh.wake.leftWing.grid[:, 0])
            for i in range(nGridSpan):
                polyLine = vtk.vtkPolyLine()
                polyLine.GetPointIds().SetNumberOfIds(nGridWake)
                for j in range(nGridWake):
                    polyLine.GetPointIds().SetId(j, m + self._mesh.wake.leftWing.grid[i, j])
                lines.InsertNextCell(polyLine)
            
            # Right wing
            nrw = len(self._mesh.wake.rightWing.vertices[:, 0])
            for i in range(nrw):
                points.InsertNextPoint(self._mesh.wake.rightWing.vertices[i, 0], self._mesh.wake.rightWing.vertices[i, 1], self._mesh.wake.rightWing.vertices[i, 2])

            nGridWake = len(self._mesh.wake.rightWing.grid[0, :])
            nGridSpan = len(self._mesh.wake.rightWing.grid[:, 0])
            for i in range(nGridSpan):
                polyLine = vtk.vtkPolyLine()
                polyLine.GetPointIds().SetNumberOfIds(nGridWake)
                for j in range(nGridWake):
                    polyLine.GetPointIds().SetId(j, m + nlw + self._mesh.wake.rightWing.grid[i, j])
                lines.InsertNextCell(polyLine)
            
            # Tail
            nt = len(self._mesh.wake.tail.vertices[:, 0])
            for i in range(nt):
                points.InsertNextPoint(self._mesh.wake.tail.vertices[i, 0], self._mesh.wake.tail.vertices[i, 1], self._mesh.wake.tail.vertices[i, 2])

            nGridWake = len(self._mesh.wake.tail.grid[0, :])
            nGridSpan = len(self._mesh.wake.tail.grid[:, 0])
            for i in range(nGridSpan):
                polyLine = vtk.vtkPolyLine()
                polyLine.GetPointIds().SetNumberOfIds(nGridWake)
                for j in range(nGridWake):
                    polyLine.GetPointIds().SetId(j, m + nlw + nrw + self._mesh.wake.tail.grid[i, j])
                lines.InsertNextCell(polyLine)

        
        pd.SetPoints(points)
        pd.SetPolys(cells)
        if wake: pd.SetLines(lines)

        #----------------------------------#
        #             VERTICES             #
        #----------------------------------#

        if self._solver.done:
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

            for i in range(self._solver.verticesParams.sigma.shape[0]):
                sigma.InsertNextTuple1(self._solver.verticesParams.sigma[i])
                doublet.InsertNextTuple1(self._solver.verticesParams.doublet[i])
                velNorm.InsertNextTuple1(self._solver.verticesParams.vel_norm[i])
                velField.InsertNextTuple3(self._solver.verticesParams.vel_field[i, 0], self._solver.verticesParams.vel_field[i, 1], self._solver.verticesParams.vel_field[i, 2])
                cp.InsertNextTuple1(self._solver.verticesParams.cp[i])
                transpiration.InsertNextTuple1(self._solver.verticesParams.transpiration[i])
                delta.InsertNextTuple1(self._solver.verticesParams.delta[i])
                A.InsertNextTuple1(self._solver.verticesParams.A[i])
                B.InsertNextTuple1(self._solver.verticesParams.B[i])
                Psi.InsertNextTuple1(self._solver.verticesParams.Psi[i])
                Ctau1.InsertNextTuple1(self._solver.verticesParams.Ctau1[i])
                Ctau2.InsertNextTuple1(self._solver.verticesParams.Ctau2[i])
                tauWall.InsertNextTuple3(self._solver.verticesParams.tau_wall[i, 0], self._solver.verticesParams.tau_wall[i, 1], self._solver.verticesParams.tau_wall[i, 2])
                amplification.InsertNextTuple1((self._solver.verticesParams.Ctau1[i] ** 2 + self._solver.verticesParams.Ctau2[i] ** 2) ** 0.5)

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
        writer.SetFileName(path + '/{}.vtp'.format(self._name))
        writer.SetInputData(pd)
        writer.Write()

        # Show using paraview
        os.system('paraview {}'.format(path + '/{}.vtp'.format(self._name)))

        return