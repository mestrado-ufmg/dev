import os
from numpy import array, zeros
import vtk

from pybird.mesh.mesh import Mesh
from pybird.solver.solver import Solver

class View:

    def __init__(self, mesh: Mesh, solver: Solver, name: str) -> None:
        self._name = name
        self._mesh = mesh
        self._solver = solver
        return
    
    def paraview(self, path: str = '.', showWake: bool = False) -> None:

        print('- Preparing view')
        
        pd = vtk.vtkPolyData()

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
        if showWake:

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
        if showWake: pd.SetLines(lines)

        # Solver
        if self._solver.sigma is not None:
            sigma = vtk.vtkFloatArray()
            sigma.SetName('sigma')

            doublet = vtk.vtkFloatArray()
            doublet.SetName('doublet')

            velNorm = vtk.vtkFloatArray()
            velNorm.SetName('velNorm')

            velField = vtk.vtkFloatArray()
            velField.SetNumberOfComponents(3)
            velField.SetName('velField')

            cp = vtk.vtkFloatArray()
            cp.SetName('cp')

            transpiration = vtk.vtkFloatArray()
            transpiration.SetName('transpiration')

            n = len(self._mesh.vertices[:, 0])
            sigmaVerticesValues = zeros(n)
            doubletVerticesValues = zeros(n)
            velNormVerticesValues = zeros(n)
            velFieldVerticesValues = zeros((n, 3))
            cpVerticesValues = zeros(n)
            transpirationVerticesValues = zeros(n)
            nVerticesValues = zeros(n)
            

            for i in range(len(self._mesh.faces[:, 0])):
                sigmaVerticesValues[self._mesh.faces[i, 0]] = sigmaVerticesValues[self._mesh.faces[i, 0]] + self._solver.sigma[i]
                sigmaVerticesValues[self._mesh.faces[i, 1]] = sigmaVerticesValues[self._mesh.faces[i, 1]] + self._solver.sigma[i]
                sigmaVerticesValues[self._mesh.faces[i, 2]] = sigmaVerticesValues[self._mesh.faces[i, 2]] + self._solver.sigma[i]
                
                doubletVerticesValues[self._mesh.faces[i, 0]] = doubletVerticesValues[self._mesh.faces[i, 0]] + self._solver.doublet[i]
                doubletVerticesValues[self._mesh.faces[i, 1]] = doubletVerticesValues[self._mesh.faces[i, 1]] + self._solver.doublet[i]
                doubletVerticesValues[self._mesh.faces[i, 2]] = doubletVerticesValues[self._mesh.faces[i, 2]] + self._solver.doublet[i]
                
                velNormVerticesValues[self._mesh.faces[i, 0]] = velNormVerticesValues[self._mesh.faces[i, 0]] + self._solver.velNorm[i]
                velNormVerticesValues[self._mesh.faces[i, 1]] = velNormVerticesValues[self._mesh.faces[i, 1]] + self._solver.velNorm[i]
                velNormVerticesValues[self._mesh.faces[i, 2]] = velNormVerticesValues[self._mesh.faces[i, 2]] + self._solver.velNorm[i]

                velFieldVerticesValues[self._mesh.faces[i, 0], 0] = velFieldVerticesValues[self._mesh.faces[i, 0], 0] + self._solver.velField[i, 0]
                velFieldVerticesValues[self._mesh.faces[i, 1], 0] = velFieldVerticesValues[self._mesh.faces[i, 1], 0] + self._solver.velField[i, 0]
                velFieldVerticesValues[self._mesh.faces[i, 2], 0] = velFieldVerticesValues[self._mesh.faces[i, 2], 0] + self._solver.velField[i, 0]
                velFieldVerticesValues[self._mesh.faces[i, 0], 1] = velFieldVerticesValues[self._mesh.faces[i, 0], 1] + self._solver.velField[i, 1]
                velFieldVerticesValues[self._mesh.faces[i, 1], 1] = velFieldVerticesValues[self._mesh.faces[i, 1], 1] + self._solver.velField[i, 1]
                velFieldVerticesValues[self._mesh.faces[i, 2], 1] = velFieldVerticesValues[self._mesh.faces[i, 2], 1] + self._solver.velField[i, 1]
                velFieldVerticesValues[self._mesh.faces[i, 0], 2] = velFieldVerticesValues[self._mesh.faces[i, 0], 2] + self._solver.velField[i, 2]
                velFieldVerticesValues[self._mesh.faces[i, 1], 2] = velFieldVerticesValues[self._mesh.faces[i, 1], 2] + self._solver.velField[i, 2]
                velFieldVerticesValues[self._mesh.faces[i, 2], 2] = velFieldVerticesValues[self._mesh.faces[i, 2], 2] + self._solver.velField[i, 2]

                cpVerticesValues[self._mesh.faces[i, 0]] = cpVerticesValues[self._mesh.faces[i, 0]] + self._solver.cp[i]
                cpVerticesValues[self._mesh.faces[i, 1]] = cpVerticesValues[self._mesh.faces[i, 1]] + self._solver.cp[i]
                cpVerticesValues[self._mesh.faces[i, 2]] = cpVerticesValues[self._mesh.faces[i, 2]] + self._solver.cp[i]

                nVel = self._solver.velField[i, 0] * self._mesh.e3[i, 0] + self._solver.velField[i, 1] * self._mesh.e3[i, 1] + self._solver.velField[i, 2] * self._mesh.e3[i, 2]
                transpirationVerticesValues[self._mesh.faces[i, 0]] = transpirationVerticesValues[self._mesh.faces[i, 0]] + nVel
                transpirationVerticesValues[self._mesh.faces[i, 1]] = transpirationVerticesValues[self._mesh.faces[i, 1]] + nVel
                transpirationVerticesValues[self._mesh.faces[i, 2]] = transpirationVerticesValues[self._mesh.faces[i, 2]] + nVel

                nVerticesValues[self._mesh.faces[i, 0]] = nVerticesValues[self._mesh.faces[i, 0]] + 1
                nVerticesValues[self._mesh.faces[i, 1]] = nVerticesValues[self._mesh.faces[i, 1]] + 1
                nVerticesValues[self._mesh.faces[i, 2]] = nVerticesValues[self._mesh.faces[i, 2]] + 1
            
            sigmaValues = sigmaVerticesValues / (nVerticesValues + 1e-8)
            doubletValues = doubletVerticesValues / (nVerticesValues + 1e-8)
            velNormValues = velNormVerticesValues / (nVerticesValues + 1e-8)
            velxFieldValues = velFieldVerticesValues[:, 0] / (nVerticesValues + 1e-8)
            velyFieldValues = velFieldVerticesValues[:, 1] / (nVerticesValues + 1e-8)
            velzFieldValues = velFieldVerticesValues[:, 2] / (nVerticesValues + 1e-8)
            cpValues = cpVerticesValues / (nVerticesValues + 1e-8)
            transpirationValues = transpirationVerticesValues / (nVerticesValues + 1e-8)

            for i in range(n):
                sigma.InsertNextTuple1(sigmaValues[i])
                doublet.InsertNextTuple1(doubletValues[i])
                velNorm.InsertNextTuple1(velNormValues[i])
                velField.InsertNextTuple3(velxFieldValues[i], velyFieldValues[i], velzFieldValues[i])
                cp.InsertNextTuple1(cpValues[i])
                transpiration.InsertNextTuple1(transpirationValues[i])

            pd.GetPointData().AddArray(sigma)
            pd.GetPointData().AddArray(doublet)
            pd.GetPointData().AddArray(velNorm)
            pd.GetPointData().AddArray(velField)
            pd.GetPointData().AddArray(cp)
            pd.GetPointData().AddArray(transpiration)

        writer = vtk.vtkXMLPolyDataWriter()
        writer.SetFileName(path + '/{}.vtp'.format(self._name))
        writer.SetInputData(pd)
        writer.Write()

        # Show using paraview
        os.system('paraview {}'.format(path + '/{}.vtp'.format(self._name)))

        return