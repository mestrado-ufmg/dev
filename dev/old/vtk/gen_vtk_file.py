import numpy as np
import vtk

if __name__ == '__main__':

    f = open('./dev/vtk/Example-2.csv')

    pd = vtk.vtkPolyData()
    points = vtk.vtkPoints()
    cells = vtk.vtkCellArray()
    stress = vtk.vtkFloatArray()
    stress.SetName('Stress')

    line = f.readline()
    for line in iter(lambda: f.readline(), ""):
        if 'None' in line: continue
        if 'Faces' in line: break
        v = line.split(',')
        points.InsertNextPoint(float(v[1]), float(v[2]), float(v[3]))
        stress.InsertNextTuple1(float(v[4]))
    
    for line in iter(lambda: f.readline(), ""):
        # v = line.split(',')
        # cell = vtk.vtkTriangle()
        # Ids = cell.GetPointIds()
        # for kId in range(len(v)):
        #     Ids.SetId(kId,int(v[kId]))
        # cells.InsertNextCell(cell)
        v = line.split(',')
        cell = vtk.vtkTriangle()
        cell.GetPointIds().SetId(0, int(v[0]))
        cell.GetPointIds().SetId(1, int(v[1]))
        cell.GetPointIds().SetId(2, int(v[2]))
        cells.InsertNextCell(cell)
    f.close()

    print(cells)

    pd.SetPoints(points)
    pd.SetPolys(cells)
    pd.GetPointData().AddArray(stress)

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName('Example-2.vtp')
    writer.SetInputData(pd)

    writer.Write()

    # # f = open('Example-2.csv')

    # vertices = np.array([
    #     [0, 0, 0],
    #     [1, 0, 0],
    #     [0, 1, 0],
    #     [0, 0, 1],
    # ], dtype=np.float32)

    # faces = np.array([
    #     [0, 1, 2],
    #     [0, 1, 3],
    # ], dtype=np.int32)

    # edges = np.array([
    #     [0, 1],
    #     [1, 2],
    #     [2, 0],
    #     [0, 3],
    #     [1, 3],
    # ], dtype=np.int32)

    # stresses = np.array([
    #     1,
    #     2,
    #     3,
    #     4,
    # ], dtype=np.float32)

    # pd = vtk.vtkPolyData()
    # points = vtk.vtkPoints()
    # cells = vtk.vtkCellArray()
    # connectivity = vtk.vtkIntArray()
    # # connectivity.SetName('Connectivity')
    # stress = vtk.vtkFloatArray()
    # stress.SetName('Stress')

    # for i in len(vertices[:, 0]):
    #     points.InsertNextPoint(vertices[i, 0], vertices[i, 1], vertices[i, 2])
    #     stress.InsertNextTuple1(stresses[i])

    # for line in iter(lambda: f.readline(), ""):
    #     v = line.split(',')
    #     cell = vtk.vtkTriangle()
    #     Ids = cell.GetPointIds()
    #     for kId in range(len(v)):
    #         Ids.SetId(kId,int(v[kId]))
    #     cells.InsertNextCell(cell)
    # # f.close()

    # pd.SetPoints(points)
    # pd.SetPolys(cells)
    # pd.GetPointData().AddArray(stress)
    # # pd.GetPointData().AddArray(connectivity)

    # writer = vtk.vtkXMLPolyDataWriter()
    # writer.SetFileName('Example-2.vtp')
    # writer.SetInputData(pd)

    # writer.Write()