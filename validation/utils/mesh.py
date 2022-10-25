from dataclasses import dataclass
import sys
from typing import Callable, List
import gmsh
import numpy as np
from scipy.spatial.transform import Rotation as R

@dataclass
class MeshData:
    vertices: np.ndarray
    edges: np.ndarray
    faces: np.ndarray
    facesCenter: np.ndarray
    e1: np.ndarray
    e2: np.ndarray
    e3: np.ndarray
    controlPoints: np.ndarray
    facesAreas: np.ndarray
    facesMaxDistance: np.ndarray
    p1Local: np.ndarray
    p2Local: np.ndarray
    p3Local: np.ndarray
    wake_vertices: np.ndarray = None
    wake_grid: np.ndarray = None
    wake_faces: np.ndarray = None

def unary(a: np.ndarray) -> np.ndarray:
    """Returns a unary vector"""
    lenght = np.linalg.norm(a)
    if -1e-8 < lenght < 1e-8:
        return np.zeros(len(a))
    
    return a / lenght

def _find_vertices(vertices: np.ndarray, edges: np.ndarray, faces: np.ndarray, func: Callable[[np.ndarray], bool], pInitial: np.ndarray, pFinal: np.ndarray) -> List[np.ndarray]:
    """Find the vertices that are contained in the curve"""

    verticesOut = []
    edgesOut = []
    facesOut = []

    # Find the start point tag
    for i in range(len(vertices[:, 0])):
        if (np.allclose(vertices[i, :], pInitial)):
            pTag = i
            verticesOut.append(pTag)
            break

    # Loop until the last point is reached
    while True:

        # Find edges tags that contain pTag
        edgesTags = np.argwhere((edges[:, 0] == pTag) | (edges[:, 1] == pTag))
        edgesTags = edgesTags.reshape(len(edgesTags))

        # Find the correct edge
        d = {"edges": [], "edgesTag": [], "dist": []}
        for i in edgesTags:
            if i not in edgesOut:
                edge = edges[i, :]
                auxTag = edge[0] if edge[0] != pTag else edge[1]
                p = vertices[auxTag, :]
                d["edges"].append(edge)
                d["edgesTag"].append(i)
                d["dist"].append(func(p))

        index = d["dist"].index(min(d["dist"]))
        edge = d["edges"][index]
        edgesOut.append(d["edgesTag"][index])
        pNewTag = edge[0] if edge[0] != pTag else edge[1]

        # Find the two faces that contain the correct edge
        facesTag = np.argwhere(((faces[:, 0] == pTag) | (faces[:, 1] == pTag) | (faces[:, 2] == pTag)) & ((faces[:, 0] == pNewTag) | (faces[:, 1] == pNewTag) | (faces[:, 2] == pNewTag)))

        # Save the vertices, edge and face
        pTag = pNewTag
        verticesOut.append(pNewTag)
        facesOut.append([facesTag[0, 0], facesTag[1, 0]] if facesTag[0, 0] < facesTag[1, 0] else [facesTag[1, 0], facesTag[0, 0]])

        # Check if new point is equal pFinal
        if np.allclose(vertices[pNewTag, :], pFinal):
            break

    return [np.asarray(verticesOut).astype(int), np.asarray(edgesOut).reshape(len(edgesOut)).astype(int), np.asarray(facesOut).astype(int)]

def _create_arrays(vertices: np.ndarray, edges: np.ndarray, faces: np.ndarray, facesTags: np.ndarray, edgesTags: np.ndarray, e1Base: np.ndarray, e2Base: np.ndarray, e1U: np.ndarray, e2U: np.ndarray, alignAll: bool = False) -> np.ndarray:
    
    def _wake_array(i: int) -> np.ndarray:

        # Edge array
        edge = edges[edgesTags[i], :]
        t = vertices[edge[1], :] - vertices[edge[0], :]

        # Normals
        face = faces[facesTags[i, 0], :]
        n1 = np.cross(vertices[face[1], :] - vertices[face[0], :], vertices[face[2], :] - vertices[face[0], :])

        face = faces[facesTags[i, 1], :]
        n2 = np.cross(vertices[face[1], :] - vertices[face[0], :], vertices[face[2], :] - vertices[face[0], :])

        # Wake array
        w1 = unary(np.cross(t, n1))
        w2 = unary(np.cross(t, n2))

        # Correct direction
        w1 = w1 if w1[0] <= 0 else -w1
        w2 = w2 if w2[0] <= 0 else -w2

        return unary(w1 + w2)
    
    n = len(facesTags[:, 0])
    arrays = np.zeros((n + 1, 3))

    for i in range(n + 1):
        
        if i == 0 or i == n:
            w = _wake_array(0 if i == 0 else n - 1)
        else:
            w = unary(_wake_array(i - 1) + _wake_array(i))
        
        arrays[i, :] = w

    return arrays

def _create_grid_backup(vertices: np.ndarray, verticesTags: np.ndarray, verticesArray: np.ndarray, u: np.ndarray, nWake: int, ds: float, func: Callable[[float], float]) -> List[np.ndarray]:
    
    nSurf = len(verticesTags)
    gridTags = np.zeros((nSurf, nWake))
    gridVertices = np.zeros((int(nSurf * nWake), 3))

    count = 0
    for i in range(nSurf):

        for j in range(nWake):

            if j == 0:
                gridVertices[count, :] = vertices[verticesTags[i], :]
                gridTags[i, j] = count
            else:
                factor = func(ds * j)
                newVertice = gridVertices[count - 1] + ds * unary(factor * verticesArray[i, :] + (1 - factor) * u)
                gridVertices[count, :] = newVertice
                gridTags[i, j] = count

            count += 1

    return [gridVertices, gridTags.astype(int)]

def _create_grid(vertices: np.ndarray, verticesTags: np.ndarray, verticesArray: np.ndarray, u: np.ndarray, nWake: int, ds: float, func: Callable[[float], float], wake_lenght: float) -> List[np.ndarray]:
    
    nSurf = len(verticesTags)
    gridTags = np.zeros((nSurf, nWake))
    gridVertices = np.zeros((int(nSurf * nWake), 3))

    count = 0
    for i in range(nSurf):
        
        lenght = 0.0

        for j in range(nWake):

            if j == 0:
                gridVertices[count, :] = vertices[verticesTags[i], :]
                gridTags[i, j] = count
            elif j == nWake - 1:
                newVertice = gridVertices[count - 1] + (wake_lenght - lenght) * unary(u)
                gridVertices[count, :] = newVertice
                gridTags[i, j] = count
            else:
                factor = func(ds * j)
                newVertice = gridVertices[count - 1] + ds * unary(factor * verticesArray[i, :] + (1 - factor) * u)
                gridVertices[count, :] = newVertice
                gridTags[i, j] = count
            
            lenght = lenght + ds

            count += 1

    return [gridVertices, gridTags.astype(int)]

def process_foil(foil: np.ndarray, chord: float) -> List[np.ndarray]:

    foil = foil * chord

    foil[:, 0] = -foil[:, 0]
    foil[:, 0] = foil[:, 0] - 0.5 * (max(foil[:, 0]) + min(foil[:, 0]))
    index = np.argmax(foil[:, 0])
    return [
        foil[0, :],         # Trailing point
        foil[index, :],     # Leading edge point
        foil[1:index, :],    # Upper surface points
        foil[index + 1:-1, :],    # Lower surface points
    ]

def unary(a: np.ndarray):

    if a.ndim == 2:
        n = a.shape[0]
        for i in range(n):
            a[i, :] = a[i, :] / np.linalg.norm(a[i, :])
        return a

    return a / np.linalg.norm(a)

def gen_surface_mesh(foil: np.ndarray, span: float, chord: float, cell_size: float, le_ratio: float, te_ratio: float) -> List:

    lc = cell_size
    te, le, upper, lower = process_foil(foil, chord)
    n = 200
    coeffs = np.linspace(-1. + 2. / (n - 1), 1. - 2. / (n - 1), num=200)

    gmsh.initialize(sys.argv)
    gmsh.model.add("model")

    # leading and trailing edge points
    le_e_point = gmsh.model.geo.addPoint(le[0], le[1], -0.5 * span, lc)
    le_d_point = gmsh.model.geo.addPoint(le[0], le[1], 0.5 * span, lc)
    te_e_point = gmsh.model.geo.addPoint(te[0], te[1], -0.5 * span, lc)
    te_d_point = gmsh.model.geo.addPoint(te[0], te[1], 0.5 * span, lc)

    # Upper and lower points
    upper_e_points = []
    for i in range(upper.shape[0]):
        upper_e_points.append(gmsh.model.geo.addPoint(upper[i, 0], upper[i, 1], -0.5 * span, lc))
    
    lower_e_points = []
    for i in range(lower.shape[0]):
        lower_e_points.append(gmsh.model.geo.addPoint(lower[i, 0], lower[i, 1], -0.5 * span, lc))
    
    upper_d_points = []
    for i in range(upper.shape[0]):
        upper_d_points.append(gmsh.model.geo.addPoint(upper[i, 0], upper[i, 1], 0.5 * span, lc))
    
    lower_d_points = []
    for i in range(lower.shape[0]):
        lower_d_points.append(gmsh.model.geo.addPoint(lower[i, 0], lower[i, 1], 0.5 * span, lc))
    
    # Curves
    upper_e_curve = gmsh.model.geo.addPolyline([te_e_point] + [upper_e_points[i] for i in range(len(upper_e_points))] + [le_e_point])
    lower_e_curve = gmsh.model.geo.addPolyline([le_e_point] + [lower_e_points[i] for i in range(len(lower_e_points))] + [te_e_point])
    upper_d_curve = gmsh.model.geo.addPolyline([te_d_point] + [upper_d_points[i] for i in range(len(upper_d_points))] + [le_d_point])
    lower_d_curve = gmsh.model.geo.addPolyline([le_d_point] + [lower_d_points[i] for i in range(len(lower_d_points))] + [te_d_point])
    line_te = gmsh.model.geo.addPolyline([te_e_point] + [gmsh.model.geo.addPoint(te[0], te[1], 0.5 * span * coeff, lc) for coeff in coeffs] + [te_d_point])
    line_le = gmsh.model.geo.addPolyline([le_e_point] + [gmsh.model.geo.addPoint(le[0], le[1], 0.5 * span * coeff, lc) for coeff in coeffs] + [le_d_point])

    # Curve loop
    upper_cl = gmsh.model.geo.addCurveLoop([line_te, upper_d_curve, -line_le, -upper_e_curve])
    lower_cl = gmsh.model.geo.addCurveLoop([line_le, lower_d_curve, -line_te, -lower_e_curve])

    # Surfaces
    gmsh.model.geo.addSurfaceFilling([upper_cl])
    gmsh.model.geo.addSurfaceFilling([lower_cl])

    # Refinement
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", [line_le])
    gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", 200)

    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", le_ratio * lc)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", lc)
    gmsh.model.mesh.field.setNumber(2, "DistMin", 0.02 * chord)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.5 * chord)

    gmsh.model.mesh.field.add("Distance", 3)
    gmsh.model.mesh.field.setNumbers(3, "CurvesList", [line_te])
    gmsh.model.mesh.field.setNumber(3, "NumPointsPerCurve", 200)

    gmsh.model.mesh.field.add("Threshold", 4)
    gmsh.model.mesh.field.setNumber(4, "InField", 3)
    gmsh.model.mesh.field.setNumber(4, "SizeMin", te_ratio * lc)
    gmsh.model.mesh.field.setNumber(4, "SizeMax", lc)
    gmsh.model.mesh.field.setNumber(4, "DistMin", 0.01 * chord)
    gmsh.model.mesh.field.setNumber(4, "DistMax", 0.5 * chord)

    gmsh.model.mesh.field.add("Min", 5)
    gmsh.model.mesh.field.setNumbers(5, "FieldsList", [2, 4])

    gmsh.model.mesh.field.setAsBackgroundMesh(5)

    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

    gmsh.model.geo.synchronize()
    gmsh.option.setNumber("Mesh.Smoothing", 100)

    gmsh.option.setNumber("Mesh.Algorithm", 5)
    gmsh.model.mesh.generate(2)

    # if "-nopopup" not in sys.argv:
    #     gmsh.fltk.initialize()
    #     while gmsh.fltk.isAvailable():
    #         gmsh.fltk.wait()
    
    # Convert mesh
    verticesData = gmsh.model.mesh.getNodes()
    verticesArray = verticesData[1].reshape((len(verticesData[0]), 3))

    gmsh.model.mesh.createEdges()
    edgesData = gmsh.model.mesh.getAllEdges()
    edgesArray = edgesData[1].reshape((len(edgesData[0]), 2)).astype(int) - 1

    gmsh.model.mesh.createFaces()
    facesData = gmsh.model.mesh.getAllFaces(3)
    facesArray = facesData[1].reshape((len(facesData[0]), 3)).astype(int) - 1

    # Remove extra vertices
    indexes = -np.ones(verticesArray.shape[0])
    indexes[facesArray[:, 0]] = 1; indexes[facesArray[:, 1]] = 1; indexes[facesArray[:, 2]] = 1

    count = 0
    for i in range(verticesArray.shape[0]):
        if indexes[i] == 1:
            indexes[i] = count
            count += 1

    for i in range(facesArray.shape[0]):
        facesArray[i, 0] = indexes[facesArray[i, 0]]
        facesArray[i, 1] = indexes[facesArray[i, 1]]
        facesArray[i, 2] = indexes[facesArray[i, 2]]
        
    for i in range(edgesArray.shape[0]):
        edgesArray[i, 0] = indexes[edgesArray[i, 0]]
        edgesArray[i, 1] = indexes[edgesArray[i, 1]]
        
    removeIndexes = np.argwhere(indexes == -1).reshape(-1)
        
    vertices = np.delete(verticesArray, removeIndexes, axis=0)
    edges = edgesArray
    faces = facesArray

    # Faces center
    facesCenter = (1 / 3) * (vertices[faces[:, 0], :] + vertices[faces[:, 1], :] + vertices[faces[:, 2], :])

    # Base vectors
    e3 = unary(np.cross(vertices[faces[:, 1], :] - vertices[faces[:, 0], :], vertices[faces[:, 2], :] - vertices[faces[:, 0], :]))
    e1 = unary(vertices[faces[:, 1], :] - facesCenter)
    e2 = unary(np.cross(e3, e1))
        
    # Control points
    controlPoints = facesCenter + 1e-8 * e3

    # Faces areas
    auxVec = np.cross(vertices[faces[:, 1]] - vertices[faces[:, 0]], vertices[faces[:, 2]] - vertices[faces[:, 0]])
    auxNorm = (auxVec[:, 0] ** 2 + auxVec[:, 1] ** 2 + auxVec[:, 2] ** 2) ** 0.5
    facesAreas = 0.5 * auxNorm

    # Faces max distance
    facesMaxDistance = 10 * (4 * facesAreas / np.pi) ** 0.5

    p1 = vertices[faces[:, 0], :] - facesCenter
    p2 = vertices[faces[:, 1], :] - facesCenter
    p3 = vertices[faces[:, 2], :] - facesCenter

    n = faces.shape[0]

    p1Local = np.empty((n, 2), dtype=np.double)
    p2Local = np.empty((n, 2), dtype=np.double)
    p3Local = np.empty((n, 2), dtype=np.double)

    p1Local[:, 0], p1Local[:, 1] = (p1 * e1).sum(axis=1), (p1 * e2).sum(axis=1)
    p2Local[:, 0], p2Local[:, 1] = (p2 * e1).sum(axis=1), (p2 * e2).sum(axis=1)
    p3Local[:, 0], p3Local[:, 1] = (p3 * e1).sum(axis=1), (p3 * e2).sum(axis=1)

    # Store data
    data = MeshData(
        vertices=vertices,
        edges=edges,
        faces=faces,
        facesCenter=facesCenter,
        e1=e1,
        e2=e2,
        e3=e3,
        controlPoints=controlPoints,
        facesAreas=facesAreas,
        facesMaxDistance=facesMaxDistance,
        p1Local=p1Local,
        p2Local=p2Local,
        p3Local=p3Local,
    )

    return [data, te]

def gen_mesh(foil: np.ndarray, span: float, chord: float, cell_size: float, le_ratio: float, te_ratio: float, alpha: float, wake_lenght: float, wake_accom_dist: float) -> MeshData:

    # Surface mesh
    data: MeshData = None
    data, te_point = gen_surface_mesh(foil, span, chord, cell_size, le_ratio, te_ratio)

    # Wake mesh
    ds = 5e-2
    nWake = int(wake_accom_dist / ds) + 10

    e1 = np.array([1, 0, 0])
    e2 = np.array([0, 1, 0])
    e3 = np.array([0, 0, 1])

    x = np.array([-1, 0, 0])
    y = np.array([0, -1, 0])
    z = np.array([0, 0, 1])

    if alpha is not None:
        r = R.from_rotvec(-np.deg2rad(alpha) * z)
        x = r.apply(x)
        z = r.apply(z)

    p_initial = np.array([te_point[0], te_point[1], - 0.5 * span])
    p_final = np.array([te_point[0], te_point[1], 0.5 * span])
    vec_unary = (p_final - p_initial) / np.linalg.norm(p_final - p_initial)

    def contain_point(point: np.ndarray) -> bool:
        """Verify if the point is contained by the trailing edge"""
        vec = (point - p_initial) / np.linalg.norm(point - p_initial)
        angle = np.arccos(np.dot(vec, vec_unary))
        return angle
    
    a = 1e-1 * wake_accom_dist / (1 - 1e-1)
    def func0(x: float) -> float:
        """Calculate the multiplication factor of the wake"""
        return a / (a + x)

    vertices_tags, edges_tags, faces_tags = _find_vertices(data.vertices, data.edges, data.faces, contain_point, p_initial, p_final)
    curve_array = _create_arrays(data.vertices, data.edges, data.faces, faces_tags, edges_tags, e1, e3, x, y)    
    grid_vertices_list, grid_vertices_tags = _create_grid(data.vertices, vertices_tags, curve_array, x, nWake, ds, func0, wake_lenght)

    data.wake_vertices = grid_vertices_list
    data.wake_grid = grid_vertices_tags
    data.wake_faces = faces_tags

    return data

def __test_process_foil() -> None:

    import matplotlib.pyplot as plt

    foil = np.loadtxt('./foils/NACA0012.dat')
    chord = 0.2
    te, le, upper, lower = process_foil(foil, chord)

    plt.figure()
    plt.scatter(te[0], te[1], label='trailing edge')
    plt.scatter(le[0], le[1], label='leading edge')
    plt.scatter(upper[:, 0], upper[:, 1], label='upper')
    plt.scatter(lower[:, 0], lower[:, 1], label='lower')
    plt.grid()
    plt.axis('equal')
    plt.legend()
    plt.show()

    return

def __test_gen_mesh() -> None:

    foil = np.loadtxt('./foils/NACA0012.dat')
    chord = 0.2
    span = 0.5

    gen_surface_mesh(foil, span, chord)

    return

if __name__ == '__main__':
    __test_gen_mesh()