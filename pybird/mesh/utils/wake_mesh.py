from math import acos
from typing import Callable, List
from numpy import argwhere, array, asarray, cross, deg2rad, dot, ndarray, allclose, zeros
from scipy.spatial.transform import Rotation as R

from pybird.geo.geo import Geo
from pybird.models.enums import TailShape
from pybird.models.types import Point
from pybird.geo.utils import bezier, vector, circle
from pybird.models.wake_model import WakeModel

def _find_vertices(vertices: ndarray, edges: ndarray, faces: ndarray, func: Callable[[Point], bool], pInitial: Point, pFinal: Point) -> List[ndarray]:
    """Find the vertices that are contained in the curve"""

    verticesOut = []
    edgesOut = []
    facesOut = []

    # Find the start point tag
    for i in range(len(vertices[:, 0])):
        if (allclose(vertices[i, :], pInitial)):
            pTag = i
            verticesOut.append(pTag)
            break

    # Loop until the last point is reached
    while True:

        # Find edges tags that contain pTag
        edgesTags = argwhere((edges[:, 0] == pTag) | (edges[:, 1] == pTag))
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
        facesTag = argwhere(((faces[:, 0] == pTag) | (faces[:, 1] == pTag) | (faces[:, 2] == pTag)) & ((faces[:, 0] == pNewTag) | (faces[:, 1] == pNewTag) | (faces[:, 2] == pNewTag)))

        # Save the vertices, edge and face
        pTag = pNewTag
        verticesOut.append(pNewTag)
        facesOut.append([facesTag[0, 0], facesTag[1, 0]] if facesTag[0, 0] < facesTag[1, 0] else [facesTag[1, 0], facesTag[0, 0]])

        # Check if new point is equal pFinal
        if allclose(vertices[pNewTag, :], pFinal):
            break

    return [asarray(verticesOut).astype(int), asarray(edgesOut).reshape(len(edgesOut)).astype(int), asarray(facesOut).astype(int)]

def _create_arrays(vertices: ndarray, edges: ndarray, faces: ndarray, facesTags: ndarray, edgesTags: ndarray, e1Base: Point, e2Base: Point, e1U: Point, e2U: Point, alignAll: bool = False) -> ndarray:
    
    def _wake_array(i: int) -> Point:

        # Edge array
        edge = edges[edgesTags[i], :]
        t = vertices[edge[1], :] - vertices[edge[0], :]

        # Normals
        face = faces[facesTags[i, 0], :]
        n1 = cross(vertices[face[1], :] - vertices[face[0], :], vertices[face[2], :] - vertices[face[0], :])

        face = faces[facesTags[i, 1], :]
        n2 = cross(vertices[face[1], :] - vertices[face[0], :], vertices[face[2], :] - vertices[face[0], :])

        # Wake array
        w1 = vector.unary(cross(t, n1))
        w2 = vector.unary(cross(t, n2))

        # Correct direction
        w1 = w1 if w1[0] <= 0 else -w1
        w2 = w2 if w2[0] <= 0 else -w2

        return vector.unary(w1 + w2)
    
    n = len(facesTags[:, 0])
    arrays = zeros((n + 1, 3))

    for i in range(n + 1):
        
        if i == 0 or i == n:
            w = _wake_array(0 if i == 0 else n - 1)
        else:
            w = vector.unary(_wake_array(i - 1) + _wake_array(i))
        
        if alignAll or i > n / 2:
            e1CompU = dot(w, e1U)
            e2CompU = dot(w, e2U)
            planeVecU = vector.unary(e1CompU * e1U + e2CompU * e2U)
            w = planeVecU
        else:
            e1CompB = dot(w, e1Base)
            e2CompB = dot(w, e2Base)
            planeVecB = vector.unary(e1CompB * e1Base + e2CompB * e2Base)

            e1CompU = dot(w, e1U)
            e2CompU = dot(w, e2U)
            planeVecU = vector.unary(e1CompU * e1U + e2CompU * e2U)

            factor = 2 * i / n
            w = vector.unary((factor ** 2) * planeVecU + (1 - factor ** 2) * planeVecB)
        
        arrays[i, :] = w

    return arrays

def _create_grid(vertices: ndarray, verticesTags: ndarray, verticesArray: ndarray, u: ndarray, nWake: int, ds: float, func: Callable[[float], float]) -> List[ndarray]:
    
    nSurf = len(verticesTags)
    gridTags = zeros((nSurf, nWake))
    gridVertices = zeros((int(nSurf * nWake), 3))

    count = 0
    for i in range(nSurf):

        for j in range(nWake):

            if j == 0:
                gridVertices[count, :] = vertices[verticesTags[i], :]
                gridTags[i, j] = count
            else:
                factor = func(ds * j)
                newVertice = gridVertices[count - 1] + ds * vector.unary(factor * verticesArray[i, :] + (1 - factor) * u)
                gridVertices[count, :] = newVertice
                gridTags[i, j] = count

            count += 1

    return [gridVertices, gridTags.astype(int)]

def build_wake(vertices: ndarray, edges: ndarray, faces: ndarray, geo: Geo, wake_dist: float, accom_dist: float, alpha: float, beta: float) -> List[WakeModel]:

    accom_dist = accom_dist if accom_dist is not None and accom_dist > 3 else 3
    a = 1e-2 * accom_dist / (1 - 1e-2)
    def func0(x: float) -> float:
        """Calculate the multiplication factor of the wake"""
        return a / (a + x)

    def func1(x: Point) -> bool:
        """Verify if Point x is contained in left wing trailing edge"""

        dist1 = bezier.point_is_contained([geo.p12e, geo.c20e, geo.c21e, geo.p13e], x)
        dist2 = bezier.point_is_contained([geo.p11e, geo.c18e, geo.c19e, geo.p12e], x)
        dist3 = bezier.point_is_contained([geo.p10e, geo.c16e, geo.c17e, geo.p11e], x)
        dist4 = bezier.point_is_contained([geo.p9e, geo.c14e, geo.c15e, geo.p10e], x)
        dist5 = bezier.point_is_contained([geo.p8e, geo.c12e, geo.c13e, geo.p9e], x)
        dist6 = bezier.point_is_contained([geo.p7e, geo.c10e, geo.c11e, geo.p8e], x)
        
        return min([dist1, dist2, dist3, dist4, dist5, dist6])
    
    def func2(x: Point) -> bool:
        """Verify if Point x is contained in right wing trailing edge"""

        dist1 = bezier.point_is_contained([geo.p12d, geo.c20d, geo.c21d, geo.p13d], x)
        dist2 = bezier.point_is_contained([geo.p11d, geo.c18d, geo.c19d, geo.p12d], x)
        dist3 = bezier.point_is_contained([geo.p10d, geo.c16d, geo.c17d, geo.p11d], x)
        dist4 = bezier.point_is_contained([geo.p9d, geo.c14d, geo.c15d, geo.p10d], x)
        dist5 = bezier.point_is_contained([geo.p8d, geo.c12d, geo.c13d, geo.p9d], x)
        dist6 = bezier.point_is_contained([geo.p7d, geo.c10d, geo.c11d, geo.p8d], x)
        
        return min([dist1, dist2, dist3, dist4, dist5, dist6])
    
    if geo.data.tail.shape == TailShape.rounded:
        circleCenter = circle.find_center(geo.p47e, geo.p51, geo.p47d)
        circleArc = circle.circle_arc(circleCenter, geo.p47e, geo.p47d, 200, removeEdges=False)
    
    def func3(x: Point) -> bool:
        """Verify if Point x is contained in tail trailing edge"""

        if geo.data.tail.shape == TailShape.rounded:
            dist = ((circleArc[:, 0] - x[0]) ** 2 + (circleArc[:, 1] - x[1]) ** 2 + (circleArc[:, 2] - x[2]) ** 2) ** 0.5
            return min(dist)
        else:
            d1 = vector.norm(cross(x - geo.p47e, x - geo.p51)) / vector.norm(geo.p51 - geo.p47e)
            d2 = vector.norm(cross(x - geo.p47d, x - geo.p51)) / vector.norm(geo.p51 - geo.p47d)
            return min([d1, d2])

    # Wake parameters
    wake_dist = wake_dist if wake_dist is not None else 30 * (geo.data.wing.h1 + geo.data.wing.h7)
    ds = 5e-2
    nWake = int(wake_dist / ds)

    # Freestream
    x = array([-1, 0, 0])
    y = array([0, -1, 0])
    z = array([0, 0, 1])

    if alpha is not None:
        r = R.from_rotvec(-deg2rad(alpha) * y)
        x = r.apply(x)
        z = r.apply(z)

    if beta is not None:
        r = R.from_rotvec(-deg2rad(beta) * z)
        x = r.apply(x)
        y = r.apply(y)

    # Wing
    n = vector.unary(cross(geo.p13e - geo.p1e, geo.p14e - geo.p1e))
    e1 = vector.unary(geo.p13e - geo.p1e)
    e2 = cross(n, e1)
    curve_vertices_tags, curve_edges_tags, curve_faces_tags = _find_vertices(vertices, edges, faces, func1, geo.p13e, geo.p7e)
    curve_array = _create_arrays(vertices, edges, faces, curve_faces_tags, curve_edges_tags, e1, e2, x, z, alignAll=False)
    grid_vertices_list, grid_vertices_tags = _create_grid(vertices, curve_vertices_tags, curve_array, x, nWake, ds, func0)
    
    leftWing = WakeModel(
        vertices=grid_vertices_list,
        grid=grid_vertices_tags,
        faces=curve_faces_tags,
    )

    n = vector.unary(cross(geo.p13d - geo.p1d, geo.p15d - geo.p1d))
    e1 = vector.unary(geo.p13d - geo.p1d)
    e2 = cross(n, e1)
    curve_vertices_tags, curve_edges_tags, curve_faces_tags = _find_vertices(vertices, edges, faces, func2, geo.p13d, geo.p7d)
    curve_array = _create_arrays(vertices, edges, faces, curve_faces_tags, curve_edges_tags, e1, e2, x, z, alignAll=False)
    grid_vertices_list, grid_vertices_tags = _create_grid(vertices, curve_vertices_tags, curve_array, x, nWake, ds, func0)
    
    rightWing = WakeModel(
        vertices=grid_vertices_list,
        grid=grid_vertices_tags,
        faces=curve_faces_tags,
    )

    # Wake
    curve_vertices_tags, curve_edges_tags, curve_faces_tags = _find_vertices(vertices, edges, faces, func3, geo.p47e, geo.p47d)
    curve_array = _create_arrays(vertices, edges, faces, curve_faces_tags, curve_edges_tags, None, None, x, z, alignAll=True)
    grid_vertices_list, grid_vertices_tags = _create_grid(vertices, curve_vertices_tags, curve_array, x, nWake, ds, func0)
    
    tail = WakeModel(
        vertices=grid_vertices_list,
        grid=grid_vertices_tags,
        faces=curve_faces_tags,
    )
    
    return [leftWing, rightWing, tail]