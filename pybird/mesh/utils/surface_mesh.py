import sys
from typing import List
import gmsh
from numpy import argwhere, isin, linspace, ndarray

from pybird.geo.geo import Geo
from pybird.geo.utils import bezier

class SurfaceMesh:

    def __init__(self, geo: Geo,
                       size: float = None,
                       n_wing_le: float = None,
                       n_wing_te: float = None,
                       n_head: float = None,
                       n_tail_le: float = None,
                       n_tail_te: float = None) -> None:
        
        size = size if size is not None else (geo.data.wing.h1 + geo.data.wing.h7) / 10
        n_wing_le = n_wing_le if n_wing_le is not None else 3
        n_wing_te = n_wing_te if n_wing_te is not None else 2
        n_head = n_head if n_head is not None else 3
        n_tail_le = n_tail_le if n_tail_le is not None else 3
        n_tail_te = n_tail_te if n_tail_te is not None else 2

        self.geo = geo
        self.cellSize = size
        self.wingLEcellSize = self.cellSize / n_wing_le
        self.wingTEcellSize = self.cellSize / n_wing_te
        self.headCellSize = self.cellSize / n_head
        self.tailLEcellSize = self.cellSize / n_tail_le
        self.tailTEcellSize = self.cellSize / n_tail_te
        return
    
    def build(self) -> None:

        gmsh.initialize(sys.argv)
        gmsh.model.add("model")
        gmsh.option.setNumber('General.Verbosity', 1)

        self.lc = self.cellSize

        self._wing()
        self._body()
        self._head()
        self._tail()

        gmsh.model.geo.synchronize()
        gmsh.option.setNumber("Mesh.Smoothing", 100)

        self._refine()

        gmsh.option.setNumber("Mesh.Algorithm", 5)
        gmsh.model.mesh.generate(2)

        # if "-nopopup" not in sys.argv:
        #     gmsh.fltk.initialize()
        #     while gmsh.fltk.isAvailable():
        #         gmsh.fltk.wait()

        # Convert mesh
        vertices = gmsh.model.mesh.getNodes()
        verticesArray = vertices[1].reshape((len(vertices[0]), 3))

        gmsh.model.mesh.createEdges()
        edges = gmsh.model.mesh.getAllEdges()
        edgesArray = edges[1].reshape((len(edges[0]), 2)).astype(int) - 1

        gmsh.model.mesh.createFaces()
        faces = gmsh.model.mesh.getAllFaces(3)
        facesArray = faces[1].reshape((len(faces[0]), 3)).astype(int) - 1

        # Create physical group

        # Left wing
        self.physicalGroup1 = gmsh.model.addPhysicalGroup( # Wing first section
            2,
            [
                self.surface1e, self.surface2e,
                self.surface3e, self.surface4e,
                self.surface13e, self.surface14e,
                self.surface15e, self.surface16e,
            ]
        )
        self.physicalGroup2 = gmsh.model.addPhysicalGroup( # Wing middle section
            2,
            [
                self.surface5e, self.surface6e,
                self.surface17e, self.surface18e,
            ]
        )
        self.physicalGroup3 = gmsh.model.addPhysicalGroup( # Wing tip
            2,
            [
                self.surface7e, self.surface8e,
                self.surface9e, self.surface10e,
                self.surface11e, self.surface12e,
                self.surface19e, self.surface20e,
                self.surface21e, self.surface22e,
                self.surface23e, self.surface24e,
            ]
        )

        # Left wing
        self.physicalGroup4 = gmsh.model.addPhysicalGroup( # Wing first section
            2,
            [
                self.surface1d, self.surface2d,
                self.surface3d, self.surface4d,
                self.surface13d, self.surface14d,
                self.surface15d, self.surface16d,
            ]
        )
        self.physicalGroup5 = gmsh.model.addPhysicalGroup( # Wing middle section
            2,
            [
                self.surface5d, self.surface6d,
                self.surface17d, self.surface18d,
            ]
        )
        self.physicalGroup6 = gmsh.model.addPhysicalGroup( # Wing tip
            2,
            [
                self.surface7d, self.surface8d,
                self.surface9d, self.surface10d,
                self.surface11d, self.surface12d,
                self.surface19d, self.surface20d,
                self.surface21d, self.surface22d,
                self.surface23d, self.surface24d,
            ]
        )

        # Body
        self.physicalGroup7 = gmsh.model.addPhysicalGroup(
            2,
            [
                self.surface25e, self.surface26e, self.surface27e, self.surface28e, self.surface29e, self.surface30e, self.surface31e, self.surface32e,
                self.surface25d, self.surface26d, self.surface27d, self.surface28d, self.surface29d, self.surface30d, self.surface31d, self.surface32d,
            ]
        )

        # Head
        self.physicalGroup8 = gmsh.model.addPhysicalGroup(
            2,
            [
                self.surface33, self.surface34, self.surface35, self.surface36, self.surface37, self.surface38, self.surface39, self.surface40,
            ]
        )

        # Til
        self.physicalGroup9 = gmsh.model.addPhysicalGroup(
            2,
            [
                self.surface41e, self.surface42e, self.surface43e, self.surface44e, self.surface45e, self.surface46e, self.surface47e, self.surface48e,
                self.surface41d, self.surface42d, self.surface43d, self.surface44d, self.surface45d, self.surface46d, self.surface47d, self.surface48d,
            ]
        )

        physicalGroupVertices1 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup1)[0] - 1
        physicalGroupVertices2 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup2)[0] - 1
        physicalGroupVertices3 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup3)[0] - 1
        leftWingFirstSectionFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices1)
        leftWingSecondSectionFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices2)
        leftWingThirdSectionFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices3)

        physicalGroupVertices4 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup4)[0] - 1
        physicalGroupVertices5 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup5)[0] - 1
        physicalGroupVertices6 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup6)[0] - 1
        rightWingFirstSectionFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices4)
        rightWingSecondSectionFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices5)
        rightWingThirdSectionFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices6)

        physicalGroupVertices7 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup7)[0] - 1
        bodyFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices7)

        physicalGroupVertices8 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup8)[0] - 1
        headFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices8)

        physicalGroupVertices9 = gmsh.model.mesh.getNodesForPhysicalGroup(2, self.physicalGroup9)[0] - 1
        tailFacesTags = self._get_physical_group_faces(facesArray, physicalGroupVertices9)

        gmsh.finalize()

        return verticesArray, edgesArray, facesArray, leftWingFirstSectionFacesTags, leftWingSecondSectionFacesTags, leftWingThirdSectionFacesTags, rightWingFirstSectionFacesTags, rightWingSecondSectionFacesTags, rightWingThirdSectionFacesTags, bodyFacesTags, headFacesTags, tailFacesTags
    
    def _get_physical_group_faces(self, faces: ndarray, vertices: ndarray) -> ndarray:
        
        facesList = argwhere(isin(faces[:, 0], vertices) & isin(faces[:, 1], vertices) & isin(faces[:, 2], vertices))
        facesList = facesList.reshape(len(facesList))

        return facesList
    
    def _wing(self) -> None:

        # Surface points
        self.p1e = gmsh.model.geo.addPoint(self.geo.p1e[0], self.geo.p1e[1], self.geo.p1e[2], self.lc)
        self.p2e = gmsh.model.geo.addPoint(self.geo.p2e[0], self.geo.p2e[1], self.geo.p2e[2], self.lc)
        self.p3e = gmsh.model.geo.addPoint(self.geo.p3e[0], self.geo.p3e[1], self.geo.p3e[2], self.lc)
        self.p4e = gmsh.model.geo.addPoint(self.geo.p4e[0], self.geo.p4e[1], self.geo.p4e[2], self.lc)
        self.p5e = gmsh.model.geo.addPoint(self.geo.p5e[0], self.geo.p5e[1], self.geo.p5e[2], self.lc)
        self.p6e = gmsh.model.geo.addPoint(self.geo.p6e[0], self.geo.p6e[1], self.geo.p6e[2], self.lc)
        self.p7e = gmsh.model.geo.addPoint(self.geo.p7e[0], self.geo.p7e[1], self.geo.p7e[2], self.lc)
        self.p8e = gmsh.model.geo.addPoint(self.geo.p8e[0], self.geo.p8e[1], self.geo.p8e[2], self.lc)
        self.p9e = gmsh.model.geo.addPoint(self.geo.p9e[0], self.geo.p9e[1], self.geo.p9e[2], self.lc)
        self.p10e = gmsh.model.geo.addPoint(self.geo.p10e[0], self.geo.p10e[1], self.geo.p10e[2], self.lc)
        self.p11e = gmsh.model.geo.addPoint(self.geo.p11e[0], self.geo.p11e[1], self.geo.p11e[2], self.lc)
        self.p12e = gmsh.model.geo.addPoint(self.geo.p12e[0], self.geo.p12e[1], self.geo.p12e[2], self.lc)
        self.p13e = gmsh.model.geo.addPoint(self.geo.p13e[0], self.geo.p13e[1], self.geo.p13e[2], self.lc)
        self.p14e = gmsh.model.geo.addPoint(self.geo.p14e[0], self.geo.p14e[1], self.geo.p14e[2], self.lc)
        self.p15e = gmsh.model.geo.addPoint(self.geo.p15e[0], self.geo.p15e[1], self.geo.p15e[2], self.lc)
        self.p16e = gmsh.model.geo.addPoint(self.geo.p16e[0], self.geo.p16e[1], self.geo.p16e[2], self.lc)
        self.p17e = gmsh.model.geo.addPoint(self.geo.p17e[0], self.geo.p17e[1], self.geo.p17e[2], self.lc)
        self.p18e = gmsh.model.geo.addPoint(self.geo.p18e[0], self.geo.p18e[1], self.geo.p18e[2], self.lc)
        self.p19e = gmsh.model.geo.addPoint(self.geo.p19e[0], self.geo.p19e[1], self.geo.p19e[2], self.lc)
        self.p20e = gmsh.model.geo.addPoint(self.geo.p20e[0], self.geo.p20e[1], self.geo.p20e[2], self.lc)
        self.p21e = gmsh.model.geo.addPoint(self.geo.p21e[0], self.geo.p21e[1], self.geo.p21e[2], self.lc)
        self.p22e = gmsh.model.geo.addPoint(self.geo.p22e[0], self.geo.p22e[1], self.geo.p22e[2], self.lc)
        self.p23e = gmsh.model.geo.addPoint(self.geo.p23e[0], self.geo.p23e[1], self.geo.p23e[2], self.lc)
        self.p24e = gmsh.model.geo.addPoint(self.geo.p24e[0], self.geo.p24e[1], self.geo.p24e[2], self.lc)
        self.p25e = gmsh.model.geo.addPoint(self.geo.p25e[0], self.geo.p25e[1], self.geo.p25e[2], self.lc)

        # Curves points
        self.curve13ePoints = []
        for i in range(len(self.geo.curve13e[:, 0])):
            self.curve13ePoints.append(gmsh.model.geo.addPoint(self.geo.curve13e[i, 0], self.geo.curve13e[i, 1], self.geo.curve13e[i, 2], self.lc))

        self.curve14ePoints = []
        for i in range(len(self.geo.curve14e[:, 0])):
            self.curve14ePoints.append(gmsh.model.geo.addPoint(self.geo.curve14e[i, 0], self.geo.curve14e[i, 1], self.geo.curve14e[i, 2], self.lc))
        
        self.curve15ePoints = []
        for i in range(len(self.geo.curve15e[:, 0])):
            self.curve15ePoints.append(gmsh.model.geo.addPoint(self.geo.curve15e[i, 0], self.geo.curve15e[i, 1], self.geo.curve15e[i, 2], self.lc))
        
        self.curve16ePoints = []
        for i in range(len(self.geo.curve16e[:, 0])):
            self.curve16ePoints.append(gmsh.model.geo.addPoint(self.geo.curve16e[i, 0], self.geo.curve16e[i, 1], self.geo.curve16e[i, 2], self.lc))
        
        self.curve17ePoints = []
        for i in range(len(self.geo.curve17e[:, 0])):
            self.curve17ePoints.append(gmsh.model.geo.addPoint(self.geo.curve17e[i, 0], self.geo.curve17e[i, 1], self.geo.curve17e[i, 2], self.lc))
        
        self.curve18ePoints = []
        for i in range(len(self.geo.curve18e[:, 0])):
            self.curve18ePoints.append(gmsh.model.geo.addPoint(self.geo.curve18e[i, 0], self.geo.curve18e[i, 1], self.geo.curve18e[i, 2], self.lc))
        
        self.curve19ePoints = []
        for i in range(len(self.geo.curve19e[:, 0])):
            self.curve19ePoints.append(gmsh.model.geo.addPoint(self.geo.curve19e[i, 0], self.geo.curve19e[i, 1], self.geo.curve19e[i, 2], self.lc))
        
        self.curve20ePoints = []
        for i in range(len(self.geo.curve20e[:, 0])):
            self.curve20ePoints.append(gmsh.model.geo.addPoint(self.geo.curve20e[i, 0], self.geo.curve20e[i, 1], self.geo.curve20e[i, 2], self.lc))
        
        self.curve21ePoints = []
        for i in range(len(self.geo.curve21e[:, 0])):
            self.curve21ePoints.append(gmsh.model.geo.addPoint(self.geo.curve21e[i, 0], self.geo.curve21e[i, 1], self.geo.curve21e[i, 2], self.lc))
        
        self.curve22ePoints = []
        for i in range(len(self.geo.curve22e[:, 0])):
            self.curve22ePoints.append(gmsh.model.geo.addPoint(self.geo.curve22e[i, 0], self.geo.curve22e[i, 1], self.geo.curve22e[i, 2], self.lc))
        
        self.curve23ePoints = []
        for i in range(len(self.geo.curve23e[:, 0])):
            self.curve23ePoints.append(gmsh.model.geo.addPoint(self.geo.curve23e[i, 0], self.geo.curve23e[i, 1], self.geo.curve23e[i, 2], self.lc))
        
        self.curve24ePoints = []
        for i in range(len(self.geo.curve24e[:, 0])):
            self.curve24ePoints.append(gmsh.model.geo.addPoint(self.geo.curve24e[i, 0], self.geo.curve24e[i, 1], self.geo.curve24e[i, 2], self.lc))
        
        self.curve25ePoints = []
        for i in range(len(self.geo.curve25e[:, 0])):
            self.curve25ePoints.append(gmsh.model.geo.addPoint(self.geo.curve25e[i, 0], self.geo.curve25e[i, 1], self.geo.curve25e[i, 2], self.lc))
        
        self.curve26ePoints = []
        for i in range(len(self.geo.curve26e[:, 0])):
            self.curve26ePoints.append(gmsh.model.geo.addPoint(self.geo.curve26e[i, 0], self.geo.curve26e[i, 1], self.geo.curve26e[i, 2], self.lc))
        
        self.curve27ePoints = []
        for i in range(len(self.geo.curve27e[:, 0])):
            self.curve27ePoints.append(gmsh.model.geo.addPoint(self.geo.curve27e[i, 0], self.geo.curve27e[i, 1], self.geo.curve27e[i, 2], self.lc))
        
        self.curve28ePoints = []
        for i in range(len(self.geo.curve28e[:, 0])):
            self.curve28ePoints.append(gmsh.model.geo.addPoint(self.geo.curve28e[i, 0], self.geo.curve28e[i, 1], self.geo.curve28e[i, 2], self.lc))
        
        self.curve29ePoints = []
        for i in range(len(self.geo.curve29e[:, 0])):
            self.curve29ePoints.append(gmsh.model.geo.addPoint(self.geo.curve29e[i, 0], self.geo.curve29e[i, 1], self.geo.curve29e[i, 2], self.lc))
        
        self.curve30ePoints = []
        for i in range(len(self.geo.curve30e[:, 0])):
            self.curve30ePoints.append(gmsh.model.geo.addPoint(self.geo.curve30e[i, 0], self.geo.curve30e[i, 1], self.geo.curve30e[i, 2], self.lc))
        
        self.curve31ePoints = []
        for i in range(len(self.geo.curve31e[:, 0])):
            self.curve31ePoints.append(gmsh.model.geo.addPoint(self.geo.curve31e[i, 0], self.geo.curve31e[i, 1], self.geo.curve31e[i, 2], self.lc))
        
        self.curve32ePoints = []
        for i in range(len(self.geo.curve32e[:, 0])):
            self.curve32ePoints.append(gmsh.model.geo.addPoint(self.geo.curve32e[i, 0], self.geo.curve32e[i, 1], self.geo.curve32e[i, 2], self.lc))
        
        self.curve33ePoints = []
        for i in range(len(self.geo.curve33e[:, 0])):
            self.curve33ePoints.append(gmsh.model.geo.addPoint(self.geo.curve33e[i, 0], self.geo.curve33e[i, 1], self.geo.curve33e[i, 2], self.lc))
        
        self.curve34ePoints = []
        for i in range(len(self.geo.curve34e[:, 0])):
            self.curve34ePoints.append(gmsh.model.geo.addPoint(self.geo.curve34e[i, 0], self.geo.curve34e[i, 1], self.geo.curve34e[i, 2], self.lc))
        
        self.curve35ePoints = []
        for i in range(len(self.geo.curve35e[:, 0])):
            self.curve35ePoints.append(gmsh.model.geo.addPoint(self.geo.curve35e[i, 0], self.geo.curve35e[i, 1], self.geo.curve35e[i, 2], self.lc))
        
        self.curve36ePoints = []
        for i in range(len(self.geo.curve36e[:, 0])):
            self.curve36ePoints.append(gmsh.model.geo.addPoint(self.geo.curve36e[i, 0], self.geo.curve36e[i, 1], self.geo.curve36e[i, 2], self.lc))
        
        self.curve37ePoints = []
        for i in range(len(self.geo.curve37e[:, 0])):
            self.curve37ePoints.append(gmsh.model.geo.addPoint(self.geo.curve37e[i, 0], self.geo.curve37e[i, 1], self.geo.curve37e[i, 2], self.lc))
        
        self.curve38ePoints = []
        for i in range(len(self.geo.curve38e[:, 0])):
            self.curve38ePoints.append(gmsh.model.geo.addPoint(self.geo.curve38e[i, 0], self.geo.curve38e[i, 1], self.geo.curve38e[i, 2], self.lc))
        
        self.curve39ePoints = []
        for i in range(len(self.geo.curve39e[:, 0])):
            self.curve39ePoints.append(gmsh.model.geo.addPoint(self.geo.curve39e[i, 0], self.geo.curve39e[i, 1], self.geo.curve39e[i, 2], self.lc))
        
        self.curve40ePoints = []
        for i in range(len(self.geo.curve40e[:, 0])):
            self.curve40ePoints.append(gmsh.model.geo.addPoint(self.geo.curve40e[i, 0], self.geo.curve40e[i, 1], self.geo.curve40e[i, 2], self.lc))
        
        self.curve41ePoints = []
        for i in range(len(self.geo.curve41e[:, 0])):
            self.curve41ePoints.append(gmsh.model.geo.addPoint(self.geo.curve41e[i, 0], self.geo.curve41e[i, 1], self.geo.curve41e[i, 2], self.lc))
        
        self.curve42ePoints = []
        for i in range(len(self.geo.curve42e[:, 0])):
            self.curve42ePoints.append(gmsh.model.geo.addPoint(self.geo.curve42e[i, 0], self.geo.curve42e[i, 1], self.geo.curve42e[i, 2], self.lc))
        
        self.curve43ePoints = []
        for i in range(len(self.geo.curve43e[:, 0])):
            self.curve43ePoints.append(gmsh.model.geo.addPoint(self.geo.curve43e[i, 0], self.geo.curve43e[i, 1], self.geo.curve43e[i, 2], self.lc))
        
        self.curve44ePoints = []
        for i in range(len(self.geo.curve44e[:, 0])):
            self.curve44ePoints.append(gmsh.model.geo.addPoint(self.geo.curve44e[i, 0], self.geo.curve44e[i, 1], self.geo.curve44e[i, 2], self.lc))
        
        self.curve45ePoints = []
        for i in range(len(self.geo.curve45e[:, 0])):
            self.curve45ePoints.append(gmsh.model.geo.addPoint(self.geo.curve45e[i, 0], self.geo.curve45e[i, 1], self.geo.curve45e[i, 2], self.lc))
        
        self.curve46ePoints = []
        for i in range(len(self.geo.curve46e[:, 0])):
            self.curve46ePoints.append(gmsh.model.geo.addPoint(self.geo.curve46e[i, 0], self.geo.curve46e[i, 1], self.geo.curve46e[i, 2], self.lc))
        
        self.curve47ePoints = []
        for i in range(len(self.geo.curve47e[:, 0])):
            self.curve47ePoints.append(gmsh.model.geo.addPoint(self.geo.curve47e[i, 0], self.geo.curve47e[i, 1], self.geo.curve47e[i, 2], self.lc))
        
        self.curve48ePoints = []
        for i in range(len(self.geo.curve48e[:, 0])):
            self.curve48ePoints.append(gmsh.model.geo.addPoint(self.geo.curve48e[i, 0], self.geo.curve48e[i, 1], self.geo.curve48e[i, 2], self.lc))

        # Curves
        self.curve1e = gmsh.model.geo.addPolyline([self.p1e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p1e, self.geo.c1e, self.geo.p2e], linspace(0, 1, num=50))[1:49]] + [self.p2e])
        self.curve2e = gmsh.model.geo.addPolyline([self.p2e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p2e, self.geo.c2e, self.geo.c3e, self.geo.p3e], linspace(0, 1, num=50))[1:49]] + [self.p3e])
        self.curve3e = gmsh.model.geo.addPolyline([self.p3e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p3e, self.geo.c4e, self.geo.c5e, self.geo.p4e], linspace(0, 1, num=50))[1:49]] + [self.p4e])
        self.curve4e = gmsh.model.geo.addLine(self.p4e, self.p5e)
        self.curve5e = gmsh.model.geo.addPolyline([self.p5e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p5e, self.geo.c6e, self.geo.c7e, self.geo.p6e], linspace(0, 1, num=50))[1:49]] + [self.p6e])
        self.curve6e = gmsh.model.geo.addPolyline([self.p6e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p6e, self.geo.c8e, self.geo.c9e, self.geo.p7e], linspace(0, 1, num=50))[1:49]] + [self.p7e])
        self.curve7e = gmsh.model.geo.addPolyline([self.p7e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p7e, self.geo.c10e, self.geo.c11e, self.geo.p8e], linspace(0, 1, num=50))[1:49]] + [self.p8e])
        self.curve8e = gmsh.model.geo.addPolyline([self.p8e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p8e, self.geo.c12e, self.geo.c13e, self.geo.p9e], linspace(0, 1, num=50))[1:49]] + [self.p9e])
        self.curve9e = gmsh.model.geo.addPolyline([self.p9e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p9e, self.geo.c14e, self.geo.c15e, self.geo.p10e], linspace(0, 1, num=50))[1:49]] + [self.p10e])
        self.curve10e = gmsh.model.geo.addPolyline([self.p10e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p10e, self.geo.c16e, self.geo.c17e, self.geo.p11e], linspace(0, 1, num=50))[1:49]] + [self.p11e])
        self.curve11e = gmsh.model.geo.addPolyline([self.p11e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p11e, self.geo.c18e, self.geo.c19e, self.geo.p12e], linspace(0, 1, num=50))[1:49]] + [self.p12e])
        self.curve12e = gmsh.model.geo.addPolyline([self.p12e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p12e, self.geo.c20e, self.geo.c21e, self.geo.p13e], linspace(0, 1, num=50))[1:49]] + [self.p13e])

        self.curve13e = gmsh.model.geo.addPolyline([self.p14e] + self.curve13ePoints + [self.p1e])
        self.curve14e = gmsh.model.geo.addPolyline([self.p15e] + self.curve14ePoints + [self.p1e])
        self.curve15e = gmsh.model.geo.addPolyline([self.p13e] + self.curve15ePoints + [self.p14e])
        self.curve16e = gmsh.model.geo.addPolyline([self.p13e] + self.curve16ePoints + [self.p15e])
        self.curve17e = gmsh.model.geo.addPolyline([self.p14e] + self.curve17ePoints + [self.p16e])
        self.curve18e = gmsh.model.geo.addPolyline([self.p15e] + self.curve18ePoints + [self.p17e])
        self.curve19e = gmsh.model.geo.addPolyline([self.p16e] + self.curve19ePoints + [self.p2e])
        self.curve20e = gmsh.model.geo.addPolyline([self.p17e] + self.curve20ePoints + [self.p2e])
        self.curve21e = gmsh.model.geo.addPolyline([self.p12e] + self.curve21ePoints + [self.p16e])
        self.curve22e = gmsh.model.geo.addPolyline([self.p12e] + self.curve22ePoints + [self.p17e])
        self.curve23e = gmsh.model.geo.addPolyline([self.p16e] + self.curve23ePoints + [self.p18e])
        self.curve24e = gmsh.model.geo.addPolyline([self.p17e] + self.curve24ePoints + [self.p19e])
        self.curve25e = gmsh.model.geo.addPolyline([self.p18e] + self.curve25ePoints + [self.p3e])
        self.curve26e = gmsh.model.geo.addPolyline([self.p19e] + self.curve26ePoints + [self.p3e])
        self.curve27e = gmsh.model.geo.addPolyline([self.p11e] + self.curve27ePoints + [self.p18e])
        self.curve28e = gmsh.model.geo.addPolyline([self.p11e] + self.curve28ePoints + [self.p19e])
        self.curve29e = gmsh.model.geo.addPolyline([self.p18e] + self.curve29ePoints + [self.p20e])
        self.curve30e = gmsh.model.geo.addPolyline([self.p19e] + self.curve30ePoints + [self.p21e])
        self.curve31e = gmsh.model.geo.addPolyline([self.p20e] + self.curve31ePoints + [self.p4e])
        self.curve32e = gmsh.model.geo.addPolyline([self.p21e] + self.curve32ePoints + [self.p4e])
        self.curve33e = gmsh.model.geo.addPolyline([self.p10e] + self.curve33ePoints + [self.p20e])
        self.curve34e = gmsh.model.geo.addPolyline([self.p10e] + self.curve34ePoints + [self.p21e])
        self.curve35e = gmsh.model.geo.addPolyline([self.p20e] + self.curve35ePoints + [self.p22e])
        self.curve36e = gmsh.model.geo.addPolyline([self.p21e] + self.curve36ePoints + [self.p23e])
        self.curve37e = gmsh.model.geo.addPolyline([self.p22e] + self.curve37ePoints + [self.p5e])
        self.curve38e = gmsh.model.geo.addPolyline([self.p23e] + self.curve38ePoints + [self.p5e])
        self.curve39e = gmsh.model.geo.addPolyline([self.p9e] + self.curve39ePoints + [self.p22e])
        self.curve40e = gmsh.model.geo.addPolyline([self.p9e] + self.curve40ePoints + [self.p23e])
        self.curve41e = gmsh.model.geo.addPolyline([self.p22e] + self.curve41ePoints + [self.p24e])
        self.curve42e = gmsh.model.geo.addPolyline([self.p23e] + self.curve42ePoints + [self.p25e])
        self.curve43e = gmsh.model.geo.addPolyline([self.p24e] + self.curve43ePoints + [self.p6e])
        self.curve44e = gmsh.model.geo.addPolyline([self.p25e] + self.curve44ePoints + [self.p6e])
        self.curve45e = gmsh.model.geo.addPolyline([self.p8e] + self.curve45ePoints + [self.p24e])
        self.curve46e = gmsh.model.geo.addPolyline([self.p8e] + self.curve46ePoints + [self.p25e])
        self.curve47e = gmsh.model.geo.addPolyline([self.p24e] + self.curve47ePoints + [self.p7e])
        self.curve48e = gmsh.model.geo.addPolyline([self.p25e] + self.curve48ePoints + [self.p7e])

        # Upper surface
        self.curveLoop1e = gmsh.model.geo.addCurveLoop([self.curve1e, -self.curve19e, -self.curve17e, self.curve13e])
        self.curveLoop2e = gmsh.model.geo.addCurveLoop([self.curve17e, -self.curve21e, self.curve12e, self.curve15e])
        self.curveLoop3e = gmsh.model.geo.addCurveLoop([self.curve2e, -self.curve25e, -self.curve23e, self.curve19e])
        self.curveLoop4e = gmsh.model.geo.addCurveLoop([self.curve23e, -self.curve27e, self.curve11e, self.curve21e])
        self.curveLoop5e = gmsh.model.geo.addCurveLoop([self.curve3e, -self.curve31e, -self.curve29e, self.curve25e])
        self.curveLoop6e = gmsh.model.geo.addCurveLoop([self.curve29e, -self.curve33e, self.curve10e, self.curve27e])
        self.curveLoop7e = gmsh.model.geo.addCurveLoop([self.curve4e, -self.curve37e, -self.curve35e, self.curve31e])
        self.curveLoop8e = gmsh.model.geo.addCurveLoop([self.curve35e, -self.curve39e, self.curve9e, self.curve33e])
        self.curveLoop9e = gmsh.model.geo.addCurveLoop([self.curve5e, -self.curve43e, -self.curve41e, self.curve37e])
        self.curveLoop10e = gmsh.model.geo.addCurveLoop([self.curve41e, -self.curve45e, self.curve8e, self.curve39e])
        self.curveLoop11e = gmsh.model.geo.addCurveLoop([self.curve6e, -self.curve47e, self.curve43e])
        self.curveLoop12e = gmsh.model.geo.addCurveLoop([self.curve47e, self.curve7e, self.curve45e])

        self.surface1e = gmsh.model.geo.addSurfaceFilling([self.curveLoop1e])
        self.surface2e = gmsh.model.geo.addSurfaceFilling([self.curveLoop2e])
        self.surface3e = gmsh.model.geo.addSurfaceFilling([self.curveLoop3e])
        self.surface4e = gmsh.model.geo.addSurfaceFilling([self.curveLoop4e])
        self.surface5e = gmsh.model.geo.addSurfaceFilling([self.curveLoop5e])
        self.surface6e = gmsh.model.geo.addSurfaceFilling([self.curveLoop6e])
        self.surface7e = gmsh.model.geo.addSurfaceFilling([self.curveLoop7e])
        self.surface8e = gmsh.model.geo.addSurfaceFilling([self.curveLoop8e])
        self.surface9e = gmsh.model.geo.addSurfaceFilling([self.curveLoop9e])
        self.surface10e = gmsh.model.geo.addSurfaceFilling([self.curveLoop10e])
        self.surface11e = gmsh.model.geo.addSurfaceFilling([self.curveLoop11e])
        self.surface12e = gmsh.model.geo.addSurfaceFilling([self.curveLoop12e])

        # Lower surface
        self.curveLoop13e = gmsh.model.geo.addCurveLoop([-self.curve1e, -self.curve14e, self.curve18e, self.curve20e])
        self.curveLoop14e = gmsh.model.geo.addCurveLoop([-self.curve18e, -self.curve16e, -self.curve12e, self.curve22e])
        self.curveLoop15e = gmsh.model.geo.addCurveLoop([-self.curve2e, -self.curve20e, self.curve24e, self.curve26e])
        self.curveLoop16e = gmsh.model.geo.addCurveLoop([-self.curve24e, -self.curve22e, -self.curve11e, self.curve28e])
        self.curveLoop17e = gmsh.model.geo.addCurveLoop([-self.curve3e, -self.curve26e, self.curve30e, self.curve32e])
        self.curveLoop18e = gmsh.model.geo.addCurveLoop([-self.curve30e, -self.curve28e, -self.curve10e, self.curve34e])
        self.curveLoop19e = gmsh.model.geo.addCurveLoop([-self.curve4e, -self.curve32e, self.curve36e, self.curve38e])
        self.curveLoop20e = gmsh.model.geo.addCurveLoop([-self.curve36e, -self.curve34e, -self.curve9e, self.curve40e])
        self.curveLoop21e = gmsh.model.geo.addCurveLoop([-self.curve5e, -self.curve38e, self.curve42e, self.curve44e])
        self.curveLoop22e = gmsh.model.geo.addCurveLoop([-self.curve42e, -self.curve40e, -self.curve8e, self.curve46e])
        self.curveLoop23e = gmsh.model.geo.addCurveLoop([-self.curve6e, -self.curve44e, self.curve48e])
        self.curveLoop24e = gmsh.model.geo.addCurveLoop([-self.curve48e, -self.curve46e, -self.curve7e])

        self.surface13e = gmsh.model.geo.addSurfaceFilling([self.curveLoop13e])
        self.surface14e = gmsh.model.geo.addSurfaceFilling([self.curveLoop14e])
        self.surface15e = gmsh.model.geo.addSurfaceFilling([self.curveLoop15e])
        self.surface16e = gmsh.model.geo.addSurfaceFilling([self.curveLoop16e])
        self.surface17e = gmsh.model.geo.addSurfaceFilling([self.curveLoop17e])
        self.surface18e = gmsh.model.geo.addSurfaceFilling([self.curveLoop18e])
        self.surface19e = gmsh.model.geo.addSurfaceFilling([self.curveLoop19e])
        self.surface20e = gmsh.model.geo.addSurfaceFilling([self.curveLoop20e])
        self.surface21e = gmsh.model.geo.addSurfaceFilling([self.curveLoop21e])
        self.surface22e = gmsh.model.geo.addSurfaceFilling([self.curveLoop22e])
        self.surface23e = gmsh.model.geo.addSurfaceFilling([self.curveLoop23e])
        self.surface24e = gmsh.model.geo.addSurfaceFilling([self.curveLoop24e])

        # Right wing
        self.p1d = gmsh.model.geo.addPoint(self.geo.p1d[0], self.geo.p1d[1], self.geo.p1d[2], self.lc)
        self.p2d = gmsh.model.geo.addPoint(self.geo.p2d[0], self.geo.p2d[1], self.geo.p2d[2], self.lc)
        self.p3d = gmsh.model.geo.addPoint(self.geo.p3d[0], self.geo.p3d[1], self.geo.p3d[2], self.lc)
        self.p4d = gmsh.model.geo.addPoint(self.geo.p4d[0], self.geo.p4d[1], self.geo.p4d[2], self.lc)
        self.p5d = gmsh.model.geo.addPoint(self.geo.p5d[0], self.geo.p5d[1], self.geo.p5d[2], self.lc)
        self.p6d = gmsh.model.geo.addPoint(self.geo.p6d[0], self.geo.p6d[1], self.geo.p6d[2], self.lc)
        self.p7d = gmsh.model.geo.addPoint(self.geo.p7d[0], self.geo.p7d[1], self.geo.p7d[2], self.lc)
        self.p8d = gmsh.model.geo.addPoint(self.geo.p8d[0], self.geo.p8d[1], self.geo.p8d[2], self.lc)
        self.p9d = gmsh.model.geo.addPoint(self.geo.p9d[0], self.geo.p9d[1], self.geo.p9d[2], self.lc)
        self.p10d = gmsh.model.geo.addPoint(self.geo.p10d[0], self.geo.p10d[1], self.geo.p10d[2], self.lc)
        self.p11d = gmsh.model.geo.addPoint(self.geo.p11d[0], self.geo.p11d[1], self.geo.p11d[2], self.lc)
        self.p12d = gmsh.model.geo.addPoint(self.geo.p12d[0], self.geo.p12d[1], self.geo.p12d[2], self.lc)
        self.p13d = gmsh.model.geo.addPoint(self.geo.p13d[0], self.geo.p13d[1], self.geo.p13d[2], self.lc)
        self.p14d = gmsh.model.geo.addPoint(self.geo.p14d[0], self.geo.p14d[1], self.geo.p14d[2], self.lc)
        self.p15d = gmsh.model.geo.addPoint(self.geo.p15d[0], self.geo.p15d[1], self.geo.p15d[2], self.lc)
        self.p16d = gmsh.model.geo.addPoint(self.geo.p16d[0], self.geo.p16d[1], self.geo.p16d[2], self.lc)
        self.p17d = gmsh.model.geo.addPoint(self.geo.p17d[0], self.geo.p17d[1], self.geo.p17d[2], self.lc)
        self.p18d = gmsh.model.geo.addPoint(self.geo.p18d[0], self.geo.p18d[1], self.geo.p18d[2], self.lc)
        self.p19d = gmsh.model.geo.addPoint(self.geo.p19d[0], self.geo.p19d[1], self.geo.p19d[2], self.lc)
        self.p20d = gmsh.model.geo.addPoint(self.geo.p20d[0], self.geo.p20d[1], self.geo.p20d[2], self.lc)
        self.p21d = gmsh.model.geo.addPoint(self.geo.p21d[0], self.geo.p21d[1], self.geo.p21d[2], self.lc)
        self.p22d = gmsh.model.geo.addPoint(self.geo.p22d[0], self.geo.p22d[1], self.geo.p22d[2], self.lc)
        self.p23d = gmsh.model.geo.addPoint(self.geo.p23d[0], self.geo.p23d[1], self.geo.p23d[2], self.lc)
        self.p24d = gmsh.model.geo.addPoint(self.geo.p24d[0], self.geo.p24d[1], self.geo.p24d[2], self.lc)
        self.p25d = gmsh.model.geo.addPoint(self.geo.p25d[0], self.geo.p25d[1], self.geo.p25d[2], self.lc)

        self.curve13dPoints = []
        for i in range(len(self.geo.curve13d[:, 0])):
            self.curve13dPoints.append(gmsh.model.geo.addPoint(self.geo.curve13d[i, 0], self.geo.curve13d[i, 1], self.geo.curve13d[i, 2], self.lc))

        self.curve14dPoints = []
        for i in range(len(self.geo.curve14d[:, 0])):
            self.curve14dPoints.append(gmsh.model.geo.addPoint(self.geo.curve14d[i, 0], self.geo.curve14d[i, 1], self.geo.curve14d[i, 2], self.lc))
        
        self.curve15dPoints = []
        for i in range(len(self.geo.curve15d[:, 0])):
            self.curve15dPoints.append(gmsh.model.geo.addPoint(self.geo.curve15d[i, 0], self.geo.curve15d[i, 1], self.geo.curve15d[i, 2], self.lc))
        
        self.curve16dPoints = []
        for i in range(len(self.geo.curve16d[:, 0])):
            self.curve16dPoints.append(gmsh.model.geo.addPoint(self.geo.curve16d[i, 0], self.geo.curve16d[i, 1], self.geo.curve16d[i, 2], self.lc))
        
        self.curve17dPoints = []
        for i in range(len(self.geo.curve17d[:, 0])):
            self.curve17dPoints.append(gmsh.model.geo.addPoint(self.geo.curve17d[i, 0], self.geo.curve17d[i, 1], self.geo.curve17d[i, 2], self.lc))
        
        self.curve18dPoints = []
        for i in range(len(self.geo.curve18d[:, 0])):
            self.curve18dPoints.append(gmsh.model.geo.addPoint(self.geo.curve18d[i, 0], self.geo.curve18d[i, 1], self.geo.curve18d[i, 2], self.lc))
        
        self.curve19dPoints = []
        for i in range(len(self.geo.curve19d[:, 0])):
            self.curve19dPoints.append(gmsh.model.geo.addPoint(self.geo.curve19d[i, 0], self.geo.curve19d[i, 1], self.geo.curve19d[i, 2], self.lc))
        
        self.curve20dPoints = []
        for i in range(len(self.geo.curve20d[:, 0])):
            self.curve20dPoints.append(gmsh.model.geo.addPoint(self.geo.curve20d[i, 0], self.geo.curve20d[i, 1], self.geo.curve20d[i, 2], self.lc))
        
        self.curve21dPoints = []
        for i in range(len(self.geo.curve21d[:, 0])):
            self.curve21dPoints.append(gmsh.model.geo.addPoint(self.geo.curve21d[i, 0], self.geo.curve21d[i, 1], self.geo.curve21d[i, 2], self.lc))
        
        self.curve22dPoints = []
        for i in range(len(self.geo.curve22d[:, 0])):
            self.curve22dPoints.append(gmsh.model.geo.addPoint(self.geo.curve22d[i, 0], self.geo.curve22d[i, 1], self.geo.curve22d[i, 2], self.lc))
        
        self.curve23dPoints = []
        for i in range(len(self.geo.curve23d[:, 0])):
            self.curve23dPoints.append(gmsh.model.geo.addPoint(self.geo.curve23d[i, 0], self.geo.curve23d[i, 1], self.geo.curve23d[i, 2], self.lc))
        
        self.curve24dPoints = []
        for i in range(len(self.geo.curve24d[:, 0])):
            self.curve24dPoints.append(gmsh.model.geo.addPoint(self.geo.curve24d[i, 0], self.geo.curve24d[i, 1], self.geo.curve24d[i, 2], self.lc))
        
        self.curve25dPoints = []
        for i in range(len(self.geo.curve25d[:, 0])):
            self.curve25dPoints.append(gmsh.model.geo.addPoint(self.geo.curve25d[i, 0], self.geo.curve25d[i, 1], self.geo.curve25d[i, 2], self.lc))
        
        self.curve26dPoints = []
        for i in range(len(self.geo.curve26d[:, 0])):
            self.curve26dPoints.append(gmsh.model.geo.addPoint(self.geo.curve26d[i, 0], self.geo.curve26d[i, 1], self.geo.curve26d[i, 2], self.lc))
        
        self.curve27dPoints = []
        for i in range(len(self.geo.curve27d[:, 0])):
            self.curve27dPoints.append(gmsh.model.geo.addPoint(self.geo.curve27d[i, 0], self.geo.curve27d[i, 1], self.geo.curve27d[i, 2], self.lc))
        
        self.curve28dPoints = []
        for i in range(len(self.geo.curve28d[:, 0])):
            self.curve28dPoints.append(gmsh.model.geo.addPoint(self.geo.curve28d[i, 0], self.geo.curve28d[i, 1], self.geo.curve28d[i, 2], self.lc))
        
        self.curve29dPoints = []
        for i in range(len(self.geo.curve29d[:, 0])):
            self.curve29dPoints.append(gmsh.model.geo.addPoint(self.geo.curve29d[i, 0], self.geo.curve29d[i, 1], self.geo.curve29d[i, 2], self.lc))
        
        self.curve30dPoints = []
        for i in range(len(self.geo.curve30d[:, 0])):
            self.curve30dPoints.append(gmsh.model.geo.addPoint(self.geo.curve30d[i, 0], self.geo.curve30d[i, 1], self.geo.curve30d[i, 2], self.lc))
        
        self.curve31dPoints = []
        for i in range(len(self.geo.curve31d[:, 0])):
            self.curve31dPoints.append(gmsh.model.geo.addPoint(self.geo.curve31d[i, 0], self.geo.curve31d[i, 1], self.geo.curve31d[i, 2], self.lc))
        
        self.curve32dPoints = []
        for i in range(len(self.geo.curve32d[:, 0])):
            self.curve32dPoints.append(gmsh.model.geo.addPoint(self.geo.curve32d[i, 0], self.geo.curve32d[i, 1], self.geo.curve32d[i, 2], self.lc))
        
        self.curve33dPoints = []
        for i in range(len(self.geo.curve33d[:, 0])):
            self.curve33dPoints.append(gmsh.model.geo.addPoint(self.geo.curve33d[i, 0], self.geo.curve33d[i, 1], self.geo.curve33d[i, 2], self.lc))
        
        self.curve34dPoints = []
        for i in range(len(self.geo.curve34d[:, 0])):
            self.curve34dPoints.append(gmsh.model.geo.addPoint(self.geo.curve34d[i, 0], self.geo.curve34d[i, 1], self.geo.curve34d[i, 2], self.lc))
        
        self.curve35dPoints = []
        for i in range(len(self.geo.curve35d[:, 0])):
            self.curve35dPoints.append(gmsh.model.geo.addPoint(self.geo.curve35d[i, 0], self.geo.curve35d[i, 1], self.geo.curve35d[i, 2], self.lc))
        
        self.curve36dPoints = []
        for i in range(len(self.geo.curve36d[:, 0])):
            self.curve36dPoints.append(gmsh.model.geo.addPoint(self.geo.curve36d[i, 0], self.geo.curve36d[i, 1], self.geo.curve36d[i, 2], self.lc))
        
        self.curve37dPoints = []
        for i in range(len(self.geo.curve37d[:, 0])):
            self.curve37dPoints.append(gmsh.model.geo.addPoint(self.geo.curve37d[i, 0], self.geo.curve37d[i, 1], self.geo.curve37d[i, 2], self.lc))
        
        self.curve38dPoints = []
        for i in range(len(self.geo.curve38d[:, 0])):
            self.curve38dPoints.append(gmsh.model.geo.addPoint(self.geo.curve38d[i, 0], self.geo.curve38d[i, 1], self.geo.curve38d[i, 2], self.lc))
        
        self.curve39dPoints = []
        for i in range(len(self.geo.curve39d[:, 0])):
            self.curve39dPoints.append(gmsh.model.geo.addPoint(self.geo.curve39d[i, 0], self.geo.curve39d[i, 1], self.geo.curve39d[i, 2], self.lc))
        
        self.curve40dPoints = []
        for i in range(len(self.geo.curve40d[:, 0])):
            self.curve40dPoints.append(gmsh.model.geo.addPoint(self.geo.curve40d[i, 0], self.geo.curve40d[i, 1], self.geo.curve40d[i, 2], self.lc))
        
        self.curve41dPoints = []
        for i in range(len(self.geo.curve41d[:, 0])):
            self.curve41dPoints.append(gmsh.model.geo.addPoint(self.geo.curve41d[i, 0], self.geo.curve41d[i, 1], self.geo.curve41d[i, 2], self.lc))
        
        self.curve42dPoints = []
        for i in range(len(self.geo.curve42d[:, 0])):
            self.curve42dPoints.append(gmsh.model.geo.addPoint(self.geo.curve42d[i, 0], self.geo.curve42d[i, 1], self.geo.curve42d[i, 2], self.lc))
        
        self.curve43dPoints = []
        for i in range(len(self.geo.curve43d[:, 0])):
            self.curve43dPoints.append(gmsh.model.geo.addPoint(self.geo.curve43d[i, 0], self.geo.curve43d[i, 1], self.geo.curve43d[i, 2], self.lc))
        
        self.curve44dPoints = []
        for i in range(len(self.geo.curve44d[:, 0])):
            self.curve44dPoints.append(gmsh.model.geo.addPoint(self.geo.curve44d[i, 0], self.geo.curve44d[i, 1], self.geo.curve44d[i, 2], self.lc))
        
        self.curve45dPoints = []
        for i in range(len(self.geo.curve45d[:, 0])):
            self.curve45dPoints.append(gmsh.model.geo.addPoint(self.geo.curve45d[i, 0], self.geo.curve45d[i, 1], self.geo.curve45d[i, 2], self.lc))
        
        self.curve46dPoints = []
        for i in range(len(self.geo.curve46d[:, 0])):
            self.curve46dPoints.append(gmsh.model.geo.addPoint(self.geo.curve46d[i, 0], self.geo.curve46d[i, 1], self.geo.curve46d[i, 2], self.lc))
        
        self.curve47dPoints = []
        for i in range(len(self.geo.curve47d[:, 0])):
            self.curve47dPoints.append(gmsh.model.geo.addPoint(self.geo.curve47d[i, 0], self.geo.curve47d[i, 1], self.geo.curve47d[i, 2], self.lc))
        
        self.curve48dPoints = []
        for i in range(len(self.geo.curve48d[:, 0])):
            self.curve48dPoints.append(gmsh.model.geo.addPoint(self.geo.curve48d[i, 0], self.geo.curve48d[i, 1], self.geo.curve48d[i, 2], self.lc))

        self.curve1d = gmsh.model.geo.addPolyline([self.p1d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p1d, self.geo.c1d, self.geo.p2d], linspace(0, 1, num=50))[1:49]] + [self.p2d])
        self.curve2d = gmsh.model.geo.addPolyline([self.p2d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p2d, self.geo.c2d, self.geo.c3d, self.geo.p3d], linspace(0, 1, num=50))[1:49]] + [self.p3d])
        self.curve3d = gmsh.model.geo.addPolyline([self.p3d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p3d, self.geo.c4d, self.geo.c5d, self.geo.p4d], linspace(0, 1, num=50))[1:49]] + [self.p4d])
        self.curve4d = gmsh.model.geo.addLine(self.p4d, self.p5d)
        self.curve5d = gmsh.model.geo.addPolyline([self.p5d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p5d, self.geo.c6d, self.geo.c7d, self.geo.p6d], linspace(0, 1, num=50))[1:49]] + [self.p6d])
        self.curve6d = gmsh.model.geo.addPolyline([self.p6d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p6d, self.geo.c8d, self.geo.c9d, self.geo.p7d], linspace(0, 1, num=50))[1:49]] + [self.p7d])
        self.curve7d = gmsh.model.geo.addPolyline([self.p7d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p7d, self.geo.c10d, self.geo.c11d, self.geo.p8d], linspace(0, 1, num=50))[1:49]] + [self.p8d])
        self.curve8d = gmsh.model.geo.addPolyline([self.p8d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p8d, self.geo.c12d, self.geo.c13d, self.geo.p9d], linspace(0, 1, num=50))[1:49]] + [self.p9d])
        self.curve9d = gmsh.model.geo.addPolyline([self.p9d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p9d, self.geo.c14d, self.geo.c15d, self.geo.p10d], linspace(0, 1, num=50))[1:49]] + [self.p10d])
        self.curve10d = gmsh.model.geo.addPolyline([self.p10d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p10d, self.geo.c16d, self.geo.c17d, self.geo.p11d], linspace(0, 1, num=50))[1:49]] + [self.p11d])
        self.curve11d = gmsh.model.geo.addPolyline([self.p11d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p11d, self.geo.c18d, self.geo.c19d, self.geo.p12d], linspace(0, 1, num=50))[1:49]] + [self.p12d])
        self.curve12d = gmsh.model.geo.addPolyline([self.p12d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p12d, self.geo.c20d, self.geo.c21d, self.geo.p13d], linspace(0, 1, num=50))[1:49]] + [self.p13d])

        self.curve13d = gmsh.model.geo.addPolyline([self.p14d] + self.curve13dPoints + [self.p1d])
        self.curve14d = gmsh.model.geo.addPolyline([self.p15d] + self.curve14dPoints + [self.p1d])
        self.curve15d = gmsh.model.geo.addPolyline([self.p13d] + self.curve15dPoints + [self.p14d])
        self.curve16d = gmsh.model.geo.addPolyline([self.p13d] + self.curve16dPoints + [self.p15d])
        self.curve17d = gmsh.model.geo.addPolyline([self.p14d] + self.curve17dPoints + [self.p16d])
        self.curve18d = gmsh.model.geo.addPolyline([self.p15d] + self.curve18dPoints + [self.p17d])
        self.curve19d = gmsh.model.geo.addPolyline([self.p16d] + self.curve19dPoints + [self.p2d])
        self.curve20d = gmsh.model.geo.addPolyline([self.p17d] + self.curve20dPoints + [self.p2d])
        self.curve21d = gmsh.model.geo.addPolyline([self.p12d] + self.curve21dPoints + [self.p16d])
        self.curve22d = gmsh.model.geo.addPolyline([self.p12d] + self.curve22dPoints + [self.p17d])
        self.curve23d = gmsh.model.geo.addPolyline([self.p16d] + self.curve23dPoints + [self.p18d])
        self.curve24d = gmsh.model.geo.addPolyline([self.p17d] + self.curve24dPoints + [self.p19d])
        self.curve25d = gmsh.model.geo.addPolyline([self.p18d] + self.curve25dPoints + [self.p3d])
        self.curve26d = gmsh.model.geo.addPolyline([self.p19d] + self.curve26dPoints + [self.p3d])
        self.curve27d = gmsh.model.geo.addPolyline([self.p11d] + self.curve27dPoints + [self.p18d])
        self.curve28d = gmsh.model.geo.addPolyline([self.p11d] + self.curve28dPoints + [self.p19d])
        self.curve29d = gmsh.model.geo.addPolyline([self.p18d] + self.curve29dPoints + [self.p20d])
        self.curve30d = gmsh.model.geo.addPolyline([self.p19d] + self.curve30dPoints + [self.p21d])
        self.curve31d = gmsh.model.geo.addPolyline([self.p20d] + self.curve31dPoints + [self.p4d])
        self.curve32d = gmsh.model.geo.addPolyline([self.p21d] + self.curve32dPoints + [self.p4d])
        self.curve33d = gmsh.model.geo.addPolyline([self.p10d] + self.curve33dPoints + [self.p20d])
        self.curve34d = gmsh.model.geo.addPolyline([self.p10d] + self.curve34dPoints + [self.p21d])
        self.curve35d = gmsh.model.geo.addPolyline([self.p20d] + self.curve35dPoints + [self.p22d])
        self.curve36d = gmsh.model.geo.addPolyline([self.p21d] + self.curve36dPoints + [self.p23d])
        self.curve37d = gmsh.model.geo.addPolyline([self.p22d] + self.curve37dPoints + [self.p5d])
        self.curve38d = gmsh.model.geo.addPolyline([self.p23d] + self.curve38dPoints + [self.p5d])
        self.curve39d = gmsh.model.geo.addPolyline([self.p9d] + self.curve39dPoints + [self.p22d])
        self.curve40d = gmsh.model.geo.addPolyline([self.p9d] + self.curve40dPoints + [self.p23d])
        self.curve41d = gmsh.model.geo.addPolyline([self.p22d] + self.curve41dPoints + [self.p24d])
        self.curve42d = gmsh.model.geo.addPolyline([self.p23d] + self.curve42dPoints + [self.p25d])
        self.curve43d = gmsh.model.geo.addPolyline([self.p24d] + self.curve43dPoints + [self.p6d])
        self.curve44d = gmsh.model.geo.addPolyline([self.p25d] + self.curve44dPoints + [self.p6d])
        self.curve45d = gmsh.model.geo.addPolyline([self.p8d] + self.curve45dPoints + [self.p24d])
        self.curve46d = gmsh.model.geo.addPolyline([self.p8d] + self.curve46dPoints + [self.p25d])
        self.curve47d = gmsh.model.geo.addPolyline([self.p24d] + self.curve47dPoints + [self.p7d])
        self.curve48d = gmsh.model.geo.addPolyline([self.p25d] + self.curve48dPoints + [self.p7d])

        self.curveLoop1d = gmsh.model.geo.addCurveLoop([-self.curve1d, -self.curve13d, self.curve17d, self.curve19d])
        self.curveLoop2d = gmsh.model.geo.addCurveLoop([-self.curve17d, -self.curve15d, -self.curve12d, self.curve21d])
        self.curveLoop3d = gmsh.model.geo.addCurveLoop([-self.curve2d, -self.curve19d, self.curve23d, self.curve25d])
        self.curveLoop4d = gmsh.model.geo.addCurveLoop([-self.curve23d, -self.curve21d, -self.curve11d, self.curve27d])
        self.curveLoop5d = gmsh.model.geo.addCurveLoop([-self.curve3d, -self.curve25d, self.curve29d, self.curve31d])
        self.curveLoop6d = gmsh.model.geo.addCurveLoop([-self.curve29d, -self.curve27d, -self.curve10d, self.curve33d])
        self.curveLoop7d = gmsh.model.geo.addCurveLoop([-self.curve4d, -self.curve31d, self.curve35d, self.curve37d])
        self.curveLoop8d = gmsh.model.geo.addCurveLoop([-self.curve35d, -self.curve33d, -self.curve9d, self.curve39d])
        self.curveLoop9d = gmsh.model.geo.addCurveLoop([-self.curve5d, -self.curve37d, self.curve41d, self.curve43d])
        self.curveLoop10d = gmsh.model.geo.addCurveLoop([-self.curve41d, -self.curve39d, -self.curve8d, self.curve45d])
        self.curveLoop11d = gmsh.model.geo.addCurveLoop([-self.curve6d, -self.curve43d, self.curve47d])
        self.curveLoop12d = gmsh.model.geo.addCurveLoop([-self.curve47d, -self.curve45d, -self.curve7d])

        self.surface1d = gmsh.model.geo.addSurfaceFilling([self.curveLoop1d])
        self.surface2d = gmsh.model.geo.addSurfaceFilling([self.curveLoop2d])
        self.surface3d = gmsh.model.geo.addSurfaceFilling([self.curveLoop3d])
        self.surface4d = gmsh.model.geo.addSurfaceFilling([self.curveLoop4d])
        self.surface5d = gmsh.model.geo.addSurfaceFilling([self.curveLoop5d])
        self.surface6d = gmsh.model.geo.addSurfaceFilling([self.curveLoop6d])
        self.surface7d = gmsh.model.geo.addSurfaceFilling([self.curveLoop7d])
        self.surface8d = gmsh.model.geo.addSurfaceFilling([self.curveLoop8d])
        self.surface9d = gmsh.model.geo.addSurfaceFilling([self.curveLoop9d])
        self.surface10d = gmsh.model.geo.addSurfaceFilling([self.curveLoop10d])
        self.surface11d = gmsh.model.geo.addSurfaceFilling([self.curveLoop11d])
        self.surface12d = gmsh.model.geo.addSurfaceFilling([self.curveLoop12d])

        self.curveLoop13d = gmsh.model.geo.addCurveLoop([self.curve1d, -self.curve20d, -self.curve18d, self.curve14d])
        self.curveLoop14d = gmsh.model.geo.addCurveLoop([self.curve18d, -self.curve22d, self.curve12d, self.curve16d])
        self.curveLoop15d = gmsh.model.geo.addCurveLoop([self.curve2d, -self.curve26d, -self.curve24d, self.curve20d])
        self.curveLoop16d = gmsh.model.geo.addCurveLoop([self.curve24d, -self.curve28d, self.curve11d, self.curve22d])
        self.curveLoop17d = gmsh.model.geo.addCurveLoop([self.curve3d, -self.curve32d, -self.curve30d, self.curve26d])
        self.curveLoop18d = gmsh.model.geo.addCurveLoop([self.curve30d, -self.curve34d, self.curve10d, self.curve28d])
        self.curveLoop19d = gmsh.model.geo.addCurveLoop([self.curve4d, -self.curve38d, -self.curve36d, self.curve32d])
        self.curveLoop20d = gmsh.model.geo.addCurveLoop([self.curve36d, -self.curve40d, self.curve9d, self.curve34d])
        self.curveLoop21d = gmsh.model.geo.addCurveLoop([self.curve5d, -self.curve44d, -self.curve42d, self.curve38d])
        self.curveLoop22d = gmsh.model.geo.addCurveLoop([self.curve42d, -self.curve46d, self.curve8d, self.curve40d])
        self.curveLoop23d = gmsh.model.geo.addCurveLoop([self.curve6d, -self.curve48d, self.curve44d])
        self.curveLoop24d = gmsh.model.geo.addCurveLoop([self.curve48d, self.curve7d, self.curve46d])

        self.surface13d = gmsh.model.geo.addSurfaceFilling([self.curveLoop13d])
        self.surface14d = gmsh.model.geo.addSurfaceFilling([self.curveLoop14d])
        self.surface15d = gmsh.model.geo.addSurfaceFilling([self.curveLoop15d])
        self.surface16d = gmsh.model.geo.addSurfaceFilling([self.curveLoop16d])
        self.surface17d = gmsh.model.geo.addSurfaceFilling([self.curveLoop17d])
        self.surface18d = gmsh.model.geo.addSurfaceFilling([self.curveLoop18d])
        self.surface19d = gmsh.model.geo.addSurfaceFilling([self.curveLoop19d])
        self.surface20d = gmsh.model.geo.addSurfaceFilling([self.curveLoop20d])
        self.surface21d = gmsh.model.geo.addSurfaceFilling([self.curveLoop21d])
        self.surface22d = gmsh.model.geo.addSurfaceFilling([self.curveLoop22d])
        self.surface23d = gmsh.model.geo.addSurfaceFilling([self.curveLoop23d])
        self.surface24d = gmsh.model.geo.addSurfaceFilling([self.curveLoop24d])

        # Create physical group
        # self.surface1e, self.surface2e,
        #         self.surface3e, self.surface4e,
        #         self.surface13e, self.surface14e,
        #         self.surface15e, self.surface16e,

        # self.physicalGroup2 = gmsh.model.addPhysicalGroup( # Wing middle section
        #     2,
        #     [
        #         self.surface5e, self.surface6e,
        #         self.surface17e, self.surface18e,
        #     ]
        # )
        # self.physicalGroup3 = gmsh.model.addPhysicalGroup( # Wing tip
        #     2,
        #     [
        #         self.surface7e, self.surface8e,
        #         self.surface9e, self.surface10e,
        #         self.surface11e, self.surface12e,
        #         self.surface19e, self.surface20e,
        #         self.surface21e, self.surface22e,
        #         self.surface23e, self.surface24e,
        #     ]
        # )

        return
    
    def _body(self) -> None:

        # Surface points
        self.p26e = gmsh.model.geo.addPoint(self.geo.p26e[0], self.geo.p26e[1], self.geo.p26e[2], self.lc)
        self.p27e = gmsh.model.geo.addPoint(self.geo.p27e[0], self.geo.p27e[1], self.geo.p27e[2], self.lc)
        self.p26d = gmsh.model.geo.addPoint(self.geo.p26d[0], self.geo.p26d[1], self.geo.p26d[2], self.lc)
        self.p27d = gmsh.model.geo.addPoint(self.geo.p27d[0], self.geo.p27d[1], self.geo.p27d[2], self.lc)
        self.p28 = gmsh.model.geo.addPoint(self.geo.p28[0], self.geo.p28[1], self.geo.p28[2], self.lc)
        self.p29 = gmsh.model.geo.addPoint(self.geo.p29[0], self.geo.p29[1], self.geo.p29[2], self.lc)
        self.p30 = gmsh.model.geo.addPoint(self.geo.p30[0], self.geo.p30[1], self.geo.p30[2], self.lc)
        self.p31 = gmsh.model.geo.addPoint(self.geo.p31[0], self.geo.p31[1], self.geo.p31[2], self.lc)
        self.p32 = gmsh.model.geo.addPoint(self.geo.p32[0], self.geo.p32[1], self.geo.p32[2], self.lc)
        self.p33 = gmsh.model.geo.addPoint(self.geo.p33[0], self.geo.p33[1], self.geo.p33[2], self.lc)
        self.p34 = gmsh.model.geo.addPoint(self.geo.p34[0], self.geo.p34[1], self.geo.p34[2], self.lc)
        self.p35 = gmsh.model.geo.addPoint(self.geo.p35[0], self.geo.p35[1], self.geo.p35[2], self.lc)
        self.p36 = gmsh.model.geo.addPoint(self.geo.p36[0], self.geo.p36[1], self.geo.p36[2], self.lc)
        self.p37 = gmsh.model.geo.addPoint(self.geo.p37[0], self.geo.p37[1], self.geo.p37[2], self.lc)
        
        # Curves
        self.curve49e = gmsh.model.geo.addPolyline([self.p26e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p26e, self.geo.c22e, self.geo.c23e, self.geo.p1e], linspace(0, 1, num=30))[1:29]] + [self.p1e])
        self.curve50e = gmsh.model.geo.addPolyline([self.p13e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p13e, self.geo.c24e, self.geo.c25e, self.geo.p27e], linspace(0, 1, num=30))[1:29]] + [self.p27e])
        self.curve49d = gmsh.model.geo.addPolyline([self.p26d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p26d, self.geo.c22d, self.geo.c23d, self.geo.p1d], linspace(0, 1, num=30))[1:29]] + [self.p1d])
        self.curve50d = gmsh.model.geo.addPolyline([self.p13d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p13d, self.geo.c24d, self.geo.c25d, self.geo.p27d], linspace(0, 1, num=30))[1:29]] + [self.p27d])
        self.curve51 = gmsh.model.geo.addPolyline([self.p28] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p28, self.geo.c26, self.geo.c27, self.geo.p29], linspace(0, 1, num=30))[1:29]] + [self.p29])
        self.curve52 = gmsh.model.geo.addPolyline([self.p29] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p29, self.geo.c28, self.geo.c29, self.geo.p30], linspace(0, 1, num=30))[1:29]] + [self.p30])
        self.curve53 = gmsh.model.geo.addPolyline([self.p30] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p30, self.geo.c30, self.geo.c31, self.geo.p31], linspace(0, 1, num=30))[1:29]] + [self.p31])
        self.curve54 = gmsh.model.geo.addPolyline([self.p31] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p31, self.geo.c32, self.geo.c33, self.geo.p32], linspace(0, 1, num=30))[1:29]] + [self.p32])
        self.curve55 = gmsh.model.geo.addPolyline([self.p33] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p33, self.geo.c34, self.geo.c35, self.geo.p34], linspace(0, 1, num=30))[1:29]] + [self.p34])
        self.curve56 = gmsh.model.geo.addPolyline([self.p34] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p34, self.geo.c36, self.geo.c37, self.geo.p35], linspace(0, 1, num=30))[1:29]] + [self.p35])
        self.curve57 = gmsh.model.geo.addPolyline([self.p35] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p35, self.geo.c38, self.geo.c39, self.geo.p36], linspace(0, 1, num=30))[1:29]] + [self.p36])
        self.curve58 = gmsh.model.geo.addPolyline([self.p36] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p36, self.geo.c40, self.geo.c41, self.geo.p37], linspace(0, 1, num=30))[1:29]] + [self.p37])

        centerPoint = 0.5 * (self.geo.p26e + self.geo.p26d)
        centerTag = gmsh.model.geo.addPoint(centerPoint[0], centerPoint[1], centerPoint[2], self.lc)
        
        self.curve59e = gmsh.model.geo.addCircleArc(self.p28, centerTag, self.p26e)
        self.curve60e = gmsh.model.geo.addCircleArc(self.p26e, centerTag, self.p33)
        self.curve59d = gmsh.model.geo.addCircleArc(self.p28, centerTag, self.p26d)
        self.curve60d = gmsh.model.geo.addCircleArc(self.p26d, centerTag, self.p33)

        self.curve61e = gmsh.model.geo.addPolyline([self.p29] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p29, self.geo.c42e, self.geo.c43e, self.geo.p1e], linspace(0, 1, num=30))[1:29]] + [self.p1e])
        self.curve62e = gmsh.model.geo.addPolyline([self.p1e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p1e, self.geo.c44e, self.geo.c45e, self.geo.p34], linspace(0, 1, num=30))[1:29]] + [self.p34])
        self.curve61d = gmsh.model.geo.addPolyline([self.p29] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p29, self.geo.c42d, self.geo.c43d, self.geo.p1d], linspace(0, 1, num=30))[1:29]] + [self.p1d])
        self.curve62d = gmsh.model.geo.addPolyline([self.p1d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p1d, self.geo.c44d, self.geo.c45d, self.geo.p34], linspace(0, 1, num=30))[1:29]] + [self.p34])
        
        self.curve63e = gmsh.model.geo.addPolyline([self.p30] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p30, self.geo.c46e, self.geo.c47e, self.geo.p14e], linspace(0, 1, num=30))[1:29]] + [self.p14e])
        self.curve64e = gmsh.model.geo.addPolyline([self.p15e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p15e, self.geo.c48e, self.geo.c49e, self.geo.p35], linspace(0, 1, num=30))[1:29]] + [self.p35])
        self.curve63d = gmsh.model.geo.addPolyline([self.p30] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p30, self.geo.c46d, self.geo.c47d, self.geo.p14d], linspace(0, 1, num=30))[1:29]] + [self.p14d])
        self.curve64d = gmsh.model.geo.addPolyline([self.p15d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p15d, self.geo.c48d, self.geo.c49d, self.geo.p35], linspace(0, 1, num=30))[1:29]] + [self.p35])

        self.curve65e = gmsh.model.geo.addPolyline([self.p31] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p31, self.geo.c50e, self.geo.c51e, self.geo.p13e], linspace(0, 1, num=30))[1:29]] + [self.p13e])
        self.curve66e = gmsh.model.geo.addPolyline([self.p13e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p13e, self.geo.c52e, self.geo.c53e, self.geo.p36], linspace(0, 1, num=30))[1:29]] + [self.p36])
        self.curve65d = gmsh.model.geo.addPolyline([self.p31] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p31, self.geo.c50d, self.geo.c51d, self.geo.p13d], linspace(0, 1, num=30))[1:29]] + [self.p13d])
        self.curve66d = gmsh.model.geo.addPolyline([self.p13d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.cubic([self.geo.p13d, self.geo.c52d, self.geo.c53d, self.geo.p36], linspace(0, 1, num=30))[1:29]] + [self.p36])

        centerPoint = 0.5 * (self.geo.p27e + self.geo.p27d)
        centerTag = gmsh.model.geo.addPoint(centerPoint[0], centerPoint[1], centerPoint[2], self.lc)

        self.curve67e = gmsh.model.geo.addEllipseArc(self.p32, centerTag, self.p27e, self.p27e)
        self.curve68e = gmsh.model.geo.addEllipseArc(self.p27e, centerTag, self.p27e, self.p37)
        self.curve67d = gmsh.model.geo.addEllipseArc(self.p32, centerTag, self.p27d, self.p27d)
        self.curve68d = gmsh.model.geo.addEllipseArc(self.p27d, centerTag, self.p27d, self.p37)
        
        # Curve loops
        self.curveLoop25e = gmsh.model.geo.addCurveLoop([self.curve59e, self.curve49e, -self.curve61e, -self.curve51])
        self.curveLoop26e = gmsh.model.geo.addCurveLoop([-self.curve13e, -self.curve63e, -self.curve52, self.curve61e])
        self.curveLoop27e = gmsh.model.geo.addCurveLoop([-self.curve15e, -self.curve65e, -self.curve53, self.curve63e])
        self.curveLoop28e = gmsh.model.geo.addCurveLoop([self.curve50e, -self.curve67e, -self.curve54, self.curve65e])
        self.curveLoop25d = gmsh.model.geo.addCurveLoop([-self.curve59d, self.curve51, self.curve61d, -self.curve49d])
        self.curveLoop26d = gmsh.model.geo.addCurveLoop([self.curve13d, -self.curve61d, self.curve52, self.curve63d])
        self.curveLoop27d = gmsh.model.geo.addCurveLoop([self.curve15d, -self.curve63d, self.curve53, self.curve65d])
        self.curveLoop28d = gmsh.model.geo.addCurveLoop([-self.curve50d, -self.curve65d, self.curve54, self.curve67d])

        self.curveLoop29e = gmsh.model.geo.addCurveLoop([self.curve60e, self.curve55, -self.curve62e, -self.curve49e])
        self.curveLoop30e = gmsh.model.geo.addCurveLoop([self.curve62e, self.curve56, -self.curve64e, self.curve14e])
        self.curveLoop31e = gmsh.model.geo.addCurveLoop([self.curve64e, self.curve57, -self.curve66e, self.curve16e])
        self.curveLoop32e = gmsh.model.geo.addCurveLoop([self.curve66e, self.curve58, -self.curve68e, -self.curve50e])
        self.curveLoop29d = gmsh.model.geo.addCurveLoop([-self.curve60d, self.curve49d, self.curve62d, -self.curve55])
        self.curveLoop30d = gmsh.model.geo.addCurveLoop([-self.curve62d, -self.curve14d, self.curve64d, -self.curve56])
        self.curveLoop31d = gmsh.model.geo.addCurveLoop([-self.curve64d, -self.curve16d, self.curve66d, -self.curve57])
        self.curveLoop32d = gmsh.model.geo.addCurveLoop([-self.curve66d, self.curve50d, self.curve68d, -self.curve58])

        # Surfaces
        self.surface25e = gmsh.model.geo.addSurfaceFilling([self.curveLoop25e])
        self.surface26e = gmsh.model.geo.addSurfaceFilling([self.curveLoop26e])
        self.surface27e = gmsh.model.geo.addSurfaceFilling([self.curveLoop27e])
        self.surface28e = gmsh.model.geo.addSurfaceFilling([self.curveLoop28e])
        self.surface29e = gmsh.model.geo.addSurfaceFilling([self.curveLoop29e])
        self.surface30e = gmsh.model.geo.addSurfaceFilling([self.curveLoop30e])
        self.surface31e = gmsh.model.geo.addSurfaceFilling([self.curveLoop31e])
        self.surface32e = gmsh.model.geo.addSurfaceFilling([self.curveLoop32e])
        self.surface25d = gmsh.model.geo.addSurfaceFilling([self.curveLoop25d])
        self.surface26d = gmsh.model.geo.addSurfaceFilling([self.curveLoop26d])
        self.surface27d = gmsh.model.geo.addSurfaceFilling([self.curveLoop27d])
        self.surface28d = gmsh.model.geo.addSurfaceFilling([self.curveLoop28d])
        self.surface29d = gmsh.model.geo.addSurfaceFilling([self.curveLoop29d])
        self.surface30d = gmsh.model.geo.addSurfaceFilling([self.curveLoop30d])
        self.surface31d = gmsh.model.geo.addSurfaceFilling([self.curveLoop31d])
        self.surface32d = gmsh.model.geo.addSurfaceFilling([self.curveLoop32d])

        return

    def _head(self) -> None:

        self.p38 = gmsh.model.geo.addPoint(self.geo.p38[0], self.geo.p38[1], self.geo.p38[2], self.lc)
        self.p39 = gmsh.model.geo.addPoint(self.geo.p39[0], self.geo.p39[1], self.geo.p39[2], self.lc)
        self.p41e = gmsh.model.geo.addPoint(self.geo.p41e[0], self.geo.p41e[1], self.geo.p41e[2], self.lc)
        self.p41d = gmsh.model.geo.addPoint(self.geo.p41d[0], self.geo.p41d[1], self.geo.p41d[2], self.lc)
        self.p42 = gmsh.model.geo.addPoint(self.geo.p42[0], self.geo.p42[1], self.geo.p42[2], self.lc)

        self.curve69 = gmsh.model.geo.addPolyline([self.p28] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p28, self.geo.c54, self.geo.p38], linspace(0, 1, num=30))[1:29]] + [self.p38])
        self.curve70 = gmsh.model.geo.addPolyline([self.p33] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p33, self.geo.c55, self.geo.p39], linspace(0, 1, num=30))[1:29]] + [self.p39])
        self.curve71 = gmsh.model.geo.addPolyline([self.p26e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p26e, self.geo.c56e, self.geo.p41e], linspace(0, 1, num=30))[1:29]] + [self.p41e])
        self.curve72 = gmsh.model.geo.addPolyline([self.p26d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p26d, self.geo.c56d, self.geo.p41d], linspace(0, 1, num=30))[1:29]] + [self.p41d])
        
        center = 0.5 * (self.geo.p38 + self.geo.p39)
        centerTag = gmsh.model.geo.addPoint(center[0], center[1], center[2], self.lc)

        self.curve73 = gmsh.model.geo.addCircleArc(self.p38, centerTag, self.p41e)
        self.curve74 = gmsh.model.geo.addCircleArc(self.p41e, centerTag, self.p39)
        self.curve75 = gmsh.model.geo.addCircleArc(self.p38, centerTag, self.p41d)
        self.curve76 = gmsh.model.geo.addCircleArc(self.p41d, centerTag, self.p39)
        self.curve77 = gmsh.model.geo.addPolyline([self.p38] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p38, self.geo.c57, self.geo.p42], linspace(0, 1, num=30))[1:29]] + [self.p42])
        self.curve78 = gmsh.model.geo.addPolyline([self.p39] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p39, self.geo.c58, self.geo.p42], linspace(0, 1, num=30))[1:29]] + [self.p42])
        self.curve79 = gmsh.model.geo.addPolyline([self.p41e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p41e, self.geo.c59e, self.geo.p42], linspace(0, 1, num=30))[1:29]] + [self.p42])
        self.curve80 = gmsh.model.geo.addPolyline([self.p41d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in bezier.quadratic([self.geo.p41d, self.geo.c59d, self.geo.p42], linspace(0, 1, num=30))[1:29]] + [self.p42])

        self.curveLoop33 = gmsh.model.geo.addCurveLoop([-self.curve59e, self.curve69, self.curve73, -self.curve71])
        self.curveLoop34 = gmsh.model.geo.addCurveLoop([-self.curve60e, self.curve71, self.curve74, -self.curve70])
        self.curveLoop35 = gmsh.model.geo.addCurveLoop([self.curve60d, self.curve70, -self.curve76, -self.curve72])
        self.curveLoop36 = gmsh.model.geo.addCurveLoop([self.curve59d, self.curve72, -self.curve75, -self.curve69])

        self.surface33 = gmsh.model.geo.addSurfaceFilling([self.curveLoop33])
        self.surface34 = gmsh.model.geo.addSurfaceFilling([self.curveLoop34])
        self.surface35 = gmsh.model.geo.addSurfaceFilling([self.curveLoop35])
        self.surface36 = gmsh.model.geo.addSurfaceFilling([self.curveLoop36])

        self.curveLoop37 = gmsh.model.geo.addCurveLoop([-self.curve73, self.curve77, -self.curve79])
        self.curveLoop38 = gmsh.model.geo.addCurveLoop([-self.curve74, self.curve79, -self.curve78])
        self.curveLoop39 = gmsh.model.geo.addCurveLoop([self.curve76, self.curve78, -self.curve80])
        self.curveLoop40 = gmsh.model.geo.addCurveLoop([self.curve75, self.curve80, -self.curve77])

        self.surface37 = gmsh.model.geo.addSurfaceFilling([self.curveLoop37])
        self.surface38 = gmsh.model.geo.addSurfaceFilling([self.curveLoop38])
        self.surface39 = gmsh.model.geo.addSurfaceFilling([self.curveLoop39])
        self.surface40 = gmsh.model.geo.addSurfaceFilling([self.curveLoop40])

        return
    
    def _tail(self) -> None:

        self.p47e = gmsh.model.geo.addPoint(self.geo.p47e[0], self.geo.p47e[1], self.geo.p47e[2], self.lc)
        self.p49e = gmsh.model.geo.addPoint(self.geo.p49e[0], self.geo.p49e[1], self.geo.p49e[2], self.lc)
        self.p51 = gmsh.model.geo.addPoint(self.geo.p51[0], self.geo.p51[1], self.geo.p51[2], self.lc)
        self.p49d = gmsh.model.geo.addPoint(self.geo.p49d[0], self.geo.p49d[1], self.geo.p49d[2], self.lc)
        self.p47d = gmsh.model.geo.addPoint(self.geo.p47d[0], self.geo.p47d[1], self.geo.p47d[2], self.lc)
        self.p43e = gmsh.model.geo.addPoint(self.geo.p43e[0], self.geo.p43e[1], self.geo.p43e[2], self.lc)
        self.p44e = gmsh.model.geo.addPoint(self.geo.p44e[0], self.geo.p44e[1], self.geo.p44e[2], self.lc)
        self.p43d = gmsh.model.geo.addPoint(self.geo.p43d[0], self.geo.p43d[1], self.geo.p43d[2], self.lc)
        self.p44d = gmsh.model.geo.addPoint(self.geo.p44d[0], self.geo.p44d[1], self.geo.p44d[2], self.lc)
        self.p45 = gmsh.model.geo.addPoint(self.geo.p45[0], self.geo.p45[1], self.geo.p45[2], self.lc)
        self.p46 = gmsh.model.geo.addPoint(self.geo.p46[0], self.geo.p46[1], self.geo.p46[2], self.lc)

        self.curve81e = gmsh.model.geo.addLine(self.p27e, self.p47e)
        self.curve81d = gmsh.model.geo.addLine(self.p27d, self.p47d)
        self.curve82e = gmsh.model.geo.addPolyline([self.p47e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve82e] + [self.p49e])
        self.curve82d = gmsh.model.geo.addPolyline([self.p47d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve82d] + [self.p49d])
        self.curve83e = gmsh.model.geo.addPolyline([self.p49e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve83e] + [self.p51])
        self.curve83d = gmsh.model.geo.addPolyline([self.p49d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve83d] + [self.p51])
        
        self.curve84e = gmsh.model.geo.addPolyline([self.p43e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve84e] + [self.p27e])
        self.curve85e = gmsh.model.geo.addPolyline([self.p44e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve85e] + [self.p27e])
        self.curve86e = gmsh.model.geo.addPolyline([self.p49e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve86e] + [self.p43e])
        self.curve87e = gmsh.model.geo.addPolyline([self.p49e] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve87e] + [self.p44e])
        self.curve84d = gmsh.model.geo.addPolyline([self.p43d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve84d] + [self.p27d])
        self.curve85d = gmsh.model.geo.addPolyline([self.p44d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve85d] + [self.p27d])
        self.curve86d = gmsh.model.geo.addPolyline([self.p49d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve86d] + [self.p43d])
        self.curve87d = gmsh.model.geo.addPolyline([self.p49d] + [gmsh.model.geo.addPoint(i[0], i[1], i[2], self.lc) for i in self.geo.curve87d] + [self.p44d])
        self.curve88e = gmsh.model.geo.addLine(self.p43e, self.p47e)
        self.curve89e = gmsh.model.geo.addLine(self.p44e, self.p47e)
        self.curve88d = gmsh.model.geo.addLine(self.p43d, self.p47d)
        self.curve89d = gmsh.model.geo.addLine(self.p44d, self.p47d)

        self.curve90e = gmsh.model.geo.addLine(self.p45, self.p43e)
        self.curve91e = gmsh.model.geo.addLine(self.p46, self.p44e)
        self.curve90d = gmsh.model.geo.addLine(self.p45, self.p43d)
        self.curve91d = gmsh.model.geo.addLine(self.p46, self.p44d)

        self.curve92 = gmsh.model.geo.addLine(self.p32, self.p45)
        self.curve93 = gmsh.model.geo.addLine(self.p37, self.p46)
        self.curve94 = gmsh.model.geo.addLine(self.p45, self.p51)
        self.curve95 = gmsh.model.geo.addLine(self.p46, self.p51)

        self.curveLoop41e = gmsh.model.geo.addCurveLoop([self.curve81e, -self.curve88e, self.curve84e])
        self.curveLoop42e = gmsh.model.geo.addCurveLoop([self.curve82e, self.curve86e, self.curve88e])
        self.curveLoop43e = gmsh.model.geo.addCurveLoop([self.curve67e, -self.curve84e, -self.curve90e, -self.curve92])
        self.curveLoop44e = gmsh.model.geo.addCurveLoop([self.curve90e, -self.curve86e, self.curve83e, -self.curve94])

        self.curveLoop41d = gmsh.model.geo.addCurveLoop([-self.curve81d, -self.curve84d, self.curve88d])
        self.curveLoop42d = gmsh.model.geo.addCurveLoop([-self.curve82d, -self.curve88d, -self.curve86d])
        self.curveLoop43d = gmsh.model.geo.addCurveLoop([-self.curve67d, self.curve92, self.curve90d, self.curve84d])
        self.curveLoop44d = gmsh.model.geo.addCurveLoop([-self.curve90d, self.curve94, -self.curve83d, self.curve86d])

        self.surface41e = gmsh.model.geo.addSurfaceFilling([self.curveLoop41e])
        self.surface42e = gmsh.model.geo.addSurfaceFilling([self.curveLoop42e])
        self.surface43e = gmsh.model.geo.addSurfaceFilling([self.curveLoop43e])
        self.surface44e = gmsh.model.geo.addSurfaceFilling([self.curveLoop44e])

        self.surface41d = gmsh.model.geo.addSurfaceFilling([self.curveLoop41d])
        self.surface42d = gmsh.model.geo.addSurfaceFilling([self.curveLoop42d])
        self.surface43d = gmsh.model.geo.addSurfaceFilling([self.curveLoop43d])
        self.surface44d = gmsh.model.geo.addSurfaceFilling([self.curveLoop44d])

        self.curveLoop45e = gmsh.model.geo.addCurveLoop([self.curve89e, -self.curve81e, -self.curve85e])
        self.curveLoop46e = gmsh.model.geo.addCurveLoop([-self.curve82e, -self.curve89e, -self.curve87e])
        self.curveLoop47e = gmsh.model.geo.addCurveLoop([self.curve95, -self.curve83e, self.curve87e, -self.curve91e])
        self.curveLoop48e = gmsh.model.geo.addCurveLoop([self.curve68e, self.curve93, self.curve91e, self.curve85e])

        self.curveLoop45d = gmsh.model.geo.addCurveLoop([-self.curve89d, self.curve85d, self.curve81d])
        self.curveLoop46d = gmsh.model.geo.addCurveLoop([self.curve82d, self.curve87d, self.curve89d])
        self.curveLoop47d = gmsh.model.geo.addCurveLoop([-self.curve95, self.curve91d, -self.curve87d, self.curve83d])
        self.curveLoop48d = gmsh.model.geo.addCurveLoop([-self.curve68d, -self.curve85d, -self.curve91d, -self.curve93])

        self.surface45e = gmsh.model.geo.addSurfaceFilling([self.curveLoop45e])
        self.surface46e = gmsh.model.geo.addSurfaceFilling([self.curveLoop46e])
        self.surface47e = gmsh.model.geo.addSurfaceFilling([self.curveLoop47e])
        self.surface48e = gmsh.model.geo.addSurfaceFilling([self.curveLoop48e])

        self.surface45d = gmsh.model.geo.addSurfaceFilling([self.curveLoop45d])
        self.surface46d = gmsh.model.geo.addSurfaceFilling([self.curveLoop46d])
        self.surface47d = gmsh.model.geo.addSurfaceFilling([self.curveLoop47d])
        self.surface48d = gmsh.model.geo.addSurfaceFilling([self.curveLoop48d])

        return

    def _refine(self) -> None:

        chord = self.geo.data.wing.h1 + self.geo.data.wing.h7
        wingMinDist = chord * 0.05
        wingMaxDist = chord * 0.30
        tailMinDist = self.geo.data.tail.h20 * 0.1
        tailMaxDist = self.geo.data.tail.h20 * 0.40

        # Wing leading edge
        gmsh.model.mesh.field.add("Distance", 1)
        gmsh.model.mesh.field.setNumbers(1, "CurvesList", [self.curve1e, self.curve2e, self.curve3e, self.curve4e, self.curve5e, self.curve6e, self.curve1d, self.curve2d, self.curve3d, self.curve4d, self.curve5d, self.curve6d, self.curve49e, self.curve49d])
        gmsh.model.mesh.field.setNumber(1, "NumPointsPerCurve", 200)

        gmsh.model.mesh.field.add("Threshold", 2)
        gmsh.model.mesh.field.setNumber(2, "InField", 1)
        gmsh.model.mesh.field.setNumber(2, "SizeMin", self.wingLEcellSize)
        gmsh.model.mesh.field.setNumber(2, "SizeMax", self.cellSize)
        gmsh.model.mesh.field.setNumber(2, "DistMin", wingMinDist)
        gmsh.model.mesh.field.setNumber(2, "DistMax", wingMaxDist)

        # Wing trailing edge
        gmsh.model.mesh.field.add("Distance", 3)
        gmsh.model.mesh.field.setNumbers(3, "CurvesList", [self.curve7e, self.curve8e, self.curve9e, self.curve10e, self.curve11e, self.curve12e, self.curve7d, self.curve8d, self.curve9d, self.curve10d, self.curve11d, self.curve12d, self.curve50e, self.curve50d])
        gmsh.model.mesh.field.setNumber(3, "NumPointsPerCurve", 200)

        gmsh.model.mesh.field.add("Threshold", 4)
        gmsh.model.mesh.field.setNumber(4, "InField", 3)
        gmsh.model.mesh.field.setNumber(4, "SizeMin", self.wingTEcellSize)
        gmsh.model.mesh.field.setNumber(4, "SizeMax", self.cellSize)
        gmsh.model.mesh.field.setNumber(4, "DistMin", wingMinDist)
        gmsh.model.mesh.field.setNumber(4, "DistMax", wingMaxDist)

        # Head
        headMinDist = self.geo.data.head.h18
        headMaxDist = 2 * headMinDist

        gmsh.model.mesh.field.add("Distance", 5)
        gmsh.model.mesh.field.setNumbers(5, "CurvesList", [self.curve69, self.curve71, self.curve72, self.curve70, self.curve77, self.curve78, self.curve79, self.curve80])
        gmsh.model.mesh.field.setNumber(5, "NumPointsPerCurve", 200)

        gmsh.model.mesh.field.add("Threshold", 6)
        gmsh.model.mesh.field.setNumber(6, "InField", 5)
        gmsh.model.mesh.field.setNumber(6, "SizeMin", self.headCellSize)
        gmsh.model.mesh.field.setNumber(6, "SizeMax", self.cellSize)
        gmsh.model.mesh.field.setNumber(6, "DistMin", headMinDist)
        gmsh.model.mesh.field.setNumber(6, "DistMax", headMaxDist)

        # Tail leading edge
        gmsh.model.mesh.field.add("Distance", 7)
        gmsh.model.mesh.field.setNumbers(7, "CurvesList", [self.curve81e, self.curve81d])
        gmsh.model.mesh.field.setNumber(7, "NumPointsPerCurve", 200)

        gmsh.model.mesh.field.add("Threshold", 8)
        gmsh.model.mesh.field.setNumber(8, "InField", 7)
        gmsh.model.mesh.field.setNumber(8, "SizeMin", self.tailLEcellSize)
        gmsh.model.mesh.field.setNumber(8, "SizeMax", self.cellSize)
        gmsh.model.mesh.field.setNumber(8, "DistMin", tailMinDist)
        gmsh.model.mesh.field.setNumber(8, "DistMax", tailMaxDist)

        # Tail trailing edge
        gmsh.model.mesh.field.add("Distance", 9)
        gmsh.model.mesh.field.setNumbers(9, "CurvesList", [self.curve82e, self.curve83e, self.curve83d, self.curve82d])
        gmsh.model.mesh.field.setNumber(9, "NumPointsPerCurve", 200)

        gmsh.model.mesh.field.add("Threshold", 10)
        gmsh.model.mesh.field.setNumber(10, "InField", 9)
        gmsh.model.mesh.field.setNumber(10, "SizeMin", self.tailTEcellSize)
        gmsh.model.mesh.field.setNumber(10, "SizeMax", self.cellSize)
        gmsh.model.mesh.field.setNumber(10, "DistMin", tailMinDist)
        gmsh.model.mesh.field.setNumber(10, "DistMax", tailMaxDist)

        gmsh.model.mesh.field.add("Min", 11)
        gmsh.model.mesh.field.setNumbers(11, "FieldsList", [2, 4, 6, 8, 10])

        gmsh.model.mesh.field.setAsBackgroundMesh(11)

        gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)

        return