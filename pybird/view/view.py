from typing import List
import matplotlib.cm as cm

from numpy import ndarray, zeros
from pybird.mesh.mesh import Mesh
from pybird.view.utils.opengl_ui import Geometry, OpenGL_UI, OpenGL_UI_Vector_Field, Wake

class View:

    def __init__(self, mesh: Mesh) -> None:
        self._mesh = mesh
        return
    
    def surfaceParameter(self, values: ndarray) -> None:

        n = len(self._mesh.vertices[:, 0])
        verticesSum = zeros(n)
        verticesNumbs = zeros(n)

        for i in range(len(self._mesh.faces[:, 0])):
            face = self._mesh.faces[i, :]
            verticesSum[face[0]] = values[i]
            verticesSum[face[1]] = values[i]
            verticesSum[face[2]] = values[i]
            verticesNumbs[face[0]] = verticesNumbs[face[0]] + 1
            verticesNumbs[face[1]] = verticesNumbs[face[1]] + 1
            verticesNumbs[face[2]] = verticesNumbs[face[2]] + 1

        verticesValues = verticesSum / (verticesNumbs + 1e-8)

        viridis = cm.get_cmap('gnuplot', 8)
        newValues = verticesValues - min(verticesValues)
        newValues = newValues / (max(newValues) - min(newValues))
        colors = viridis(newValues)

        opengl = OpenGL_UI(
            geos=[
                Geometry(
                    vertices=self._mesh.vertices,
                    edges=self._mesh.edges,
                    faces=self._mesh.faces,
                ),
            ],
            colors=colors,
        )
        opengl.ui()

        return
    
    def faceVectors(self, vectors: List[ndarray]) -> None:
        opengl = OpenGL_UI_Vector_Field(
            geo=Geometry(
                    vertices=self._mesh.vertices,
                    edges=self._mesh.edges,
                    faces=self._mesh.faces,
            ),
            parts=[
                self._mesh.leftWingFirstSectionFacesTags,
                self._mesh.leftWingSecondSectionFacesTags,
                self._mesh.leftWingThirdSectionFacesTags,
                self._mesh.rightWingFirstSectionFacesTags,
                self._mesh.rightWingSecondSectionFacesTags,
                self._mesh.rightWingThirdSectionFacesTags,
                self._mesh.bodyFacesTags,
                self._mesh.headFacesTags,
                self._mesh.tailFacesTags,
            ],
            vectors=vectors,
            positions=self._mesh.facesCenters,
        )
        opengl.ui()
        return

    def mesh(self) -> None:
        opengl = OpenGL_UI(
            geos=[
                Geometry(
                    vertices=self._mesh.vertices,
                    edges=self._mesh.edges,
                    faces=self._mesh.faces,
                ),
            ],
            wakes=[
                Wake(
                    vertices=self._mesh.wake.leftWing.vertices,
                    ids=self._mesh.wake.leftWing.grid,
                ),
                Wake(
                    vertices=self._mesh.wake.rightWing.vertices,
                    ids=self._mesh.wake.rightWing.grid,
                ),
                Wake(
                    vertices=self._mesh.wake.tail.vertices,
                    ids=self._mesh.wake.tail.grid,
                ),
            ],
            parts=[
                self._mesh.leftWingFirstSectionFacesTags,
                self._mesh.leftWingSecondSectionFacesTags,
                self._mesh.leftWingThirdSectionFacesTags,
                self._mesh.rightWingFirstSectionFacesTags,
                self._mesh.rightWingSecondSectionFacesTags,
                self._mesh.rightWingThirdSectionFacesTags,
                self._mesh.bodyFacesTags,
                self._mesh.headFacesTags,
                self._mesh.tailFacesTags,
            ]
        )
        opengl.ui()
        return