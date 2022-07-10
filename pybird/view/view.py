from pybird.mesh.mesh import Mesh
from pybird.view.utils.opengl_ui import Geometry, OpenGL_UI, Wake

class View:

    def __init__(self, mesh: Mesh) -> None:
        self._mesh = mesh
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