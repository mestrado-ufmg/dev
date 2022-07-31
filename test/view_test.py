import os
import sys

sys.path.append('./')

import pybird

import numpy as np
import pygame
from OpenGL.GL import *
from OpenGL.GLU import *
from ctypes import *

if __name__ == '__main__':

    # Create model
    model = pybird.model('test')

    # Load data
    model.load('./data/data.case')

    # Create geometry
    model.geo.build()

    # Create mesh
    model.mesh.build(
        size=0.06,
        accom_dist=10,
        n_head=3,
        n_wing_le=3,
        n_wing_te=3,
        n_tail_le=3,
        n_tail_te=3,
    )

    # View
    pygame.init ()
    screen = pygame.display.set_mode ((400,300), pygame.OPENGL | pygame.DOUBLEBUF | pygame.RESIZABLE, 24)
    glViewport (0, 0, 400,300)
    glClearColor (0.0, 0.5, 0.5, 1.0)
    # os.environ['SDL_VIDEO_CENTERED'] = '1'
    # windowSize = (1000, 500)
    # # pygame.init ()
    # pygame.display.set_mode(windowSize, pygame.DOUBLEBUF | pygame.OPENGL | pygame.RESIZABLE)
    # pygame.display.set_caption('pyMESH')
    # # gluPerspective(45, (windowSize[0] / windowSize[1]), 0.1, 1000)

    # # glTranslatef(1, 0, 0) # camera positions
    # # self.body_system.rotZ(-45)
    # # glRotatef(-45, self.body_system.e3[0], self.body_system.e3[1], self.body_system.e3[2])
    # # self.body_system.rotX(-45)
    # # glRotatef(-45, self.body_system.e1[0], self.body_system.e1[1], self.body_system.e1[2])

    # glEnable(GL_CULL_FACE)
    # glEnable(GL_DEPTH_TEST)
    # glCullFace(GL_BACK)
    # glClearColor(0.2, 0.85, 0.85, 1.0)

    glEnableClientState (GL_VERTEX_ARRAY)

    # Create vertices array
    n = len(model.mesh.faces[:, 0])
    vertices1 = model.mesh.vertices[model.mesh.faces[:, 0], :].astype(np.float32)
    vertices2 = model.mesh.vertices[model.mesh.faces[:, 1], :].astype(np.float32)
    vertices3 = model.mesh.vertices[model.mesh.faces[:, 2], :].astype(np.float32)
    vertices = []
    for i in range(n):
        vertices.append(vertices1[i, 0])
        vertices.append(vertices1[i, 1])
        vertices.append(vertices1[i, 2])
        vertices.append(vertices2[i, 0])
        vertices.append(vertices2[i, 1])
        vertices.append(vertices2[i, 2])
        vertices.append(vertices3[i, 0])
        vertices.append(vertices3[i, 1])
        vertices.append(vertices3[i, 2])

    # vertices = [ 0.0, 1.0, 0.0,  0.0, 0.0, 0.0,  1.0, 1.0, 0.0 ]

    vbo = glGenBuffers(1)
    glBindBuffer(GL_ARRAY_BUFFER, vbo)
    glBufferData(GL_ARRAY_BUFFER, len(vertices)*4, (c_float*len(vertices))(*vertices), GL_STATIC_DRAW)

    running = True
    while running:

        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False

        glClear (GL_COLOR_BUFFER_BIT)

        glBindBuffer (GL_ARRAY_BUFFER, vbo)
        glVertexPointer (3, GL_FLOAT, 0, None)

        glDrawArrays (GL_TRIANGLES, 0, 3)

        pygame.display.flip ()