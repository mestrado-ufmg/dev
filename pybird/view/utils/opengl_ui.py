from __future__ import annotations
from dataclasses import dataclass
from enum import Enum
from math import fabs

import os
os.environ["PYGAME_HIDE_SUPPORT_PROMPT"] = "hide"

import pygame
from pygame.locals import *
from OpenGL.GL import *
from OpenGL.GLU import *
from typing import List, Tuple
import numpy as np
from scipy.spatial.transform import Rotation

'''
====================================================
 Models
====================================================
'''
class Color(Enum):
    light_blue = 1
    light_orange = 2

    @staticmethod
    def rgb(value: Color) -> Tuple(int):
        if value == Color.light_blue:
            return (0.37, 0.45, 0.64)
        elif value == Color.light_orange:
            return (0.90, 0.65, 0.50)

        return (0, 0, 0)

@dataclass
class Geometry:
    '''
        Contains the information about the mesh.
    '''
    vertices: Tuple[Tuple]
    edges: Tuple[Tuple]
    faces: Tuple[Tuple]
    const_color: Color = Color.light_blue
    colors: Tuple[Tuple] = None

@dataclass
class Wake:
    '''
        Contains the information about the mesh.
    '''
    vertices: Tuple[Tuple]
    ids: Tuple[Tuple]

'''
====================================================
 Stores
====================================================
'''
class BodySystem:
    '''
        e1, e2 and e3 are the unary vactors of the body system.
    '''
    
    def __init__(self) -> None:
        self.e1 = np.array([1, 0, 0])
        self.e2 = np.array([0, 1, 0])
        self.e3 = np.array([0, 0, 1])
        return
    
    def rotX(self, angle: float) -> None:
        angle = np.deg2rad(-angle)
        r = Rotation.from_quat((np.sin(angle / 2) * self.e1).tolist() + [np.cos(angle / 2)])
        self.e2 = r.apply(self.e2)
        self.e3 = r.apply(self.e3)
        return
    
    def rotY(self, angle: float) -> None:
        angle = np.deg2rad(-angle)
        r = Rotation.from_quat((np.sin(angle / 2) * self.e2).tolist() + [np.cos(angle / 2)])
        self.e1 = r.apply(self.e1)
        self.e3 = r.apply(self.e3)
        return
    
    def rotZ(self, angle: float) -> None:
        angle = np.deg2rad(-angle)
        r = Rotation.from_quat((np.sin(angle / 2) * self.e3).tolist() + [np.cos(angle / 2)])
        self.e1 = r.apply(self.e1)
        self.e2 = r.apply(self.e2)
        return

'''
====================================================
 External class
====================================================
'''
class OpenGL_UI:
    
    def __init__(self, geos: List[Geometry] = [], wakes: List[Wake] = [], parts: List[np.ndarray] = []) -> None:

        # Inputs
        self.geos = geos
        self.wakes = wakes
        self.parts = parts


        # Stores
        self.body_system = BodySystem()

        return
    
    def _get_color(self, face: int) -> Tuple:
        for i in range(len(self.parts)):
            if face in self.parts[i]:
                if i == 0 or i == 3:
                    return (0.3, 0.3, 0.6)
                elif i == 1 or i == 4:
                    return (0.6, 0.3, 0.3)
                elif i == 2 or i == 5:
                    return (0.3, 0.6, 0.3)
                elif i == 6:
                    return (1.0, 0.4, 0.0)
                elif i == 7:
                    return (0.7, 0.7, 0.0)
                elif i == 8:
                    return (0.8, 0.1, 0.1)

        return Color.light_blue
    
    def __geo_limits(self) -> None:
        max_x, max_y, max_z = -1, -1, -1
        for geo in self.geos:
            for vertice in geo.vertices:
                if fabs(vertice[0]) > max_x: max_x = vertice[0]
                if fabs(vertice[1]) > max_y: max_y = vertice[1]
                if fabs(vertice[2]) > max_z: max_z = vertice[2]
        self.__max_x, self.__max_y, self.__max_z = max_x, max_y, max_z
        return
    
    def __draw_surface_mesh(self) -> None:

        glBegin(GL_TRIANGLES)
        for geo in self.geos:
            for i in range(len(geo.faces)):
                color = self._get_color(i) # Color.rgb(geo.const_color) if geo.colors is None else Color.rgb(geo.colors[i])
                glColor3f(color[0], color[1], color[2])
                for vertex in geo.faces[i]:
                    glVertex3fv(geo.vertices[vertex])
        glEnd()

        glBegin(GL_LINES)
        glColor3f(0, 0, 0)
        for geo in self.geos:
            for edge in geo.edges:
                for vertex in edge:
                    glVertex3fv(geo.vertices[vertex])
        
        for wake in self.wakes:
            for filament in range(len(wake.ids[:, 0])):
                for id in range(len(wake.ids[filament, :]) - 1):
                    glVertex3fv(wake.vertices[wake.ids[filament, id], :])
                    glVertex3fv(wake.vertices[wake.ids[filament, id + 1], :])

        glEnd()

        return
    
    def __draw_wake_mesh(self) -> None:
        return

    def ui(self) -> None:

        # Necessary parameters
        self.__geo_limits()

        self.__rotate_button = False
        self.__translate_button = False
        self.__zoom_button = False

        # Window
        # pygame.init()
        os.environ['SDL_VIDEO_CENTERED'] = '1'
        windowSize = (1000, 500)
        pygame.display.set_mode(windowSize, pygame.DOUBLEBUF | pygame.OPENGL | pygame.RESIZABLE)
        pygame.display.set_caption('pyMESH')
        gluPerspective(45, (windowSize[0] / windowSize[1]), 0.1, 1000)

        glTranslatef(self.__max_x, 0, - 3 * self.__max_y) # camera positions
        self.body_system.rotZ(-45)
        glRotatef(-45, self.body_system.e3[0], self.body_system.e3[1], self.body_system.e3[2])
        self.body_system.rotX(-45)
        glRotatef(-45, self.body_system.e1[0], self.body_system.e1[1], self.body_system.e1[2])

        glEnable(GL_CULL_FACE)
        glEnable(GL_DEPTH_TEST)
        glCullFace(GL_BACK)
        glClearColor(0.85, 0.85, 0.85, 1.0) # background color

        clock = pygame.time.Clock()
        running = True
        while running:
            for event in pygame.event.get():
                if event.type == pygame.QUIT:
                    running = False
                    break
                elif event.type == pygame.KEYDOWN:
                    if event.key == pygame.K_ESCAPE:
                        running = False
                        break
                elif event.type == pygame.MOUSEMOTION:
                    if self.__rotate_button == True:
                        
                        rotx = 0.2 * event.rel[1]
                        roty = 0.2 * event.rel[0]

                        self.body_system.rotX(rotx)
                        glRotatef(rotx, self.body_system.e1[0], self.body_system.e1[1], self.body_system.e1[2])
                        self.body_system.rotY(roty)
                        glRotatef(roty, self.body_system.e2[0], self.body_system.e2[1], self.body_system.e2[2])

                    if self.__translate_button == True:
                        
                        v = 0.002 * event.rel[0] * self.body_system.e1 - 0.002 * event.rel[1] * self.body_system.e2
                        glTranslatef(v[0], v[1], v[2])

                    if self.__zoom_button == True:

                        v = 0.002 * event.rel[1] * self.body_system.e3
                        glTranslatef(v[0], v[1], v[2])
            
            for event in pygame.mouse.get_pressed():
                self.__rotate_button = True if pygame.mouse.get_pressed()[0] == 1 else False
                self.__translate_button = True if pygame.mouse.get_pressed()[1] == 1 else False
                self.__zoom_button = True if pygame.mouse.get_pressed()[2] == 1 else False
            
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

            self.__draw_surface_mesh()
            self.__draw_wake_mesh()

            glFlush()

            pygame.display.flip()
            clock.tick(100)

        pygame.quit()
        quit()