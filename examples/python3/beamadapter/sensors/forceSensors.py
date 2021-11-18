# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:35:00 2019

@author: PSC
"""

import Sofa
import numpy as np
import math

from OpenGL.GL import *
from OpenGL.GLU import *


class CollisionMonitor( Sofa.PythonScriptController ):
    def __init__(self, node, MO, numberNodes, numberDOFs, verbose = False ):
        self.MO = MO
        self.numberNodes = numberNodes
        self.numberDOFs = numberDOFs
        self.verbose = verbose
        print('init collision monitor')

    def onBeginAnimationStep(self, dt):

        self.sortedCollisionMatrix = [[0, 0, 0] for i in range (self.numberNodes)]

        self.collisionMatrix = self.MO.constraint.splitlines()
        self.collisionMatrix = [line.split() for line in self.collisionMatrix]
        self.collisionMatrix = [[float(Value) for Value in self.collisionMatrix[i]] for i in range(len(self.collisionMatrix))]

        for i in range(len(self.collisionMatrix)):
          numberElements = self.collisionMatrix[i][1]

          if numberElements == 1:
            element1 = int(self.collisionMatrix[i][2])
            self.sortedCollisionMatrix[element1] = [ self.sortedCollisionMatrix[element1][j] + self.collisionMatrix[i][3+j] for j in range(3)] 

          if numberElements == 2:
            element1 = int(self.collisionMatrix[i][2])
            element2 = int(self.collisionMatrix[i][3 + self.numberDOFs ])
            self.sortedCollisionMatrix[element1] = [ self.sortedCollisionMatrix[element1][j] + self.collisionMatrix[i][3+j] for j in range(3)] 
            self.sortedCollisionMatrix[element2] = [ self.sortedCollisionMatrix[element2][j] + self.collisionMatrix[i][4+ self.numberDOFs +j] for j in range(3)] 

          if numberElements == 3:
            element1 = int(self.collisionMatrix[i][2])
            element2 = int(self.collisionMatrix[i][3+ self.numberDOFs ])
            element3 = int(self.collisionMatrix[i][4 + 2*self.numberDOFs ])
            self.sortedCollisionMatrix[element1] = [ self.sortedCollisionMatrix[element1][j] + self.collisionMatrix[i][3+j] for j in range(3)] 
            self.sortedCollisionMatrix[element2] = [ self.sortedCollisionMatrix[element2][j] + self.collisionMatrix[i][4+ self.numberDOFs +j] for j in range(3)] 
            self.sortedCollisionMatrix[element3] = [ self.sortedCollisionMatrix[element3][j] + self.collisionMatrix[i][5+ 2*self.numberDOFs +j] for j in range(3)] 

        for i in range(len(self.sortedCollisionMatrix)):
            self.sortedCollisionMatrix[i] = np.linalg.norm( np.array(self.sortedCollisionMatrix[i]) )

        if(self.verbose):
            print(self.MO, self.sortedCollisionMatrix)


def eulerAnglesToRotationMatrix(theta) :
    theta = [-theta[0]/180 * math.pi, -theta[1]/180 * math.pi, -theta[2]/180 * math.pi]
     
    R_x = np.array([[1,         0,                  0                   ],
                    [0,         math.cos(theta[0]), -math.sin(theta[0]) ],
                    [0,         math.sin(theta[0]), math.cos(theta[0])  ]
                    ])
                        
    R_y = np.array([[math.cos(theta[1]),    0,      math.sin(theta[1])  ],
                    [0,                     1,      0                   ],
                    [-math.sin(theta[1]),   0,      math.cos(theta[1])  ]
                    ])
                 
    R_z = np.array([[math.cos(theta[2]),    -math.sin(theta[2]),    0],
                    [math.sin(theta[2]),    math.cos(theta[2]),     0],
                    [0,                     0,                      1]
                    ])
                      
    R = np.dot(R_z, np.dot( R_y, R_x ))
 
    return R

def mapColors(val, maxval):
    if val > maxval :
        red = 1
    else :
        red = 1 - (maxval - val)/ maxval
    green = 1 - red
    return ([red, green])

class drawCollisionMonitor( Sofa.PythonScriptController ):
    def __init__(self, node, mesh, MO, numberNodes, numberDOFs, verbose = False ):
        self.mesh = mesh
        self.MO = MO
        self.numberNodes = numberNodes
        self.numberDOFs = numberDOFs
        self.sortedCollisionMatrix = [[0, 0, 0] for i in range (self.numberNodes)]
        self.verbose = verbose
        # print('init collision monitor')

    def onBeginAnimationStep(self, dt):

        self.sortedCollisionMatrix = [[0, 0, 0] for i in range (self.numberNodes)]

        self.collisionMatrix = self.MO.constraint.splitlines()
        self.collisionMatrix = [line.split() for line in self.collisionMatrix]
        self.collisionMatrix = [[float(Value) for Value in self.collisionMatrix[i]] for i in range(len(self.collisionMatrix))]

        for i in range(len(self.collisionMatrix)):
          numberElements = self.collisionMatrix[i][1]

          if numberElements == 1:
            element1 = int(self.collisionMatrix[i][2])
            self.sortedCollisionMatrix[element1] = [ self.sortedCollisionMatrix[element1][j] + self.collisionMatrix[i][3+j] for j in range(3)] 

          if numberElements == 2:
            element1 = int(self.collisionMatrix[i][2])
            element2 = int(self.collisionMatrix[i][3 + self.numberDOFs ])
            self.sortedCollisionMatrix[element1] = [ self.sortedCollisionMatrix[element1][j] + self.collisionMatrix[i][3+j] for j in range(3)] 
            self.sortedCollisionMatrix[element2] = [ self.sortedCollisionMatrix[element2][j] + self.collisionMatrix[i][4+ self.numberDOFs +j] for j in range(3)] 

          if numberElements == 3:
            element1 = int(self.collisionMatrix[i][2])
            element2 = int(self.collisionMatrix[i][3+ self.numberDOFs ])
            element3 = int(self.collisionMatrix[i][4 + 2*self.numberDOFs ])
            self.sortedCollisionMatrix[element1] = [ self.sortedCollisionMatrix[element1][j] + self.collisionMatrix[i][3+j] for j in range(3)] 
            self.sortedCollisionMatrix[element2] = [ self.sortedCollisionMatrix[element2][j] + self.collisionMatrix[i][4+ self.numberDOFs +j] for j in range(3)] 
            self.sortedCollisionMatrix[element3] = [ self.sortedCollisionMatrix[element3][j] + self.collisionMatrix[i][5+ 2*self.numberDOFs +j] for j in range(3)] 

        for i in range(len(self.sortedCollisionMatrix)):
            self.sortedCollisionMatrix[i] = np.linalg.norm( np.array(self.sortedCollisionMatrix[i]) )

        if(self.verbose):
            print(self.MO, self.sortedCollisionMatrix)

    def draw(self):

        scale = np.array([3, 3, 3])
        R = eulerAnglesToRotationMatrix([0.0, 5.0, 0.0])

        self.verticies = scale * np.array(self.mesh.position)
        self.verticies = np.dot(self.verticies, R)

        self.triangles = self.mesh.triangles

        glBegin(GL_TRIANGLES)
        for triangle in self.triangles:
            for vertex in triangle:

                if self.sortedCollisionMatrix[vertex] > 0.5:
                    glColor3f(mapColors(self.sortedCollisionMatrix[vertex], 10)[0], mapColors(self.sortedCollisionMatrix[vertex], 10)[1], 0)
                    p = self.verticies[vertex]
                    glVertex3f(p[0], p[1], p[2])


        glEnd()