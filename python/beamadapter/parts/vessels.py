# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:21:48 2019

@author: PSC
"""

def createGeometry(node, stl, scale=1, rotation=[0.0, 0.0, 0.0], VISUAL = True):
    Geometry = node.createChild('Vessels')
    Geometry.createObject('MeshSTLLoader', filename=stl, flipNormals=True, triangulate=True, name='meshLoader', scale=scale, rotation=rotation)
    Geometry.createObject('Mesh', position='@meshLoader.position', triangles='@meshLoader.triangles')
    Geometry.createObject('MechanicalObject', name='DOFs1', scale=scale, rotation=rotation)
    Geometry.createObject('Triangle', moving=False, simulated=False)
    Geometry.createObject('Line', moving=False, simulated=False)
    Geometry.createObject('Point', moving=False, simulated=False)

    ### This feature is broken
    # if VISUAL :
    # 	Geometry.createObject('OglModel', color=[1, 0, 0, 0.3], src='@meshLoader', name='Visual')

    
    return(Geometry)
