# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:21:48 2019

@author: PSC
"""

def createGeometry(node, stl, scale=1, rotation=[0.0, 0.0, 0.0], VISUAL = True):
    Geometry = node.addChild('Vessels')
    Geometry.addObject('MeshSTLLoader', filename=stl, flipNormals=False, triangulate=True, name='meshLoader', scale=scale, rotation=rotation)
    Geometry.addObject('MeshTopology', position='@meshLoader.position', triangles='@meshLoader.triangles')
    Geometry.addObject('MechanicalObject', name='DOFs1', scale=1, rotation=rotation)
    Geometry.addObject('TriangleCollisionModel', moving=False, simulated=False)
    Geometry.addObject('LineCollisionModel', moving=False, simulated=False)
    Geometry.addObject('PointCollisionModel', moving=False, simulated=False)

    ### This feature is broken
    # if VISUAL :
    # 	Geometry.addObject('OglModel', color=[1, 0, 0, 0.3], src='@meshLoader', name='Visual')

    
    return(Geometry)