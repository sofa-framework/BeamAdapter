# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:35:00 2019

@author: PSC
"""

import Sofa

from beamadapter.parts import createGuide, createInstrumentsCombined, createGeometry
from beamadapter.sensors import CollisionMonitor


def createScene(rootNode):

    rootNode.createObject('RequiredPlugin', pluginName='SoftRobots')
    rootNode.createObject('RequiredPlugin', pluginName='BeamAdapter')
    # AnimationManager(rootNode)    
    
    rootNode.createObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings hideForceFields')
    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('LCPConstraintSolver', mu='0.1', tolerance='1e-10', maxIt='1000', build_lcp='false')
    rootNode.createObject('CollisionPipeline', draw='0', depth='6', verbose='1')
    rootNode.createObject('BruteForceDetection', name='N2')
    rootNode.createObject('LocalMinDistance', contactDistance='1', alarmDistance='3', name='localmindistance', angleCone='0.02')
    rootNode.createObject('CollisionResponse', name='Response', response='FrictionContact')
    rootNode.createObject('CollisionGroup', name='Group')
    
    topoLines_guide = createGuide(rootNode, straightLength=980.0, length=1000.0, 
                youngModulus=20000, numEdges=200, spireDiameter=25, 
                numEdgesCollis=[50,10], spireHeight=0.0, densityOfBeams=[30,5], 
                youngModulusExtremity=1000)
    instrumentsCombined = createInstrumentsCombined(rootNode)
    
    vessels = createGeometry(rootNode, 'mesh/FlatLCA1_5.stl', scale=3, rotation=[0.0, 5.0, 0.0] , VISUAL = False)

    CollisionMonitor(rootNode, instrumentsCombined.DOFs, 60, 6, verbose= True )
    CollisionMonitor(rootNode, vessels.DOFs1, len(vessels.meshLoader.position), 3)
    
