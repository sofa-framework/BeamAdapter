# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:28:17 2019

@author: PSC
"""

def createGuide(node, name, listening=False, straightLength=980.0, length=1000.0, 
                numEdges=200, youngModulus=20000, spireDiameter=25, 
                numEdgesCollis=[50,10], spireHeight=0.0, densityOfBeams=[30,5], youngModulusExtremity=10000):
    topoLines_guide = node.createChild('topoLines_'+name)
    topoLines_guide.createObject('WireRestShape', name=name+'RestShape', 
                                 straightLength=straightLength, length=length, 
                                 numEdges=numEdges, youngModulus=youngModulus, 
                                 spireDiameter=spireDiameter, numEdgesCollis=numEdgesCollis, 
                                 printLog=True, template='Rigid3d', spireHeight=spireHeight, 
                                 densityOfBeams=densityOfBeams, youngModulusExtremity=youngModulusExtremity)
    topoLines_guide.createObject('EdgeSetTopologyContainer', name='meshLinesGuide')
    topoLines_guide.createObject('EdgeSetTopologyModifier', name='Modifier')
    topoLines_guide.createObject('EdgeSetTopologyAlgorithms', name='TopoAlgo', template='Rigid3d')
    topoLines_guide.createObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid3d')
    topoLines_guide.createObject('MechanicalObject', name='dofTopo2', template='Rigid3d')
    
    return(topoLines_guide)
    
def createInstrumentsCombined(node, xtip=[1, 0, 0], instruments=['guide'], 
                              step=0.5, listening=True, startingPos=[0, 0, 0, 1, 0, 0, 0], 
                              rotationInstrument=[0, 0, 0], speed=0, controlledInstrument=0):
    InstrumentCombined = node.createChild('InstrumentCombined')
    InstrumentCombined.createObject('EulerImplicit', rayleighStiffness=0.2, 
                                    printLog=False, rayleighMass=0.1)
    InstrumentCombined.createObject('BTDLinearSolver', verification=False, 
                                    subpartSolve=False, verbose=False)
    InstrumentCombined.createObject('RegularGrid', name='meshLinesCombined', 
                                    zmax=1, zmin=1, nx=60, ny=1, nz=1, 
                                    xmax=1.0, xmin=0.0, ymin=0, ymax=0)
    InstrumentCombined.createObject('MechanicalObject', showIndices=False, name='DOFs', template='Rigid3d', ry=-90)
    for i in range(len(instruments)):
        InstrumentCombined.createObject('WireBeamInterpolation', WireRestShape='@../topoLines_'+instruments[i]+'/'+instruments[i]+'RestShape', 
                                    radius=0.15, printLog=True, name='Interpol'+instruments[i])
        InstrumentCombined.createObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155, 
                                    name=instruments[i]+'ForceField', interpolation='@Interpol'+instruments[i])
    # InstrumentCombined.createObject('WireBeamInterpolation', WireRestShape='@../topoLines_guide/guideRestShape', 
    #                                 radius=0.15, printLog=True, name='InterpolGuide')
    # InstrumentCombined.createObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155, 
    #                                 name='GuideForceField', interpolation='@InterpolGuide')
    InstrumentCombined.createObject('InterventionalRadiologyController', xtip=xtip, name='m_ircontroller', 
                                    instruments=['Interpol'+instruments[i] for i in range(len(instruments))], 
                                    step=step, printLog=True, 
                                    listening=listening, template='Rigid3d', startingPos=startingPos, 
                                    rotationInstrument=rotationInstrument, speed=speed, 
                                    controlledInstrument=controlledInstrument)
    InstrumentCombined.createObject('LinearSolverConstraintCorrection', wire_optimization='true', printLog=False)
    InstrumentCombined.createObject('FixedConstraint', indices=0, name='FixedConstraint')
    InstrumentCombined.createObject('RestShapeSpringsForceField', points='@m_ircontroller.indexFirstNode', 
                                    angularStiffness=1e8, stiffness=1e8)

    # Collision model
    Collis = InstrumentCombined.createChild('Collis')
    Collis.activated = 'true'
    Collis.createObject('EdgeSetTopologyContainer', name='collisEdgeSet')
    Collis.createObject('EdgeSetTopologyModifier', name='colliseEdgeModifier')
    Collis.createObject('MechanicalObject', name='CollisionDOFs')
    Collis.createObject('MultiAdaptiveBeamMapping', controller='../m_ircontroller', 
                        useCurvAbs=True, printLog=False, name='collisMap')
    Collis.createObject('Line', proximity=0.0, group=1)
    Collis.createObject('Point', proximity=0.0, group=1)

    # Visualization Guide
    VisuGuide = InstrumentCombined.createChild('VisuGuide')
    VisuGuide.activated = 'true'
    VisuGuide.createObject('MechanicalObject', name='Quads')
    VisuGuide.createObject('QuadSetTopologyContainer', name='ContainerGuide')
    VisuGuide.createObject('QuadSetTopologyModifier', name='Modifier')
    VisuGuide.createObject('QuadSetTopologyAlgorithms', name='TopoAlgo', template='Vec3d')
    VisuGuide.createObject('QuadSetGeometryAlgorithms', name='GeomAlgo', template='Vec3d')
    VisuGuide.createObject('Edge2QuadTopologicalMapping', radius='1', listening='true', 
                           input='@../../topoLines_guide/meshLinesGuide', 
                           nbPointsOnEachCircle='10', flipNormals='true', output='@ContainerGuide')
    VisuGuide.createObject('AdaptiveBeamMapping', interpolation='@../InterpolGuide', 
                           name='visuMapGuide', output='@Quads', isMechanical='false', 
                           input='@../DOFs', useCurvAbs=True, printLog=True)

    # Ogl model
    VisuOgl = VisuGuide.createChild('VisuOgl')
    VisuOgl.createObject('OglModel', color=[0.2, 0.2, 0.8], quads='@../ContainerGuide.quads', 
                         material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20', name='Visual')
    VisuOgl.createObject('IdentityMapping', input='@../Quads', output='@Visual')
    
    return (InstrumentCombined)