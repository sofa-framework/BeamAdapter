# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 14:28:17 2019

@author: PSC
"""

def createGuide(node, name, listening=False, straightLength=980.0, length=1000.0, 
                numEdges=200, youngModulus=20000, spireDiameter=25, 
                numEdgesCollis=[50,10], spireHeight=0.0, densityOfBeams=[30,5], youngModulusExtremity=10000):
    topoLines_guide = node.addChild('topoLines_'+name)
    topoLines_guide.addObject('WireRestShape', name=name+'RestShape', 
                             straightLength=straightLength, length=length, 
                             numEdges=numEdges, youngModulus=youngModulus, 
                             spireDiameter=spireDiameter, numEdgesCollis=numEdgesCollis, 
                             printLog=True, template='Rigid3d', spireHeight=spireHeight, 
                             densityOfBeams=densityOfBeams, youngModulusExtremity=youngModulusExtremity)

    topoLines_guide.addObject('MechanicalObject', name='dofTopo2', template='Rigid3d')

    topoLines_guide.addObject('EdgeSetTopologyContainer', name='meshLinesGuide')
    topoLines_guide.addObject('EdgeSetTopologyModifier', name='Modifier')
    topoLines_guide.addObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid3d')

    
    
    return(topoLines_guide)

def combineInstruments(node, xtip=[1, 0, 0], instruments=['guide'], 
                              step=0.5, listening=True, startingPos=[0, 0, 0, 1, 0, 0, 0], 
                              rotationInstrument=[0, 0, 0], speed=0, controlledInstrument=0):

    InstrumentCombined = node.addChild('InstrumentCombined')
    InstrumentCombined.addObject('EulerImplicitSolver', rayleighStiffness=0.2, 
                                    printLog=False, rayleighMass=0.1)
    InstrumentCombined.addObject('BTDLinearSolver', verification=False, 
                                    subpartSolve=False, verbose=False)
    InstrumentCombined.addObject('RegularGridTopology', name='meshLinesCombined', 
                                    zmax=1, zmin=1, nx=60, ny=1, nz=1, 
                                    xmax=1.0, xmin=0.0, ymin=0, ymax=0)
    InstrumentCombined.addObject('MechanicalObject', showIndices=False, name='DOFs', template='Rigid3d', ry=-90)
    for i in range(len(instruments)):
        InstrumentCombined.addObject('WireBeamInterpolation', WireRestShape='@../topoLines_'+instruments[i]+'/'+instruments[i]+'RestShape', 
                                    radius=0.15, printLog=False, name='Interpol'+instruments[i])
        InstrumentCombined.addObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155, 
                                    name=instruments[i]+'ForceField', interpolation='@Interpol'+instruments[i])

    InstrumentCombined.addObject('InterventionalRadiologyController', xtip=xtip, name='m_ircontroller', 
                                    instruments=['Interpol'+instruments[i] for i in range(len(instruments))], 
                                    step=step, printLog=True, 
                                    listening=listening, template='Rigid3d', startingPos=startingPos, 
                                    rotationInstrument=rotationInstrument, speed=speed, 
                                    controlledInstrument=controlledInstrument)
    InstrumentCombined.addObject('LinearSolverConstraintCorrection', wire_optimization='true', printLog=False)
    InstrumentCombined.addObject('FixedConstraint', indices=0, name='FixedConstraint')
    InstrumentCombined.addObject('RestShapeSpringsForceField', points='@m_ircontroller.indexFirstNode', 
                                    angularStiffness=1e8, stiffness=1e8)

    return (InstrumentCombined)



    
def createInstrumentsCombined(node, xtip=[1, 0, 0], instruments=['guide'], 
                              step=0.5, listening=True, startingPos=[0, 0, 0, 1, 0, 0, 0], 
                              rotationInstrument=[0, 0, 0], speed=0, controlledInstrument=0):
    # InstrumentCombined = combineInstruments(node, xtip, instruments, step, listening, startingPos, 
    #                           rotationInstrument, speed, controlledInstrument)

    InstrumentCombined = node.addChild('InstrumentCombined')
    InstrumentCombined.addObject('EulerImplicitSolver', rayleighStiffness=0.2, 
                                    printLog=False, rayleighMass=0.1)
    InstrumentCombined.addObject('BTDLinearSolver', verification=False, 
                                    subpartSolve=False, verbose=False)
    InstrumentCombined.addObject('RegularGridTopology', name='meshLinesCombined', 
                                    zmax=1, zmin=1, nx=60, ny=1, nz=1, 
                                    xmax=1.0, xmin=0.0, ymin=0, ymax=0)
    InstrumentCombined.addObject('MechanicalObject', showIndices=False, name='DOFs', template='Rigid3d', ry=-90)
    for i in range(len(instruments)):
        InstrumentCombined.addObject('WireBeamInterpolation', WireRestShape='@../topoLines_'+instruments[i]+'/'+instruments[i]+'RestShape', 
                                    radius=0.15, printLog=False, name='Interpol'+instruments[i])
        InstrumentCombined.addObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155, 
                                    name=instruments[i]+'ForceField', interpolation='@Interpol'+instruments[i])

    InstrumentCombined.addObject('InterventionalRadiologyController', xtip=xtip, name='m_ircontroller', 
                                    instruments=['Interpol'+instruments[i] for i in range(len(instruments))], 
                                    step=step, printLog=True, 
                                    listening=listening, template='Rigid3d', startingPos=startingPos, 
                                    rotationInstrument=rotationInstrument, speed=speed, 
                                    controlledInstrument=controlledInstrument)
    InstrumentCombined.addObject('LinearSolverConstraintCorrection', wire_optimization='true', printLog=False)
    InstrumentCombined.addObject('FixedConstraint', indices=0, name='FixedConstraint')
    InstrumentCombined.addObject('RestShapeSpringsForceField', points='@m_ircontroller.indexFirstNode', 
                                    angularStiffness=1e8, stiffness=1e8)

    # Collision model
    Collis = InstrumentCombined.addChild('Collis')
    Collis.activated = True
    Collis.addObject('EdgeSetTopologyContainer', name='collisEdgeSet')
    Collis.addObject('EdgeSetTopologyModifier', name='colliseEdgeModifier')
    Collis.addObject('MechanicalObject', name='CollisionDOFs')
    Collis.addObject('MultiAdaptiveBeamMapping', controller='../m_ircontroller', 
                        useCurvAbs=True, printLog=False, name='collisMap')
    Collis.addObject('LineCollisionModel', proximity=0.0, group=1)
    Collis.addObject('PointCollisionModel', proximity=0.0, group=1)

    # Visualization Guide
    VisuGuide = InstrumentCombined.addChild('VisuGuide')
    VisuGuide.activated = True
    VisuGuide.addObject('MechanicalObject', name='Quads')
    VisuGuide.addObject('QuadSetTopologyContainer', name='ContainerGuide')
    VisuGuide.addObject('QuadSetTopologyModifier', name='Modifier')
    VisuGuide.addObject('QuadSetGeometryAlgorithms', name='GeomAlgo', template='Vec3d')
    VisuGuide.addObject('Edge2QuadTopologicalMapping', radius='1', listening='true', 
                           input='@../../topoLines_guide/meshLinesGuide', 
                           nbPointsOnEachCircle='10', flipNormals='true', output='@ContainerGuide')
    VisuGuide.addObject('AdaptiveBeamMapping', interpolation='@../InterpolGuide', 
                           name='visuMapGuide', output='@Quads', isMechanical=False, 
                           input='@../DOFs', useCurvAbs=True, printLog=False)

    # Ogl model
    VisuOgl = VisuGuide.addChild('VisuOgl')
    VisuOgl.addObject('OglModel', color=[0.2, 0.2, 0.8], quads='@../ContainerGuide.quads', 
                         material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20', name='Visual')
    VisuOgl.addObject('IdentityMapping', input='@../Quads', output='@Visual')
    
    return (InstrumentCombined)


def createInstrumentsCombinedXRay(node, xtip=[1, 0, 0], instruments=['guide'], 
                              step=0.5, listening=True, startingPos=[0, 0, 0, 1, 0, 0, 0], 
                              rotationInstrument=[0, 0, 0], speed=0, controlledInstrument=0):
    
    InstrumentCombined = combineInstruments(node, xtip, instruments, step, listening, startingPos, 
                              rotationInstrument, speed, controlledInstrument)

    # Collision model
    Collis = InstrumentCombined.addChild('Collis')
    #Collis.activated = 'true'
    Collis.addObject('EdgeSetTopologyContainer', name='collisEdgeSet')
    Collis.addObject('EdgeSetTopologyModifier', name='colliseEdgeModifier')
    Collis.addObject('MechanicalObject', name='CollisionDOFs')
    Collis.addObject('MultiAdaptiveBeamMapping', controller='../m_ircontroller', 
                        useCurvAbs=True, printLog=False, name='collisMap')
    Collis.addObject('LineCollisionModel', proximity=0.0, group=1)
    Collis.addObject('PointCollisionModel', proximity=0.0, group=1)

    # Visualization Guide
    VisuGuide = InstrumentCombined.addChild('VisuGuide')
    #VisuGuide.activated = 'true'
    VisuGuide.addObject('MechanicalObject', name='Quads')
    VisuGuide.addObject('QuadSetTopologyContainer', name='ContainerGuide')
    VisuGuide.addObject('QuadSetTopologyModifier', name='Modifier')
    VisuGuide.addObject('QuadSetGeometryAlgorithms', name='GeomAlgo', template='Vec3d')
    VisuGuide.addObject('Edge2QuadTopologicalMapping', radius='1', listening=True, 
                           input='@../../topoLines_guide/meshLinesGuide', 
                           nbPointsOnEachCircle='10', flipNormals='true', output='@ContainerGuide')
    VisuGuide.addObject('AdaptiveBeamMapping', interpolation='@../InterpolGuide', 
                           name='visuMapGuide', output='@Quads', isMechanical=False, 
                           input='@../DOFs', useCurvAbs=True, printLog=True)

    TriangleTopology = VisuGuide.addChild('TriangleTopology')
    TriangleTopology.addObject('TriangleSetTopologyContainer', name='TriangleContainer')
    TriangleTopology.addObject('TriangleSetTopologyModifier', name='Modifier')
    TriangleTopology.addObject('TriangleSetGeometryAlgorithms', name='GeomAlgo', template='Vec3d')
    TriangleTopology.addObject('Quad2TriangleTopologicalMapping', input='@../ContainerGuide', output='@TriangleContainer')

    # Ogl model
    VisuOgl = VisuGuide.addChild('VisuOgl')
    VisuOgl.addObject('OglModel', color=[0.2, 0.2, 0.8], triangles='@../TriangleTopology/TriangleContainer.triangles', 
                         material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20', name='Visual')
    VisuOgl.addObject('IdentityMapping', input='@../Quads', output='@Visual')
    
    return (InstrumentCombined)