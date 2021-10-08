import Sofa

def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', pluginName='SoftRobots')
    rootNode.addObject('RequiredPlugin', pluginName='BeamAdapter')

    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('LCPConstraintSolver', mu='0.1', tolerance='1e-10', maxIt='1000', build_lcp='false')
    rootNode.addObject('BruteForceBroadPhase', name='N2')
    rootNode.addObject('BVHNarrowPhase')
    rootNode.addObject('LocalMinDistance', contactDistance='1', alarmDistance='3', name='localmindistance', angleCone='0.02')

    topoLines_guide = rootNode.addChild('topoLines_guide')
    topoLines_guide.addObject('WireRestShape', name='GuideRestShape', 
                                 straightLength=980.0, length=1000.0, 
                                 numEdges=200, youngModulus=20000, 
                                 spireDiameter=25, numEdgesCollis=[50,10], 
                                 printLog=True, template='Rigid3d', spireHeight=0.0, 
                                 densityOfBeams=[30,5], youngModulusExtremity=20000)
    topoLines_guide.addObject('EdgeSetTopologyContainer', name='meshLinesGuide')
    topoLines_guide.addObject('EdgeSetTopologyModifier', name='Modifier')
    topoLines_guide.addObject('MechanicalObject', name='dofTopo2', template='Rigid3d')


    InstrumentCombined = rootNode.addChild('InstrumentCombined')
    InstrumentCombined.addObject('EulerImplicit', rayleighStiffness=0.2, 
                                    printLog=False, rayleighMass=0.1)
    InstrumentCombined.addObject('BTDLinearSolver', verification=False, 
                                    subpartSolve=False, verbose=False)
    InstrumentCombined.addObject('RegularGrid', name='meshLinesCombined', 
                                    zmax=1, zmin=1, nx=60, ny=1, nz=1, 
                                    xmax=1.0, xmin=0.0, ymin=0, ymax=0)
    InstrumentCombined.addObject('MechanicalObject', showIndices=False, name='DOFs', template='Rigid3d', ry=-90)
    InstrumentCombined.addObject('WireBeamInterpolation', WireRestShape='@../topoLines_guide/GuideRestShape', 
                                    radius=0.9, printLog=True, name='InterpolGuide')
    InstrumentCombined.addObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155, 
                                    name='GuideForceField', interpolation='@InterpolGuide')

    InstrumentCombined.addObject('InterventionalRadiologyController', xtip=[1, 0, 0], name='m_ircontroller', 
                                    instruments='InterpolGuide', 
                                    step=0.5, printLog=True, 
                                    listening=True, template='Rigid3d', startingPos=[0, 0, 0, 1, 0, 0, 0], 
                                    rotationInstrument=[0, 0, 0], speed=0, 
                                    controlledInstrument=0)
    InstrumentCombined.addObject('LinearSolverConstraintCorrection', wire_optimization='true', printLog=False)
    InstrumentCombined.addObject('FixedConstraint', indices=0, name='FixedConstraint')
    InstrumentCombined.addObject('RestShapeSpringsForceField', points='@m_ircontroller.indexFirstNode', 
                                    angularStiffness=1e8, stiffness=1e8)
    
