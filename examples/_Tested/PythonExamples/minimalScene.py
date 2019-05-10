import Sofa

from splib.animation import AnimationManager

def createScene(rootNode):

    rootNode.createObject('RequiredPlugin', pluginName='SoftRobots')
    rootNode.createObject('RequiredPlugin', pluginName='BeamAdapter')

    AnimationManager(rootNode)  

    rootNode.createObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('LCPConstraintSolver', mu='0.1', tolerance='1e-10', maxIt='1000', build_lcp='false')
    rootNode.createObject('BruteForceDetection', name='N2')
    rootNode.createObject('LocalMinDistance', contactDistance='1', alarmDistance='3', name='localmindistance', angleCone='0.02')

    topoLines_guide = rootNode.createChild('topoLines_guide')
    topoLines_guide.createObject('WireRestShape', name='GuideRestShape', 
                                 straightLength=980.0, length=1000.0, 
                                 numEdges=200, youngModulus=20000, 
                                 spireDiameter=25, numEdgesCollis=[50,10], 
                                 printLog=True, template='Rigid3d', spireHeight=0.0, 
                                 densityOfBeams=[30,5], youngModulusExtremity=20000)
    topoLines_guide.createObject('EdgeSetTopologyContainer', name='meshLinesGuide')
    topoLines_guide.createObject('EdgeSetTopologyModifier', name='Modifier')
    topoLines_guide.createObject('EdgeSetTopologyAlgorithms', name='TopoAlgo', template='Rigid3d')
    topoLines_guide.createObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid3d')
    topoLines_guide.createObject('MechanicalObject', name='dofTopo2', template='Rigid3d')


    InstrumentCombined = rootNode.createChild('InstrumentCombined')
    InstrumentCombined.createObject('EulerImplicit', rayleighStiffness=0.2, 
                                    printLog=False, rayleighMass=0.1)
    InstrumentCombined.createObject('BTDLinearSolver', verification=False, 
                                    subpartSolve=False, verbose=False)
    InstrumentCombined.createObject('RegularGrid', name='meshLinesCombined', 
                                    zmax=1, zmin=1, nx=60, ny=1, nz=1, 
                                    xmax=1.0, xmin=0.0, ymin=0, ymax=0)
    InstrumentCombined.createObject('MechanicalObject', showIndices=False, name='DOFs', template='Rigid3d', ry=-90)
    InstrumentCombined.createObject('WireBeamInterpolation', WireRestShape='@../topoLines_guide/GuideRestShape', 
                                    radius=0.9, printLog=True, name='InterpolGuide')
    InstrumentCombined.createObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155, 
                                    name='GuideForceField', interpolation='@InterpolGuide')

    InstrumentCombined.createObject('InterventionalRadiologyController', xtip=[1, 0, 0], name='m_ircontroller', 
                                    instruments='InterpolGuide', 
                                    step=0.5, printLog=True, 
                                    listening=True, template='Rigid3d', startingPos=[0, 0, 0, 1, 0, 0, 0], 
                                    rotationInstrument=[0, 0, 0], speed=0, 
                                    controlledInstrument=0)
    InstrumentCombined.createObject('LinearSolverConstraintCorrection', wire_optimization='true', printLog=False)
    InstrumentCombined.createObject('FixedConstraint', indices=0, name='FixedConstraint')
    InstrumentCombined.createObject('RestShapeSpringsForceField', points='@m_ircontroller.indexFirstNode', 
                                    angularStiffness=1e8, stiffness=1e8)
    