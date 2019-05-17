import Sofa

from splib.animation import AnimationManager
from beamadapter.parts import createGuide, createInstrumentsCombined, createGeometry

def createScene(rootNode):

    rootNode.createObject('RequiredPlugin', pluginName='SoftRobots')
    rootNode.createObject('RequiredPlugin', pluginName='BeamAdapter')

    AnimationManager(rootNode)  

    rootNode.createObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('LCPConstraintSolver', mu='0.1', tolerance='1e-10', maxIt='1000', build_lcp='false')
    rootNode.createObject('CollisionPipeline', draw='0', depth='6', verbose='1')
    rootNode.createObject('BruteForceDetection', name='N2')
    rootNode.createObject('LocalMinDistance', contactDistance='1', alarmDistance='3', name='localmindistance', angleCone='0.02')
    rootNode.createObject('CollisionResponse', name='Response', response='FrictionContact')
    rootNode.createObject('CollisionGroup', name='Group')

    topoLines_cath = createGuide(rootNode, 'cath', straightLength=600.0, length=1000.0, 
                numEdges=200, youngModulus=10000, spireDiameter=4000.0, 
                numEdgesCollis=[40,20], spireHeight=0.0, densityOfBeams=[40,10], 
                youngModulusExtremity=10000)

    topoLines_guide = createGuide(rootNode, 'guide', straightLength=980.0, length=1000.0, 
                numEdges=200, youngModulus=20000, spireDiameter=25.0, 
                numEdgesCollis=[50,10], spireHeight=0.0, densityOfBeams=[30,5], 
                youngModulusExtremity=20000)

    topoLines_coils = createGuide(rootNode, 'coils', straightLength=540.0, length=600.0, 
                numEdges=400, youngModulus=16800, spireDiameter=7.0, 
                numEdgesCollis=[30,30], spireHeight=5.0, densityOfBeams=[40,20], 
                youngModulusExtremity=1680)

    InstrumentCombined = createInstrumentsCombined(rootNode, xtip=[1, 0, 0], instruments=['cath','guide','coils'], 
                              startingPos=[0, 0, 0, 1, 0, 0, 0], controlledInstrument=0)

    vessels = createGeometry(rootNode, 'mesh/carotids.stl', scale=1.5, rotation=[-40.0, 0.0, 0.0])

    