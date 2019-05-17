import Sofa

from splib.animation import AnimationManager
from beamadapter.parts import createGuide, createInstrumentsCombined

def createScene(rootNode):

    rootNode.createObject('RequiredPlugin', pluginName='SoftRobots')
    rootNode.createObject('RequiredPlugin', pluginName='BeamAdapter')

    AnimationManager(rootNode)  

    rootNode.createObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')
    rootNode.createObject('FreeMotionAnimationLoop')
    rootNode.createObject('LCPConstraintSolver', mu='0.1', tolerance='1e-10', maxIt='1000', build_lcp='false')
    rootNode.createObject('BruteForceDetection', name='N2')
    rootNode.createObject('LocalMinDistance', contactDistance='1', alarmDistance='3', name='localmindistance', angleCone='0.02')


    topoLines_guide = createGuide(rootNode, 'guide', straightLength=980.0, length=1000.0, 
                numEdges=200, spireDiameter=25, 
                numEdgesCollis=[50,10], spireHeight=0.0, densityOfBeams=[30,5], youngModulusExtremity=10000)
    instrumentsCombined = createInstrumentsCombined(rootNode)
