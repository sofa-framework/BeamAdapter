import Sofa

from beamadapter.parts import createGuide, createInstrumentsCombined, createGeometry

def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', pluginName='BeamAdapter SofaMiscCollision SofaConstraint SofaImplicitOdeSolver SofaGeneralLinearSolver SofaBoundaryCondition SofaDeformable SofaTopologyMapping SofaOpenglVisual SofaMeshCollision' )

    rootNode.addObject('DefaultVisualManagerLoop')
    rootNode.addObject('FreeMotionAnimationLoop')

    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields')

    rootNode.addObject('LCPConstraintSolver', mu='0.1', tolerance='1e-10', maxIt='1000', build_lcp='false')
    rootNode.addObject('DefaultPipeline', draw='0', depth='6', verbose='1')
    rootNode.addObject('BruteForceBroadPhase', name='N2')
    rootNode.addObject('BVHNarrowPhase')
    rootNode.addObject('LocalMinDistance', contactDistance='1', alarmDistance='3', name='localmindistance', angleCone='0.02')
    rootNode.addObject('DefaultContactManager', name='Response', response='FrictionContactConstraint')


    topoLines_guide = createGuide(rootNode, 'guide', straightLength=980.0, length=1000.0, 
                numEdges=200, spireDiameter=25, 
                numEdgesCollis=[50,10], spireHeight=0.0, densityOfBeams=[30,5], youngModulusExtremity=10000)
    instrumentsCombined = createInstrumentsCombined(rootNode)
    
    vessels = createGeometry(rootNode, '../mesh/carotids.stl', scale=1.5, rotation=[-40.0, 0.0, 0.0])
