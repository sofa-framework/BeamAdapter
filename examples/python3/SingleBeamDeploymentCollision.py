import sys
import Sofa

def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', pluginName='BeamAdapter SofaMeshCollision SofaBoundaryCondition SofaConstraint SofaMiscCollision SofaDeformable SofaGeneralLinearSolver SofaImplicitOdeSolver')
    rootNode.addObject('VisualStyle', displayFlags='hideVisualModels hideBehaviorModels showCollisionModels hideMappings showInteractionForceFields')

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('DefaultVisualManagerLoop')

    rootNode.addObject('LCPConstraintSolver', mu='0.1', tolerance='1e-10', maxIt='1000', build_lcp='false')
    rootNode.addObject('DefaultPipeline', draw='0', depth='6', verbose='1')
    rootNode.addObject('BruteForceBroadPhase', name='N2')
    rootNode.addObject('BVHNarrowPhase')
    rootNode.addObject('LocalMinDistance', contactDistance='1', alarmDistance='5', name='localmindistance', angleCone='0.02')
    rootNode.addObject('DefaultContactManager', name='Response', response='FrictionContact')

    topoLines = rootNode.addChild('EdgeTopology')
    topoLines.addObject('WireRestShape', name='BeamRestShape', 
                                 straightLength=980.0, length=1000.0, 
                                 numEdges=200, youngModulus=20000, 
                                 spireDiameter=25, numEdgesCollis=[50,10], 
                                 printLog=True, template='Rigid3d', spireHeight=0.0, 
                                 densityOfBeams=[30,5], youngModulusExtremity=20000)
    topoLines.addObject('EdgeSetTopologyContainer', name='meshLines')
    topoLines.addObject('EdgeSetTopologyModifier', name='Modifier')
    topoLines.addObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid3d')
    topoLines.addObject('MechanicalObject', name='dofTopo2', template='Rigid3d')


    BeamMechanics = rootNode.addChild('BeamModel')
    BeamMechanics.addObject('EulerImplicitSolver', rayleighStiffness=0.2, rayleighMass=0.1)
    BeamMechanics.addObject('BTDLinearSolver', verification=False, subpartSolve=False, verbose=False)
    BeamMechanics.addObject('RegularGridTopology', name='MeshLines', 
                                    nx=60, ny=1, nz=1,
                                    xmax=0.0, xmin=0.0, ymin=0, ymax=0, zmax=0, zmin=0,
                                    p0=[0,0,0])
    BeamMechanics.addObject('MechanicalObject', showIndices=False, name='DOFs', template='Rigid3d', ry=-90)
    BeamMechanics.addObject('WireBeamInterpolation', name='BeamInterpolation', WireRestShape='@../EdgeTopology/BeamRestShape', 
                                    radius=0.9, printLog=False)
    BeamMechanics.addObject('AdaptiveBeamForceFieldAndMass', name='BeamForceField', massDensity=0.00000155, interpolation='@BeamInterpolation')
    BeamMechanics.addObject('InterventionalRadiologyController', name='DeployController', template='Rigid3d', instruments='BeamInterpolation', 
                                    startingPos=[0, 0, 0, 0, 0, 0, 1], xtip=[0, 0, 0], printLog=True, 
                                    rotationInstrument=[0, 0, 0], step=5., speed=5., 
                                    listening=True, controlledInstrument=0)
    BeamMechanics.addObject('LinearSolverConstraintCorrection', wire_optimization='true', printLog=False)
    BeamMechanics.addObject('FixedConstraint', indices=0, name='FixedConstraint')
    BeamMechanics.addObject('RestShapeSpringsForceField', points='@DeployController.indexFirstNode', angularStiffness=1e8, stiffness=1e8)
    

    BeamCollis = BeamMechanics.addChild('CollisionModel')
    BeamCollis.activated = True
    BeamCollis.addObject('EdgeSetTopologyContainer', name='collisEdgeSet')
    BeamCollis.addObject('EdgeSetTopologyModifier', name='colliseEdgeModifier')
    BeamCollis.addObject('MechanicalObject', name='CollisionDOFs')
    BeamCollis.addObject('MultiAdaptiveBeamMapping', controller='../DeployController', useCurvAbs=True, printLog=False, name='collisMap')
    BeamCollis.addObject('LineCollisionModel', proximity=0.0)
    BeamCollis.addObject('PointCollisionModel', proximity=0.0)

    Carotids = rootNode.addChild('Carotids')
    Carotids.addObject('MeshSTLLoader', filename='../mesh/carotids.stl', flipNormals=False, triangulate=True, name='meshLoader', rotation=[10.0, 0.0, -90.0])
    Carotids.addObject('MeshTopology', position='@meshLoader.position', triangles='@meshLoader.triangles')
    Carotids.addObject('MechanicalObject', position=[0,0,400], scale=3, name='DOFs1', ry=90)
    Carotids.addObject('TriangleCollisionModel', moving=False, simulated=False)
    Carotids.addObject('LineCollisionModel', moving=False, simulated=False)


def main():
    import SofaRuntime
    import Sofa.Gui

    root = Sofa.Core.Node('root')
    createScene(root)
    Sofa.Simulation.init(root)

    Sofa.Gui.GUIManager.Init('myscene', 'qglviewer')
    Sofa.Gui.GUIManager.createGUI(root, __file__)
    Sofa.Gui.GUIManager.SetDimension(1080, 1080)
    Sofa.Gui.GUIManager.MainLoop(root)
    Sofa.Gui.GUIManager.closeGUI()


# Function used only if this script is called from a python environment
if __name__ == '__main__':
    main()
