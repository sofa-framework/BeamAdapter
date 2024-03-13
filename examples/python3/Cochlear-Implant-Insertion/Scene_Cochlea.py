# -*- coding: utf-8 -*-
# units: mm, kg, s
import Sofa

import os
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'

def createScene(rootNode):
        rootNode.addObject('RequiredPlugin',name='BeamAdapter', pluginName='BeamAdapter ')
        rootNode.addObject('RequiredPlugin',name='SofaPython3', pluginName='SofaPython3')
        rootNode.addObject('RequiredPlugin',name='SOFA Modules', pluginName='Sofa.Component.AnimationLoop Sofa.Component.Collision.Detection.Algorithm Sofa.Component.Collision.Detection.Intersection Sofa.Component.Collision.Geometry Sofa.Component.Collision.Response.Contact Sofa.Component.Constraint.Lagrangian.Correction Sofa.Component.Constraint.Lagrangian.Solver Sofa.Component.Constraint.Projective Sofa.Component.IO.Mesh Sofa.Component.LinearSolver.Direct Sofa.Component.ODESolver.Backward Sofa.Component.SolidMechanics.Spring Sofa.Component.Topology.Container.Constant Sofa.Component.Topology.Container.Dynamic Sofa.Component.Topology.Container.Grid Sofa.Component.Mapping.Linear Sofa.Component.Topology.Mapping Sofa.Component.StateContainer Sofa.Component.Visual Sofa.GL.Component.Rendering3D')
        rootNode.findData('dt').value=0.01
        rootNode.findData('gravity').value=[0,0,0]
        rootNode.addObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideMappings hideForceFields showInteractionForceFields hideWireframe')

        rootNode.addObject('DefaultVisualManagerLoop')
        rootNode.addObject('FreeMotionAnimationLoop')
        rootNode.addObject('GenericConstraintSolver', tolerance="1e-3", maxIterations="5000", resolutionMethod="UnbuiltGaussSeidel")
        rootNode.addObject('CollisionPipeline', depth="6", verbose="0", draw="1")
        rootNode.addObject('BruteForceBroadPhase', name='N2')
        rootNode.addObject('BVHNarrowPhase')
        rootNode.addObject('CollisionResponse', response="FrictionContactConstraint", responseParams="mu=0.65")
        rootNode.addObject('LocalMinDistance', name="Proximity", alarmDistance="0.6", contactDistance="0.44", angleCone="0.01")


        cochleaNode = rootNode.addChild('cochleaNode')
        cochleaNode.addObject('MeshOBJLoader', name='loader', filename='Mesh/cochleeCompleteTroueeSimpleOrientationMegaTrou.obj', flipNormals="false")
        cochleaNode.addObject('MeshTopology',src = '@loader')
        cochleaNode.addObject('MechanicalObject', name='dofs', template='Vec3d', showIndices='false', showIndicesScale='4e-5', rx='0',printLog="0")
        cochleaNode.addObject('TriangleCollisionModel', group='1')
        cochleaNode.addObject('LineCollisionModel', group='1')
        cochleaNode.addObject('PointCollisionModel', group='1')
        visuCochleaNode = cochleaNode.addChild('visuCochleaNode')
        visuCochleaNode.addObject('OglModel', name="VisualModel", color="3.0 0.5 0.0 0.9")


        topoLines_cath = rootNode.addChild('topoLines_cath')
        topoLines_cath.addObject('RodStraightSection', name='StraightSection', 
                                 length=20.0, radius="@../Proximity.contactDistance", 
                                 nbEdgesCollis=20, nbEdgesVisu=20, 
                                 youngModulus=2.5e5, massDensity=0.000005, poissonRatio=0.3)

        topoLines_cath.addObject('WireRestShape', name='catheterRestShape', template="Rigid3d", wireMaterials="@StraightSection")
        topoLines_cath.addObject('EdgeSetTopologyContainer', name="meshLinesCath")
        topoLines_cath.addObject('EdgeSetTopologyModifier', name="Modifier")
        topoLines_cath.addObject('EdgeSetGeometryAlgorithms', name="GeomAlgo", template="Rigid3d")
        topoLines_cath.addObject('MechanicalObject', template="Rigid3d", name="dofTopo1")

        RefStartingPos = rootNode.addChild('RefStartingPos')
        RefStartingPos.addObject('MechanicalObject', name="ReferencePos", template="Rigid3d", position="-3 1.5 0.3  1 0 0 0")

        InstrumentCombined = rootNode.addChild('InstrumentCombined')
        InstrumentCombined.addObject('EulerImplicitSolver', rayleighStiffness="0.01", rayleighMass="0.03", printLog=False )
        InstrumentCombined.addObject('BTDLinearSolver')
        InstrumentCombined.addObject('RegularGridTopology', name="meshLinesCombined", nx="100", ny="1", nz="1")

        InstrumentCombined.addObject('MechanicalObject', template="Rigid3d", name="DOFs" )
        InstrumentCombined.addObject('InterventionalRadiologyController', template="Rigid3d", name="m_ircontroller", printLog=False, xtip="0.1",speed ='4',   step="0.1", rotationInstrument="0", controlledInstrument="0", startingPos="@../RefStartingPos/ReferencePos.position", instruments="InterpolCatheter")
        InstrumentCombined.addObject('WireBeamInterpolation', name="InterpolCatheter", WireRestShape="@../topoLines_cath/catheterRestShape", radius="3.0", printLog=False)
        InstrumentCombined.addObject('AdaptiveBeamForceFieldAndMass', name="CatheterForceField", massDensity="0.000005", interpolation="@InterpolCatheter", printLog=False)
        InstrumentCombined.addObject('LinearSolverConstraintCorrection', printLog=False, wire_optimization="true")
        InstrumentCombined.addObject("FixedProjectiveConstraint", indices="0")
        InstrumentCombined.addObject('RestShapeSpringsForceField', name="MeasurementFF", points="@m_ircontroller.indexFirstNode",  stiffness="1e10", recompute_indices="1", angularStiffness="1e10", external_rest_shape="@../RefStartingPos/ReferencePos", external_points="0", drawSpring="1", springColor="1 0 0 1")

        CollisInstrumentCombined = InstrumentCombined.addChild('CollisInstrumentCombined')
        CollisInstrumentCombined.addObject('EdgeSetTopologyContainer', name="collisEdgeSet")
        CollisInstrumentCombined.addObject('EdgeSetTopologyModifier', name="colliseEdgeModifier")
        CollisInstrumentCombined.addObject('MechanicalObject', name="CollisionDOFs")
        CollisInstrumentCombined.addObject('MultiAdaptiveBeamMapping', name="multimapp", ircontroller="../m_ircontroller", useCurvAbs="1", printLog="false")
        CollisInstrumentCombined.addObject('LineCollisionModel' )
        CollisInstrumentCombined.addObject('PointCollisionModel')

        visuInstrumentCombined = InstrumentCombined.addChild('visuInstrumentCombined')
        visuInstrumentCombined.addObject('MechanicalObject', name="Quads")
        visuInstrumentCombined.addObject('QuadSetTopologyContainer', name="ContainerCath")
        visuInstrumentCombined.addObject('QuadSetTopologyModifier', name="Modifier" )
        visuInstrumentCombined.addObject('QuadSetGeometryAlgorithms', name="GeomAlgo", template="Vec3d")
        visuInstrumentCombined.addObject('Edge2QuadTopologicalMapping', nbPointsOnEachCircle="10", radius="@../../Proximity.contactDistance", input="@../../topoLines_cath/meshLinesCath", output="@ContainerCath", flipNormals="true",printLog=False)
        visuInstrumentCombined.addObject('AdaptiveBeamMapping', name="VisuMapCath", useCurvAbs="1", printLog=False, isMechanical="false",  interpolation="@../InterpolCatheter")


        realVisuInstrumentCombined = visuInstrumentCombined.addChild('realVisuInstrumentCombined')
        realVisuInstrumentCombined.addObject('OglModel',name="VisualCathOGL", src="@../ContainerCath", color='white')
        realVisuInstrumentCombined.addObject('IdentityMapping', input="@../Quads", output="@VisualCathOGL")
        return rootNode
