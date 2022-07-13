# -*- coding: utf-8 -*-
# units: mm, kg, s
import Sofa

import os
path = os.path.dirname(os.path.abspath(__file__))+'/mesh/'

def createScene(rootNode):
        rootNode.addObject('RequiredPlugin',name='BeamAdapter', pluginName='BeamAdapter ')
        rootNode.addObject('RequiredPlugin',name='SofaPython3', pluginName='SofaPython3')
        rootNode.addObject('RequiredPlugin',name='SOFA Modules', pluginName='Sofa.Component.AnimationLoop Sofa.Component.Collision.Detection.Algorithm Sofa.Component.Collision.Detection.Intersection Sofa.Component.Collision.Geometry Sofa.Component.Collision.Response.Contact Sofa.Component.Constraint.Lagrangian.Correction Sofa.Component.Constraint.Lagrangian.Solver Sofa.Component.Constraint.Projective Sofa.Component.IO.Mesh Sofa.Component.LinearSolver.Direct Sofa.Component.ODESolver.Backward Sofa.Component.SolidMechanics.Spring Sofa.Component.Topology.Container.Constant Sofa.Component.Topology.Container.Dynamic Sofa.Component.Topology.Container.Grid  Sofa.Component.Topology.Mapping Sofa.Component.Visual Sofa.GL.Component.Rendering3D')
        rootNode.findData('dt').value=0.01
        rootNode.findData('gravity').value=[0,0,0]
        rootNode.addObject('VisualStyle', displayFlags='showVisualModels hideBehaviorModels showCollisionModels hideMappings hideForceFields showInteractionForceFields hideWireframe')
        rootNode.addObject('DefaultVisualManagerLoop')
        rootNode.addObject('FreeMotionAnimationLoop')
        rootNode.addObject('GenericConstraintSolver', tolerance="1e-3", maxIterations="5000", unbuilt="false")
        rootNode.addObject('DefaultPipeline', depth="6", verbose="0", draw="1")
        rootNode.addObject('BruteForceBroadPhase', name='N2')
        rootNode.addObject('BVHNarrowPhase')
        rootNode.addObject('DefaultContactManager', response="FrictionContact", responseParams="mu=0.65")
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

        #membraneNode = rootNode.addChild('membraneNode')
        #membraneNode.addObject('EulerImplicitSolver', rayleighStiffness='0.0', rayleighMass='0.0')
        #membraneNode.addObject('SparseLDLSolver')
        #membraneNode.addObject('MeshObjLoader', name='loader', filename='Mesh/membraneBasilaireBetterFit.obj', flipNormals="false")
        #membraneNode.addObject('Vertex2Frame', name="frames", template="Rigid3d", position="@loader.position", normals="@loader.normals", invertNormals="false")
        ##membraneNode.addObject('Mesh',src ='@loader')
        #membraneNode.addObject('MechanicalObject', name='MO', template="Rigid3d", position="@frames.frames", showIndices="false", showIndicesScale="0.000005")
        #membraneNode.addObject('TriangleSetTopologyContainer', name="coarseTopo", src="@loader")

        #membraneNode.addObject('UniformMass', showAxisSizeFactor="0.1", totalMass="0.001")
        #membraneNode.addObject('TriangularShellForceField', name="FEM", youngModulus="1e4", poissonRatio="0.33", rayleighStiffness="0", thickness="4.0e-1", measure="Strain (norm)")
        #membraneNode.addObject('RestShapeSpringsForceField', points="7 8 10 11 12 13 14 19 20 22 23 35 36 51 52 53 56 57 58 60 61 62 64 67 80 81 82 83 89 90 91 99 100 101 114 115 116 121 122 124 126 130 131 134 152 153 159 160 164 172 173 188 189 197 198 199 201 212 213 214 215 216 217 218 220 224 225 226 234 235 236 239 243 244 245 246 248 250 254 258 259 267 271 277 284 285 287 294 296 297 298 299 300 301 302 303 310 311 317 329 330 331 332 334 336 338 339 340 341 343 347 348 350 357 359 361 364 366 367 369 370 371 372 374 376 381 383 385 386 390 396 399 403 404 405 406 407 433 435 450 453 454 456 460 464 465 470 471 472 474 478 479 480 488 491 505 506 510 512 513 514 516 520 523 525 526 529 534 535 536 554 557 558 561 562 563 564 565 567 568 569 577 586 590 592 593 598 599 600 605 610 612 616 627 630 631 632 639 647 649 650 651 652 653 654 655 662 663 668 675 687 692 693 694 695 696 697 698", stiffness="1000000", angularStiffness="100000")
        #membraneNode.addObject('LinearSolverConstraintCorrection')

        #collMembraneNode = membraneNode.addChild('collMembraneNode')
        #collMembraneNode.addObject('MechanicalObject', name="Coll_MO", template="Vec3d", src="@../loader")
        #collMembraneNode.addObject('Triangle', bothSide="1", group='1')
        #collMembraneNode.addObject('LineCollisionModel', bothSide="1", group='1')
        #collMembraneNode.addObject('PointCollisionModel', bothSide="1", group='1')
        #collMembraneNode.addObject('IdentityMapping', input="@../MO", output="@Coll_MO")

        #visuMembraneNode = membraneNode.addChild('visuMembraneNode')
        #visuMembraneNode.addObject('OglModel', name="Visual", color="red")

        topoLines_cath = rootNode.addChild('topoLines_cath')
        topoLines_cath.addObject('WireRestShape', template="Rigid3d", printLog=False, name="catheterRestShape", length="20", straightLength="20", spireDiameter="0", spireHeight="0.0", densityOfBeams="40", numEdges="20", numEdgesCollis="20", youngModulus="2.5e5", youngModulusExtremity="2.5e5", radius="@../Proximity.contactDistance")
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
        InstrumentCombined.addObject("FixedConstraint", indices="0")
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
