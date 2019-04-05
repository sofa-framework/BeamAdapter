"""
3instruments_collis1Python
is based on the scene 
/Users/duriez/Desktop/SOFA/BeamAdapter/example/3instruments_collis1.scn
but it uses the SofaPython plugin. 
Further informations on the usage of the plugin can be found in 
sofa/applications/plugins/SofaPython/doc/SofaPython.pdf
To lance the scene, type 
runSofa /Users/duriez/Desktop/SOFA/BeamAdapter/example/3instruments_collis1Python.scn

The current file has been written by the python script
scn2python.py
Author of scn2python.py: Christoph PAULUS, christoph.paulus@inria.fr
"""

import Sofa

class 3instruments_collis1 (Sofa.PythonScriptController):

    def createGraph(self,rootNode):

        # rootNode
        rootNode.createObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings hideForceFields')
        rootNode.createObject('FreeMotionAnimationLoop')
        rootNode.createObject('LCPConstraintSolver', mu='0.1', tolerance='1e-10', maxIt='1000', build_lcp='false')
        rootNode.createObject('CollisionPipeline', draw='0', depth='6', verbose='1')
        rootNode.createObject('BruteForceDetection', name='N2')
        rootNode.createObject('LocalMinDistance', contactDistance='1', alarmDistance='3', name='localmindistance', angleCone='0.02')
        rootNode.createObject('CollisionResponse', name='Response', response='FrictionContact')
        rootNode.createObject('CollisionGroup', name='Group')

        # rootNode/topoLines_cath
        topoLines_cath = rootNode.createChild('topoLines_cath')
        topoLines_cath.createObject('WireRestShape', name='catheterRestShape', straightLength='600', length='1000.0', numEdges='200', youngModulus='10000', spireDiameter='4000.0', numEdgesCollis='40 20', printLog='1', template='Rigid', spireHeight='0.0', densityOfBeams='40 10', youngModulusExtremity='10000')
        topoLines_cath.createObject('EdgeSetTopologyContainer', name='meshLinesCath')
        topoLines_cath.createObject('EdgeSetTopologyModifier', name='Modifier')
        topoLines_cath.createObject('EdgeSetTopologyAlgorithms', name='TopoAlgo', template='Rigid')
        topoLines_cath.createObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid')
        topoLines_cath.createObject('MechanicalObject', name='dofTopo1', template='Rigid')

        # rootNode/topoLines_guide
        topoLines_guide = rootNode.createChild('topoLines_guide')
        topoLines_guide.createObject('WireRestShape', name='GuideRestShape', straightLength='980.0', length='1000.0', numEdges='200', youngModulus='10000', spireDiameter='25', numEdgesCollis='50 10', printLog='1', template='Rigid', spireHeight='0.0', densityOfBeams='30 5', youngModulusExtremity='10000')
        topoLines_guide.createObject('EdgeSetTopologyContainer', name='meshLinesGuide')
        topoLines_guide.createObject('EdgeSetTopologyModifier', name='Modifier')
        topoLines_guide.createObject('EdgeSetTopologyAlgorithms', name='TopoAlgo', template='Rigid')
        topoLines_guide.createObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid')
        topoLines_guide.createObject('MechanicalObject', name='dofTopo2', template='Rigid')

        # rootNode/topoLines_coils
        topoLines_coils = rootNode.createChild('topoLines_coils')
        topoLines_coils.createObject('WireRestShape', name='CoilRestShape', straightLength='540.0', length='600.0', numEdges='400', youngModulus='168000', spireDiameter='7', numEdgesCollis='30 30', printLog='1', template='Rigid', spireHeight='5.0', densityOfBeams='40 20', youngModulusExtremity='168000')
        topoLines_coils.createObject('EdgeSetTopologyContainer', name='meshLinesCoils')
        topoLines_coils.createObject('EdgeSetTopologyModifier', name='Modifier')
        topoLines_coils.createObject('EdgeSetTopologyAlgorithms', name='TopoAlgo', template='Rigid')
        topoLines_coils.createObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid')
        topoLines_coils.createObject('MechanicalObject', name='dofTopo3', template='Rigid')

        # rootNode/InstrumentCombined
        InstrumentCombined = rootNode.createChild('InstrumentCombined')
        InstrumentCombined.createObject('EulerImplicit', rayleighStiffness='0.2', printLog='false', rayleighMass='0.1')
        InstrumentCombined.createObject('BTDLinearSolver', verification='0', subpartSolve='0', verbose='0')
        InstrumentCombined.createObject('RegularGrid', zmax='1', name='meshLinesCombined', zmin='1', nx='60', ny='1', nz='1', xmax='1.0', xmin='0.0', ymin='0', ymax='0')
        InstrumentCombined.createObject('MechanicalObject', showIndices='0', name='DOFs', template='Rigid', ry='-90')
        InstrumentCombined.createObject('WireBeamInterpolation', WireRestShape='@../topoLines_cath/catheterRestShape', radius='1', printLog='1', name='InterpolCatheter')
        InstrumentCombined.createObject('AdaptiveBeamForceFieldAndMass', massDensity='0.00000155', name='CatheterForceField', interpolation='@InterpolCatheter')
        InstrumentCombined.createObject('WireBeamInterpolation', WireRestShape='@../topoLines_guide/GuideRestShape', radius='0.9', printLog='1', name='InterpolGuide')
        InstrumentCombined.createObject('AdaptiveBeamForceFieldAndMass', massDensity='0.00000155', name='GuideForceField', interpolation='@InterpolGuide')
        InstrumentCombined.createObject('WireBeamInterpolation', WireRestShape='@../topoLines_coils/CoilRestShape', radius='0.1', printLog='1', name='InterpolCoils')
        InstrumentCombined.createObject('AdaptiveBeamForceFieldAndMass', massDensity='0.000021', name='CoilsForceField', interpolation='@InterpolCoils')
        InstrumentCombined.createObject('InterventionalRadiologyController', xtip='1 0 0', name='m_ircontroller', instruments='InterpolCatheter InterpolGuide InterpolCoils', step='3', printLog='1', template='Rigid', startingPos='0 0 0 0 -0.7071068 0 0.7071068', rotationInstrument='0 0 0', speed='0', controlledInstrument='0')
        InstrumentCombined.createObject('LinearSolverConstraintCorrection', wire_optimization='true', printLog='false')
        InstrumentCombined.createObject('FixedConstraint', indices='0', name='FixedConstraint')
        InstrumentCombined.createObject('RestShapeSpringsForceField', points='@m_ircontroller.indexFirstNode', angularStiffness='1e8', stiffness='1e8')

        # rootNode/InstrumentCombined/Collis
        Collis = InstrumentCombined.createChild('Collis')
        Collis.activated = 'true'
        Collis.createObject('EdgeSetTopologyContainer', name='collisEdgeSet')
        Collis.createObject('EdgeSetTopologyModifier', name='colliseEdgeModifier')
        Collis.createObject('MechanicalObject', name='CollisionDOFs')
        Collis.createObject('MultiAdaptiveBeamMapping', controller='../m_ircontroller', useCurvAbs='1', printLog='1', name='collisMap')
        Collis.createObject('Line', proximity='0.0', group='1')
        Collis.createObject('Point', proximity='0.0', group='1')

        # rootNode/InstrumentCombined/VisuCatheter
        VisuCatheter = InstrumentCombined.createChild('VisuCatheter')
        VisuCatheter.activated = 'true'
        VisuCatheter.createObject('MechanicalObject', name='Quads')
        VisuCatheter.createObject('QuadSetTopologyContainer', name='ContainerCath')
        VisuCatheter.createObject('QuadSetTopologyModifier', name='Modifier')
        VisuCatheter.createObject('QuadSetTopologyAlgorithms', name='TopoAlgo', template='Vec3d')
        VisuCatheter.createObject('QuadSetGeometryAlgorithms', name='GeomAlgo', template='Vec3d')
        VisuCatheter.createObject('Edge2QuadTopologicalMapping', nbPointsOnEachCircle='10', radius='2', flipNormals='true', input='@../../topoLines_cath/meshLinesCath', output='@ContainerCath')
        VisuCatheter.createObject('AdaptiveBeamMapping', interpolation='@../InterpolCatheter', name='VisuMapCath', output='@Quads', isMechanical='false', input='@../DOFs', useCurvAbs='1', printLog='1')

        # rootNode/InstrumentCombined/VisuCatheter/VisuOgl
        VisuOgl = VisuCatheter.createChild('VisuOgl')
        VisuOgl.activated = 'true'
        VisuOgl.createObject('OglModel', color='0.7 0.7 0.7', quads='@../ContainerCath.quads', material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20', name='Visual')
        VisuOgl.createObject('IdentityMapping', input='@../Quads', output='@Visual')

        # rootNode/InstrumentCombined/VisuGuide
        VisuGuide = InstrumentCombined.createChild('VisuGuide')
        VisuGuide.activated = 'true'
        VisuGuide.createObject('MechanicalObject', name='Quads')
        VisuGuide.createObject('QuadSetTopologyContainer', name='ContainerGuide')
        VisuGuide.createObject('QuadSetTopologyModifier', name='Modifier')
        VisuGuide.createObject('QuadSetTopologyAlgorithms', name='TopoAlgo', template='Vec3d')
        VisuGuide.createObject('QuadSetGeometryAlgorithms', name='GeomAlgo', template='Vec3d')
        VisuGuide.createObject('Edge2QuadTopologicalMapping', radius='1', listening='true', input='@../../topoLines_guide/meshLinesGuide', nbPointsOnEachCircle='10', flipNormals='true', output='@ContainerGuide')
        VisuGuide.createObject('AdaptiveBeamMapping', interpolation='@../InterpolGuide', name='visuMapGuide', output='@Quads', isMechanical='false', input='@../DOFs', useCurvAbs='1', printLog='1')

        # rootNode/InstrumentCombined/VisuGuide/VisuOgl
        VisuOgl = VisuGuide.createChild('VisuOgl')
        VisuOgl.createObject('OglModel', color='0.2 0.2 0.8', quads='@../ContainerGuide.quads', material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20', name='Visual')
        VisuOgl.createObject('IdentityMapping', input='@../Quads', output='@Visual')

        # rootNode/InstrumentCombined/VisuCoils
        VisuCoils = InstrumentCombined.createChild('VisuCoils')
        VisuCoils.activated = 'true'
        VisuCoils.createObject('MechanicalObject', name='Quads')
        VisuCoils.createObject('QuadSetTopologyContainer', name='ContainerCoils')
        VisuCoils.createObject('QuadSetTopologyModifier', name='Modifier')
        VisuCoils.createObject('QuadSetTopologyAlgorithms', name='TopoAlgo', template='Vec3d')
        VisuCoils.createObject('QuadSetGeometryAlgorithms', name='GeomAlgo', template='Vec3d')
        VisuCoils.createObject('Edge2QuadTopologicalMapping', radius='0.3', listening='true', input='@../../topoLines_coils/meshLinesCoils', nbPointsOnEachCircle='10', flipNormals='true', output='@ContainerCoils')
        VisuCoils.createObject('AdaptiveBeamMapping', interpolation='@../InterpolCoils', name='visuMapCoils', output='@Quads', isMechanical='false', input='@../DOFs', useCurvAbs='1', printLog='1')

        # rootNode/InstrumentCombined/VisuCoils/VisuOgl
        VisuOgl = VisuCoils.createChild('VisuOgl')
        VisuOgl.createObject('OglModel', color='0.2 0.8 0.2', quads='@../ContainerCoils.quads', material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20', name='Visual')
        VisuOgl.createObject('IdentityMapping', input='@../Quads', output='@Visual')

        # rootNode/CollisionModel
        CollisionModel = rootNode.createChild('CollisionModel')
        CollisionModel.createObject('MeshObjLoader', filename='mesh/phantom.obj', flipNormals='1', triangulate='true', name='meshLoader')
        CollisionModel.createObject('Mesh', position='@meshLoader.position', triangles='@meshLoader.triangles')
        CollisionModel.createObject('MechanicalObject', position='0 0 400', scale='3', name='DOFs1', ry='90')
        CollisionModel.createObject('Triangle', moving='0', simulated='0')
        CollisionModel.createObject('Line', moving='0', simulated='0')
        CollisionModel.createObject('Point', moving='0', simulated='0')
        CollisionModel.createObject('OglModel', color='1 0 0 0.3', scale='3', fileMesh='mesh/phantom.obj', ry='90', name='Visual')

        return 0;

    def onKeyPressed(self, c):
        return 0;

    def onKeyReleased(self, c):
        return 0;

    def onLoaded(self, node):
        return 0;

    def onMouseButtonLeft(self, mouseX,mouseY,isPressed):
        return 0;

    def onMouseButtonRight(self, mouseX,mouseY,isPressed):
        return 0;

    def onMouseButtonMiddle(self, mouseX,mouseY,isPressed):
        return 0;

    def onMouseWheel(self, mouseX,mouseY,wheelDelta):
        return 0;

    def onGUIEvent(self, strControlID,valueName,strValue):
        return 0;

    def onBeginAnimationStep(self, deltaTime):
        return 0;

    def onEndAnimationStep(self, deltaTime):
        return 0;

    def onScriptEvent(self, senderNode, eventName,data):
        return 0;

    def initGraph(self, node):
        return 0;

    def bwdInitGraph(self, node):
        return 0;

    def storeResetState(self):
        return 0;

    def reset(self):
        return 0;

    def cleanup(self):
        return 0;
