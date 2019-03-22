from stlib.scene import Scene

def addEdgeSetTopology(self):
    self.createObject("EdgeSetTopologyContainer", name='container')
    self.createObject("EdgeSetTopologyModifier", name='modifier')
    self.createObject("EdgeSetTopologyAlgorithms", name='topoAlgo') #template='Rigid'/>
    self.createObject("EdgeSetGeometryAlgorithms", name='geomAlgo') #  template='Rigid'/>
    return self

def TopoLines(node, name="TopoLines", **kwargs):
    c = node.createChild(name)
	
    c.createObject("MechanicalObject", template='Rigid3', name='dofs')
    addEdgeSetTopology(c)
    c.createObject("WireRestShape", template='Rigid3', name='wire', **kwargs)
    
    return c

def createScene(root):
    root = Scene(root, plugins=["BeamAdapter"])
    root.dt = 0.05
    root.gravity = [0,0,0]    

    root.createObject("VisualStyle", displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings hideForceFields')
    root.createObject("FreeMotionAnimationLoop")
    root.createObject("LCPConstraintSolver", mu='0.1', tolerance='1e-6', maxIt='1000', build_lcp='false')

    root.createObject("CollisionPipeline", depth='6', verbose='1', draw='0')
    root.createObject("BruteForceDetection", name='N2')
    root.createObject("LocalMinDistance", name='localmindistance', alarmDistance='3', contactDistance='1', angleCone='0.02')
    root.createObject("CollisionResponse", name='Response', response='FrictionContact')
    root.createObject("CollisionGroup", name='Group')
    
    cath = TopoLines(root, name="TopoCath", length='1000.0',  straightLength='600', spireDiameter='4000.0', 
                                                 spireHeight='0.0', densityOfBeams='40 10', numEdges='200', numEdgesCollis='40 20', 
                                                 youngModulus='10000',
                                                 youngModulusExtremity='10000')
    guide = TopoLines(root, name="TopoGuide", length='1000.0',  straightLength='980', spireDiameter='25.0', 
                                                 spireHeight='0.0', densityOfBeams='30 5', numEdges='200', numEdgesCollis='50 10', 
                                                 youngModulus='10000',
                                                 youngModulusExtremity='10000')
    coils = TopoLines(root, name="TopoCoils", length='600.0',  straightLength='540', spireDiameter='7.0', 
                                                 spireHeight='5.0', densityOfBeams='40 20', numEdges='400', numEdgesCollis='30 30', 
                                                 youngModulus='168000',
                                                 youngModulusExtremity='168000')

    i = root.createChild("Instrument")
    i.createObject("EulerImplicit", rayleighStiffness='0.2', rayleighMass='0.1', printLog='false')
    i.createObject("BTDLinearSolver", subpartSolve='0', verification='0', verbose='0')

    i.createObject("RegularGrid", name='meshLinesCombined', 
                   nx='60', ny='1', nz='1',
                   xmin='0.0', xmax='1.0',
                   ymin='0', ymax='0',
                   zmin='1', zmax='1')

    i.createObject("MechanicalObject", template='Rigid3', name='dofs', showIndices='0', ry='-90')


#
#		<WireBeamInterpolation name='InterpolGuide' WireRestShape='@../topoLines_guide/GuideRestShape' radius='0.9' printLog='0'/> 
#		<AdaptiveBeamForceFieldAndMass name='GuideForceField' interpolation='@InterpolGuide' massDensity='0.00000155'/> 	
#
#		<WireBeamInterpolation name='InterpolCoils' WireRestShape='@../topoLines_coils/CoilRestShape' radius='0.1' printLog='0'/> 
#		<AdaptiveBeamForceFieldAndMass name='CoilsForceField' interpolation='@InterpolCoils' massDensity='0.000021'/> <!-- platine E = 168000 MPa // 21000 Kg/m3-->	
#
#		
#		<InterventionalRadiologyController template='Rigid' name='m_ircontroller' printLog='0' xtip='1 0 0'  step='3' rotationInstrument='0 0 0' 
#							controlledInstrument='0' startingPos='0 0 0 0 -0.7071068 0 0.7071068' speed='0' instruments='InterpolCatheter InterpolGuide InterpolCoils' />
#
#		<LinearSolverConstraintCorrection printLog='false' wire_optimization='true'/>
#
#		<FixedConstraint name='FixedConstraint' indices='0' />
#<!--         <FixedConstraint name='FixedConstraint2' indices='@m_ircontroller.indexFirstNode' /> -->
#        <RestShapeSpringsForceField  points='@m_ircontroller.indexFirstNode' stiffness='1e8' angularStiffness='1e8'   />
#
#		<Node name='Collis' activated='true'>
#				<EdgeSetTopologyContainer name='collisEdgeSet' />
#				<EdgeSetTopologyModifier   name='colliseEdgeModifier' />
#				<MechanicalObject name='CollisionDOFs'/>
#				<MultiAdaptiveBeamMapping  name='collisMap'  controller='../m_ircontroller' useCurvAbs='1' printLog='0'/> 
#				<Line proximity='0.0' group='1'/>
#				<Point proximity='0.0' group='1'/>
#		</Node>	
#
#		<Node name='VisuCatheter' activated='true'>
#			<MechanicalObject name='Quads' />
#			<QuadSetTopologyContainer  name='ContainerCath' />
#			<QuadSetTopologyModifier   name='Modifier' />
#			<QuadSetTopologyAlgorithms name='TopoAlgo'  template='Vec3d' />
#			<QuadSetGeometryAlgorithms name='GeomAlgo'  template='Vec3d' />
#			<Edge2QuadTopologicalMapping nbPointsOnEachCircle='10' radius='2' input='@../../topoLines_cath/meshLinesCath' output='@ContainerCath' flipNormals='true'/>
#
#			<AdaptiveBeamMapping  name='VisuMapCath' useCurvAbs='1' printLog='0' interpolation='@../InterpolCatheter' input='@../DOFs' output='@Quads' isMechanical='false'  />
#
#			<Node name='VisuOgl' activated='true'>
#				<OglModel name='Visual' color='0.7 0.7 0.7'  quads='@../ContainerCath.quads' material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20'/>
#				<IdentityMapping input='@../Quads' output='@Visual'/>
#			</Node>	
#		</Node>
#	
#		<Node name='VisuGuide' activated='true'>
#			<MechanicalObject name='Quads' />
#			<QuadSetTopologyContainer  name='ContainerGuide' />
#			<QuadSetTopologyModifier   name='Modifier' />
#			<QuadSetTopologyAlgorithms name='TopoAlgo'  template='Vec3d' />
#			<QuadSetGeometryAlgorithms name='GeomAlgo'  template='Vec3d' />
#			<Edge2QuadTopologicalMapping nbPointsOnEachCircle='10' radius='1' input='@../../topoLines_guide/meshLinesGuide' output='@ContainerGuide'  flipNormals='true' listening='true'/> 
#			<AdaptiveBeamMapping  name='visuMapGuide' useCurvAbs='1' printLog='0' interpolation='@../InterpolGuide' input='@../DOFs' output='@Quads' isMechanical='false' />
#			<Node name='VisuOgl'>
#				<OglModel name='Visual' color='0.2 0.2 0.8'  material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20' quads='@../ContainerGuide.quads'/>
#				<IdentityMapping input='@../Quads' output='@Visual'/>
#			</Node>		
#		</Node> 
#	   
#		<Node name='VisuCoils' activated='true'>
#			<MechanicalObject name='Quads'/>
#			<QuadSetTopologyContainer  name='ContainerCoils' />
#			<QuadSetTopologyModifier   name='Modifier' />
#			<QuadSetTopologyAlgorithms name='TopoAlgo'  template='Vec3d' />
#			<QuadSetGeometryAlgorithms name='GeomAlgo'  template='Vec3d' />
#			<Edge2QuadTopologicalMapping nbPointsOnEachCircle='10' radius='0.3'  input='@../../topoLines_coils/meshLinesCoils' output='@ContainerCoils'  flipNormals='true' listening='true' /> 
#			<AdaptiveBeamMapping  name='visuMapCoils' useCurvAbs='1' printLog='0' interpolation='@../InterpolCoils' input='@../DOFs' output='@Quads' isMechanical='false' />
#			<Node name='VisuOgl'>
#				<OglModel name='Visual' color='0.2 0.8 0.2'  material='texture Ambient 1 0.2 0.2 0.2 0.0 Diffuse 1 1.0 1.0 1.0 1.0 Specular 1 1.0 1.0 1.0 1.0 Emissive 0 0.15 0.05 0.05 0.0 Shininess 1 20' quads='@../ContainerCoils.quads'/>
#				<IdentityMapping input='@../Quads' output='@Visual'/>
#			</Node>		
#		</Node>   
#	</Node>
#
# <!--
#	<Node name='CollidedPoints'>
#		<MechanicalObject name='DOFs1' position='0 0 40' />		
#		<FixedConstraint name='m_fixedConstraint' indices='0' />
#		<Point radius='10'/>
#	</Node> 
#	-->
#	
#
#	<Node name='CollisionModel'> 
#		<MeshObjLoader name='meshLoader' filename='mesh/phantom.obj' triangulate='true' flipNormals='1'/>
#		<Mesh position='@meshLoader.position'  triangles='@meshLoader.triangles'/>
#		<MechanicalObject name='DOFs1' position='0 0 400' scale='3' ry='90'/>	
#			<Triangle simulated='0' moving='0'/>	
#		<Line simulated='0' moving='0'/>
#		<Point simulated='0' moving='0'/>
#		<OglModel name='Visual' fileMesh='mesh/phantom.obj' color='1 0 0 0.3' scale='3' ry='90'/>
#	</Node>
#	
#	
#	
#</Node>

