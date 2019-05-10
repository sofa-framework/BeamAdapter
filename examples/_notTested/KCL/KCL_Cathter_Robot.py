import Sofa

import os
path = 'mesh/'

def createScene(rootNode):

  rootNode.createObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels hideCollisionModels hideBoundingCollisionModels hideForceFields hideInteractionForceFields hideWireframe')
  rootNode.createObject('RequiredPlugin', name="Robot Model", pluginName='BeamAdapter')
  rootNode.createObject('FreeMotionAnimationLoop')
  rootNode.createObject('GenericConstraintSolver', tolerance="0.001", maxIterations="1000")
  rootNode.createObject('CollisionPipeline', depth='6')
  rootNode.createObject('BruteForceDetection')
  rootNode.createObject('LocalMinDistance', alarmDistance='0.1', contactDistance='0.01', angleCone='0.1')
  rootNode.createObject('CollisionResponse', response='Contact')

  heart = rootNode.createChild('Heart Model')
  heart.createObject('Mesh', name='mesh', filename=path+'Heart_Phantom_Redux_Flip.obj')
  heart.createObject('MechanicalObject', name='dofs', src='@mesh')
  heart.createObject('TriangleSetTopologyContainer', src='@mesh')
  heart.createObject('TriangleSetTopologyModifier')
  heart.createObject('TriangleSetTopologyAlgorithms')
  heart.createObject('TriangleSetGeometryAlgorithms')
  heart.createObject('VisualModel', src='@mesh')

  catheter = rootNode.createChild('CatheterModel')
  catheter.createObject('WireRestShape', template='Rigid', name='catheterRestShape', length='20', straightLength='99', spireDiameter='20', spireHeight='0.0', densityOfBeams='49 1', numEdges='50',numEdgesCollis='49 1', youngModulus='1e7', youngModulusExtremity='10')
  catheter.createObject('EdgeSetTopologyContainer')
  catheter.createObject('EdgeSetTopologyModifier')
  catheter.createObject('EdgeSetTopologyAlgorithms', template='Rigid')
  catheter.createObject('EdgeSetGeometryAlgorithms', template='Rigid')
  catheter.createObject('MechanicalObject', template='Rigid', name='dofsTopo')

  instruments = rootNode.createChild('InstrumentsCombined')
  instruments.createObject('EulerImplicit', rayleighStiffness='0.01', rayleighMass='0.03')
  instruments.createObject('BTDLinearSolver')
  instruments.createObject('RegularGrid', nx='100', ny='1', nz='1', xmin='0.0', xmax='1.0', ymin='0', ymax='0', zmin='1', zmax='1')
  instruments.createObject('MechanicalObject', template='Rigid', name='dofs')
  instruments.createObject('InterventionalRadiologyController', template='Rigid', xtip='0', step='0.02', rotationInstrument='0 0', controlledInstrument='0', startingPos='-66.0 -115.0 135.0 0.0 0.7 0.0 0.7', speed='800.0', instruments='InterpolCatheter', threshold='0.1', name="IRController")
  instruments.createObject('WireBeamInterpolation', name='InterpolCatheter', WireRestShape='@../CatheterModel/catheterRestShape')
  instruments.createObject('AdaptiveBeamForceFieldAndMass', youngModulus='5', massDensity='0.000005', poisson='0.48', interpolation='@InterpolCatheter')
  instruments.createObject('LinearSolverConstraintCorrection', wire_optimization='true')
  instruments.createObject('FixedConstraint', indices='0')

  collis = instruments.createChild('Collision')
  collis.createObject('EdgeSetTopologyContainer')
  collis.createObject('MechanicalObject')
  collis.createObject('MultiAdaptiveBeamMapping', ircontroller='../IRController', useCurvAbs='1')
  collis.createObject('Line', proximity='0.1', group='1', selfCollision='1', LineActiverPath='../IRController')
  collis.createObject('Point', proximity='0.1', group='1', selfCollision='1', PointActiverPath='../IRController')

  return rootNode
