import string
import numpy as np
import random as rd
import Sofa
import Sofa.Core
import Sofa.Simulation
import SofaRuntime
from math import sin,cos, sqrt, pi
import splib
from splib.animation import AnimationManager, addAnimation

## Register all the common component in the factory.
SofaRuntime.importPlugin("SofaAllCommonComponents")


class Animation(Sofa.Core.Controller):

    def __init__(self, InterventionalRadiology, controlledInstrument, step, angularStep):
        Sofa.Core.Controller.__init__(self)
        self.InterventionalRadiology = InterventionalRadiology
        
        self.controlledInstrument = controlledInstrument
        self.step = step
        self.angularStep = angularStep
    
        return;

    def init(self):
        self.startingPos = self.InterventionalRadiology.findData('startingPos').value
        self.nbInstruments = len(self.InterventionalRadiology.findData('instruments').value)
        self.totalLength = self.InterventionalRadiology.findData('totalLengths').value
        print ("totalLen : ", self.totalLength)


    def onEvent(self, event):
        pass

    def onKeypressedEvent(self, c):
        key = c['key']
        if key=="0": #
            self.controlledInstrument = 0
        if key=="1": #
            if(self.nbInstruments >= 2) :
                self.controlledInstrument = 1
            else:
                print("ERROR number of instruments is {}".format(self.nbInstruments))
        if key=="2": #
            if(self.nbInstruments >= 3) :
                self.controlledInstrument = 2
            else:
                print("ERROR number of instruments is {}".format(self.nbInstruments))

        if ord(key)==19: # up
            #@todo check we don't go further back than 0
            with self.InterventionalRadiology.xtip.writeable() as xTip:
                if (self.totalLength[self.controlledInstrument] >= xTip[self.controlledInstrument] + self.step - self.startingPos[0]) :
                    xTip[self.controlledInstrument] += self.step


        if ord(key)==21: # down
            with self.InterventionalRadiology.xtip.writeable() as xTip:
                if (self.totalLength[self.controlledInstrument] >= xTip[self.controlledInstrument] + self.step - self.startingPos[0]) :
                    xTip[self.controlledInstrument] -= self.step


        if ord(key)==18: # left
            with self.InterventionalRadiology.rotationInstrument.writeable() as rotationInstrument:
                rotationInstrument[self.controlledInstrument] -= self.angularStep

        if ord(key)==20: # right
            with self.InterventionalRadiology.rotationInstrument.writeable() as rotationInstrument:
                rotationInstrument[self.controlledInstrument] += self.angularStep


step = 0.1
angularStep = pi/20.0

def createScene(rootNode):

    rootNode.addObject('RequiredPlugin', name='SofaPython3')
    rootNode.addObject('RequiredPlugin', pluginName='BeamAdapter')

    rootNode.dt=0.01
    rootNode.addObject('VisualStyle', displayFlags='showVisualModels showBehaviorModels showCollisionModels hideMappings showForceFields showWireframe')
    
    rootNode.addObject('OglSceneFrame', style="Arrows", alignment="TopRight")

    rootNode.addObject('FreeMotionAnimationLoop')
    rootNode.addObject('LCPConstraintSolver', mu=0.1, tolerance=3e-4, maxIt=1000, build_lcp=False)
    rootNode.addObject('CollisionPipeline', draw=False, depth=6, verbose=False)
    rootNode.addObject('BruteForceDetection', name='N2')
    rootNode.addObject('LocalMinDistance', contactDistance=0.1, alarmDistance=0.3, name='localmindistance', angleCone=0.02)
    rootNode.addObject('CollisionResponse', name='Response', response='FrictionContact')

    manager = AnimationManager(rootNode)

    needle_0 = rootNode.addChild('needle_0')
    wire_0 = needle_0.addObject('WireRestShape', name='guideRestShape',
                             straightLength=0.0, length=500.0,
                             numEdges=250, youngModulus=20000,
                             spireDiameter=25, numEdgesCollis=[50, 10],
                             printLog=False, template='Rigid3d', spireHeight=0.0,
                             densityOfBeams=[30, 5], youngModulusExtremity=1000)
    needle_0.addObject('EdgeSetTopologyContainer', name='meshLinesGuide')
    needle_0.addObject('EdgeSetTopologyModifier', name='Modifier')
    needle_0.addObject('EdgeSetTopologyAlgorithms', name='TopoAlgo', template='Rigid3d')
    needle_0.addObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid3d')
    needle_0.addObject('MechanicalObject', name='dofTopo2', template='Rigid3d')

    needle_1 = rootNode.addChild('needle_1')
    wire_1 = needle_1.addObject('WireRestShape', name='guideRestShape',
                             straightLength=500.0, length=1000.0,
                             numEdges=250, youngModulus=20000,
                             spireDiameter=20, numEdgesCollis=[50, 10],
                             printLog=False, template='Rigid3d', spireHeight=0.0,
                             densityOfBeams=[30, 5], youngModulusExtremity=1000)
    needle_1.addObject('EdgeSetTopologyContainer', name='meshLinesGuide')
    needle_1.addObject('EdgeSetTopologyModifier', name='Modifier')
    needle_1.addObject('EdgeSetTopologyAlgorithms', name='TopoAlgo', template='Rigid3d')
    needle_1.addObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid3d')
    needle_1.addObject('MechanicalObject', name='dofTopo2', template='Rigid3d')
    
    needle_2 = rootNode.addChild('needle_2')
    wire_2 = needle_2.addObject('WireRestShape', name='guideRestShape',
                             straightLength=500.0, length=1000.0,
                             numEdges=250, youngModulus=20000,
                             spireDiameter=20, numEdgesCollis=[50, 10],
                             printLog=False, template='Rigid3d', spireHeight=0.0,
                             densityOfBeams=[30, 5], youngModulusExtremity=1000)
    needle_2.addObject('EdgeSetTopologyContainer', name='meshLinesGuide')
    needle_2.addObject('EdgeSetTopologyModifier', name='Modifier')
    needle_2.addObject('EdgeSetTopologyAlgorithms', name='TopoAlgo', template='Rigid3d')
    needle_2.addObject('EdgeSetGeometryAlgorithms', name='GeomAlgo', template='Rigid3d')
    needle_2.addObject('MechanicalObject', name='dofTopo2', template='Rigid3d')


    # instrumentsCombined = createInstrumentsCombined(rootNode)
    InstrumentCombined = rootNode.addChild('InstrumentCombined')
    InstrumentCombined.addObject('EulerImplicit', rayleighStiffness=0.2, printLog=False, rayleighMass=0.1)
    InstrumentCombined.addObject('BTDLinearSolver', verification=False, subpartSolve=False, verbose=False)
    InstrumentCombined.addObject('RegularGrid', name='meshLinesCombined', zmax=1, zmin=1, nx=60, ny=1, nz=1,
                                    xmax=1.0, xmin=0.0, ymin=0, ymax=0)
    InstrumentCombined.addObject('MechanicalObject', showIndices=False, name='DOFs', template='Rigid3d', ry=-90)
    
    InstrumentCombined.addObject('WireBeamInterpolation', WireRestShape='@../needle_0/guideRestShape',
                                    radius=0.15, printLog=False, name='Interpolguide')
    InstrumentCombined.addObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155,
                                    name='guideForceField', interpolation='@Interpolguide')
    
    InstrumentCombined.addObject('WireBeamInterpolation', WireRestShape='@../needle_1/guideRestShape',
                                    radius=0.15, printLog=False, name='Interpolguide_1')
    InstrumentCombined.addObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155,
                                    name='guideForceField_1', interpolation='@Interpolguide_1')
    
    InstrumentCombined.addObject('WireBeamInterpolation', WireRestShape='@../needle_2/guideRestShape',
                                    radius=0.15, printLog=False, name='Interpolguide_2')
    InstrumentCombined.addObject('AdaptiveBeamForceFieldAndMass', massDensity=0.00000155,
                                    name='guideForceField_2', interpolation='@Interpolguide_2')
    

    intrevention = InstrumentCombined.addObject('InterventionalRadiology', xtip=[0.0, 0.1, 0.1], name='m_ircontroller',
                                    instruments=['Interpolguide','Interpolguide_1','Interpolguide_2'], step=0.5, printLog=False,
                                    listening=True, template='Rigid3d', startingPos=[0, 0, 0, 1, 0, 0, 0],
                                    rotationInstrument=[0, 0, 0.1], speed=0, controlledInstrument=0)
        
    InstrumentCombined.addObject('LinearSolverConstraintCorrection', wire_optimization='true', printLog=False)
    InstrumentCombined.addObject('FixedConstraint', indices=0, name='FixedConstraint')
    InstrumentCombined.addObject('RestShapeSpringsForceField', points='@m_ircontroller.indexFirstNode',
                                    angularStiffness=1e8, stiffness=1e8)
    InstrumentCombined.addObject(Animation(intrevention, 0, step, angularStep))

    return(rootNode)

