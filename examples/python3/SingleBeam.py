import sys
import Sofa

def createScene(rootNode):

    # rootNode
    rootNode.addObject('RequiredPlugin', pluginName='BeamAdapter Sofa.Component.Constraint.Projective Sofa.Component.LinearSolver.Direct Sofa.Component.ODESolver.Backward Sofa.Component.StateContainer Sofa.Component.Topology.Container.Constant Sofa.Component.Topology.Container.Grid Sofa.Component.Visual')
    rootNode.addObject('VisualStyle', displayFlags='showBehaviorModels showCollisionModels hideBoundingCollisionModels showForceFields')
    rootNode.addObject('DefaultAnimationLoop')
    rootNode.addObject('DefaultVisualManagerLoop')

    # rootNode/BeamModel
    BeamModel = rootNode.addChild('BeamModel', bbox="-3 -6 -3 3 3 3")
    BeamModel.addObject('EulerImplicitSolver', rayleighStiffness='0', rayleighMass='0', printLog='false')
    BeamModel.addObject('BTDLinearSolver', verbose='0')
    BeamModel.addObject('RegularGridTopology', name='MeshLines', n=[200, 1, 1], min=[0, 0, 0], max=[100, 0, 0])
    BeamModel.addObject('MechanicalObject', template='Rigid3d', name='DOFs', position='@MeshLines.positions')
    BeamModel.addObject('FixedProjectiveConstraint', name='FixedConstraint', indices='0')
    BeamModel.addObject('BeamInterpolation', name='BeamInterpolation', radius='0.1')
    BeamModel.addObject('AdaptiveBeamForceFieldAndMass', name='BeamForceField', computeMass='1', massDensity='10')

    return 0;

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
