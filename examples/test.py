import Sofa
import BeamAdapter 

def BeamModel(root):
        """Creates a BeamModel"""
        b = root.createChild("BeamModel")
        b.createObject("MechanicalObject", template="Rigid3", name="dofs")
        b.createObject("EdgeSetTopologyContainer", name="container", 
                              position=[[0.0,0.0,0.0], [1.0,0.0,0.0], [2.0,0.0,0.0], [3.0,0.0,0.0], [4.0,0.0,0.0]],
                              edges=[[0,1],[1,2],[2,3],[3,4]])
        b.createObject("WireBeamInterpolation", crossSectionShape="circular", name="interpolation", radius=2, printLog=True)
        b.init()
	return b

def createScene(root):
        beam = BeamModel(root)
        beam.dofs.showObject=True        
        beam.dofs.showObjectScale=0.7
	interpolatedValues = beam.interpolation.getValuesAt([(0,0.0),(0,0.5),(0,1.0)])
        print("Values: ", interpolatedValues);
        print("Rest Total Length: ", beam.interpolation.getRestTotalLength())
        print("Length: ", beam.interpolation.getLength(2))

