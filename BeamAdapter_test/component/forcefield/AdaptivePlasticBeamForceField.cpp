#include <gtest/gtest.h>
#include <sofa/testing/BaseSimulationTest.h>

#include <SofaSimulationGraph/DAGSimulation.h>
#include <SofaSimulationGraph/SimpleApi.h>

#include "../../../component/forcefield/AdaptivePlasticBeamForceField.h"

using namespace sofa::simulation;
using namespace sofa::simpleapi;

using sofa::plugin::beamadapter::component::forcefield::_adaptiveplasticbeamforcefield_::AdaptivePlasticBeamForceField;
using sofa::defaulttype::Rigid3dTypes;


typedef AdaptivePlasticBeamForceField<Rigid3dTypes> AdaptivePlaticBeamForceField_Rigid;

TEST(Quadrature, gaussPointInitialisation)
{
	EXPECT_MSG_NOEMIT(Error);
	setSimulation(new sofa::simulation::graph::DAGSimulation());
	auto root = getSimulation()->createNewNode("root");

	createObject(root, "RequiredPlugin", { {"pluginName", "BeamAdapter"} });

	createObject(root, "RegularGridTopology", { {"name", "grid"}, {"min", "-7.5 -7.5 0"}, {"max", "7.5 7.5 80"}, {"n", "3 3 9"} });

	auto beam = createChild(root, "beam");

	getSimulation()->unload(root);
}