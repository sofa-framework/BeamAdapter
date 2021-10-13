#include <sofa/testing/BaseSimulationTest.h>

#include <sofa/helper/BackTrace.h>

#include <SofaBaseUtils/initSofaBaseUtils.h>

#include <SofaBaseMechanics/MechanicalObject.h>
using sofa::core::topology::BaseMeshTopology ;
using sofa::core::objectmodel::Data ;

#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
using sofa::component::topology::TetrahedronSetTopologyContainer ;

using sofa::helper::WriteAccessor ;
using sofa::defaulttype::Rigid3dTypes ;

#include <SofaSimulationCommon/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML ;

#include <SofaSimulationGraph/DAGSimulation.h>
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;
using sofa::simulation::setSimulation ;
using sofa::core::objectmodel::New ;
using sofa::core::objectmodel::BaseData ;
using sofa::component::container::MechanicalObject ;

#include "../../../component/forcefield/AdaptiveBeamForceFieldAndMass.h"
using sofa::component::forcefield::AdaptiveBeamForceFieldAndMass;

#include <string>
using std::string;

namespace sofa
{

struct AdaptiveBeamForceFieldAndMassTest : public sofa::testing::BaseSimulationTest
{
    void simpleSceneTest(){
        sofa::component::initSofaBaseUtils();

        string scene =
                "<?xml version='1.0'?>"
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'>"
                "               <RequiredPlugin name='SofaBaseLinearSolver' />"
                "               <RequiredPlugin name='SofaImplicitOdeSolver' />"
                "               <RequiredPlugin name='SofaBoundaryCondition' />"
                "   			<EulerImplicitSolver rayleighStiffness='0.08' rayleighMass='0.08' printLog='false' />"
                "               <CGLinearSolver iterations='100' threshold='1e-10' tolerance='1e-15' />"
                "               <Mesh name='meshSuture' edges='0 1' />"
                "               <MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
                "               <BeamInterpolation name='Interpol' radius='0.1'/>"
                "               <AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='1.0'/>"
                "               <FixedConstraint indices='0' />"
                "</Node> " ;
        Node::SPtr root = SceneLoaderXML::loadFromMemory ( "test1", scene.c_str(), scene.size());

        ASSERT_NE(root.get(), nullptr);
        MechanicalObject<defaulttype::Rigid3dTypes>* mechanicalObject = nullptr;
        root->getTreeObject(mechanicalObject);

        ASSERT_NE(mechanicalObject, nullptr);
        EXPECT_TRUE(mechanicalObject->getName() == "DOFs") ;

        AdaptiveBeamForceFieldAndMass<defaulttype::Rigid3dTypes>* beamForceFieldMass  = nullptr;

        root->getTreeObject(beamForceFieldMass);

        ASSERT_NE(beamForceFieldMass, nullptr);
        ASSERT_NO_THROW(beamForceFieldMass->init());
        ASSERT_NO_THROW(beamForceFieldMass->reinit());
    }
};

TEST_F(AdaptiveBeamForceFieldAndMassTest, SimpleScene) {
    ASSERT_NO_THROW(this->simpleSceneTest()) ;
}

}
