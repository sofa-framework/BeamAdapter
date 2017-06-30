#include <string>
using std::string ;

#include <SofaTest/Sofa_test.h>
#include <sofa/helper/BackTrace.h>
#include <SofaBaseMechanics/MechanicalObject.h>

#include <SofaBaseLinearSolver/FullVector.h>
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

namespace sofa
{

struct AdaptiveBeamForceFieldAndMassTest : public Sofa_test<>
{
    void simpleSceneTest(){
        string scene =
                "<?xml version='1.0'?>"
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'>"
                "   			<EulerImplicit rayleighStiffness='0.08' rayleighMass='0.08' printLog='false' />"
                "               <CGLinearSolver iterations='100' threshold='1e-10' tolerance='1e-15' />"
                "               <Mesh name='meshSuture' edges='0 1' />"
                "               <MechanicalObject template='Rigid' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
                "               <BeamInterpolation name='Interpol' radius='0.1'/>"
                "               <AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='1.0'/>"
                "               <FixedConstraint indices='0' />"
                "</Node> " ;
        Node::SPtr root = SceneLoaderXML::loadFromMemory ( "test1", scene.c_str(), scene.size());

        ASSERT_NE(root.get(), nullptr);
        MechanicalObject<Rigid3>* mechanicalObject = nullptr;
        root->getTreeObject(mechanicalObject);

        ASSERT_NE(mechanicalObject, nullptr);
        EXPECT_TRUE(mechanicalObject->getName() == "DOFs") ;

        AdaptiveBeamForceFieldAndMass<Rigid3>* beamForceFieldMass  = nullptr;

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
