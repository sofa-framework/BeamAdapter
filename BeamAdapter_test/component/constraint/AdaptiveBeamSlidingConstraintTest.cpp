#include <string>
using std::string ;

#include <SofaTest/Sofa_test.h>
#include <sofa/helper/BackTrace.h>
#include <SofaBaseMechanics/MechanicalObject.h>

#include <SofaBaseLinearSolver/FullVector.h>
using sofa::core::topology::BaseMeshTopology ;
using sofa::core::objectmodel::Data ;

using sofa::helper::WriteAccessor ;
using sofa::defaulttype::Rigid3dTypes ;

#include <SofaSimulationCommon/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML ;

#include <SofaSimulationGraph/DAGSimulation.h>
using sofa::simulation::graph::DAGSimulation;
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;
using sofa::simulation::setSimulation ;
using sofa::core::objectmodel::New ;
using sofa::core::objectmodel::BaseData ;
using sofa::component::container::MechanicalObject ;

#include "component/constraint/AdaptiveBeamSlidingConstraint.h"
#include "component/WireBeamInterpolation.h"
using sofa::component::constraintset::AdaptiveBeamSlidingConstraint ;
using sofa::component::fem::WireBeamInterpolation ;


namespace sofa
{

template <typename DataTypes>
struct AdaptiveBeamSlidingConstraintTest : public Sofa_test<typename DataTypes::Real>, AdaptiveBeamSlidingConstraint<DataTypes>
{
    void normalBehavior(){
        Simulation* simu;
        setSimulation(simu = new DAGSimulation());

        typename AdaptiveBeamSlidingConstraint<DataTypes>::SPtr thisObject = New<AdaptiveBeamSlidingConstraint<DataTypes>>();
        thisObject->setName("myname");
        EXPECT_TRUE(thisObject->getName() == "myname");
        EXPECT_TRUE(thisObject->findLink("interpolation") != nullptr );

        Node::SPtr node = simu->createNewGraph("root");
        typename MechanicalObject<DataTypes>::SPtr mecaobject = New<MechanicalObject<DataTypes> >();
        typename WireBeamInterpolation<DataTypes>::SPtr interpolation = New<WireBeamInterpolation<DataTypes> >();

        interpolation->findData("name")->read("wireInterpolation");
        node->addObject(mecaobject);
        node->addObject(interpolation);
        node->addObject(thisObject);
        thisObject->findLink("interpolation")->read("@./wireInterpolation");

        EXPECT_NO_THROW( thisObject->init() );
        EXPECT_NO_THROW( thisObject->reset() );
    }

};

using testing::Types;
typedef Types<Rigid3dTypes> DataTypes;

TYPED_TEST_CASE(AdaptiveBeamSlidingConstraintTest, DataTypes);

TYPED_TEST(AdaptiveBeamSlidingConstraintTest, NormalBehavior) {
    this->normalBehavior() ;
}

}

