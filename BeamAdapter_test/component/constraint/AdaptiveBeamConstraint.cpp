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
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;
using sofa::simulation::setSimulation ;
using sofa::core::objectmodel::New ;
using sofa::core::objectmodel::BaseData ;
using sofa::component::container::MechanicalObject ;

#include "component/constraint/AdaptiveBeamSlidingConstraint.h"
using sofa::component::constraintset::AdaptiveBeamSlidingConstraint ;


namespace sofa
{

template <typename _DataTypes>
struct AdaptiveBeamSlidingConstraintTest : public Sofa_test<typename _DataTypes::Real>,
        AdaptiveBeamSlidingConstraint<_DataTypes>
{
    typedef _DataTypes DataTypes;
    typedef AdaptiveBeamSlidingConstraint<_DataTypes> ThisClass ;

    void normalBehavior(){
        Simulation* simu;
        setSimulation(simu = new sofa::simulation::graph::DAGSimulation());

        typename ThisClass::SPtr thisobject = New<ThisClass>() ;
        thisobject->setName("myname") ;
        EXPECT_TRUE(thisobject->getName() == "myname") ;
        EXPECT_TRUE( thisobject->findLink("interpolation") != nullptr ) ;
    }

};

using testing::Types;
typedef Types<Rigid3dTypes> DataTypes;

TYPED_TEST_CASE(AdaptiveBeamSlidingConstraintTest, DataTypes);

TYPED_TEST(AdaptiveBeamSlidingConstraintTest, NormalBehavior) {
    this->normalBehavior() ;
}

}

