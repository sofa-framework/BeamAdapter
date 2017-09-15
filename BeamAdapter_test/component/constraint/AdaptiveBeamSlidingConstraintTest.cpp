#include <string>
using std::string ;

#include <SofaTest/Sofa_test.h>
#include <sofa/helper/BackTrace.h>
#include <SofaBaseMechanics/MechanicalObject.h>
using sofa::core::behavior::MechanicalState;

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


#include <sofa/core/behavior/PairInteractionConstraint.h>
using sofa::core::behavior::PairInteractionConstraint ;


namespace sofa
{

template <typename DataTypes>
struct AdaptiveBeamSlidingConstraintTest : public Sofa_test<typename DataTypes::Real>, AdaptiveBeamSlidingConstraint<DataTypes>
{
    Simulation* m_simu;
    typename AdaptiveBeamSlidingConstraint<DataTypes>::SPtr m_thisObject;
    Node::SPtr m_node;

    void SetUp()
    {
        setSimulation(m_simu = new DAGSimulation());

        m_thisObject = New<AdaptiveBeamSlidingConstraint<DataTypes>>();
        m_thisObject->setName("myname");
        EXPECT_TRUE(m_thisObject->getName() == "myname");
        EXPECT_TRUE(m_thisObject->findLink("interpolation") != nullptr );

        m_node = m_simu->createNewGraph("root");
        typename MechanicalObject<DataTypes>::SPtr mecaobject = New<MechanicalObject<DataTypes> >();
        typename WireBeamInterpolation<DataTypes>::SPtr interpolation = New<WireBeamInterpolation<DataTypes> >();

        interpolation->findData("name")->read("wireInterpolation");
        m_node->addObject(mecaobject);
        m_node->addObject(interpolation);
        m_node->addObject(m_thisObject);
        m_thisObject->findLink("interpolation")->read("@./wireInterpolation");
    }

    void normalBehavior()
    {
        EXPECT_NO_THROW( m_thisObject->init() );
        EXPECT_NO_THROW( m_thisObject->reset() );
    }

    void testInternalInit()
    {

    }

};



using testing::Types;
typedef Types<Rigid3dTypes> DataTypes;

TYPED_TEST_CASE(AdaptiveBeamSlidingConstraintTest, DataTypes);

TYPED_TEST(AdaptiveBeamSlidingConstraintTest, NormalBehavior) {
    this->normalBehavior() ;
}

TYPED_TEST(AdaptiveBeamSlidingConstraintTest, testInternalInit) {
    this->testInternalInit() ;
}

}

