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
using sofa::defaulttype::Vec3Types ;

#include <SofaSimulationCommon/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML ;

#include <SofaSimulationGraph/DAGSimulation.h>
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;
using sofa::simulation::setSimulation ;
using sofa::core::objectmodel::New ;
using sofa::core::objectmodel::BaseData ;
using sofa::component::container::MechanicalObject ;




namespace sofa
{

template <typename _DataTypes>
struct BeamAdapterFirstTest : public Sofa_test<typename _DataTypes::Real>
{

    typedef _DataTypes DataTypes;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;



    void normalTests(){
        Simulation* simu;
        setSimulation(simu = new sofa::simulation::graph::DAGSimulation());

        Node::SPtr node = simu->createNewGraph("root");
        typename MechanicalObject<DataTypes>::SPtr mecaobject = New<MechanicalObject<DataTypes> >() ;
        mecaobject->init() ;

        node->addObject(mecaobject) ;

        mecaobject->setName("myname") ;
        EXPECT_TRUE(mecaobject->getName() == "myname") ;


        return ;
    }

    void simpleSceneTest(){
        string scene =
                "<?xml version='1.0'?>"
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'> "
                "   <MechanicalObject template/>
                "</Node> " ;
        EXPECT_NO_THROW(SceneLoaderXML::loadFromMemory ( "test1", scene.c_str(), scene.size())) ;
    }

    double ComputationTest()
    {
        double toto=3.0;
        return toto;

    }

};

using testing::Types;
typedef Types<Vec3Types> DataTypes;

TYPED_TEST_CASE(BeamAdapterFirstTest, DataTypes);


TYPED_TEST(BeamAdapterFirstTest, NormalBehavior) {
    ASSERT_NO_THROW(this->normalTests()) ;
}

TYPED_TEST(BeamAdapterFirstTest, SimpleScene) {
    ASSERT_NO_THROW(this->simpleSceneTest()) ;
}

TYPED_TEST(BeamAdapterFirstTest, ComputationTest) {
    ASSERT_DOUBLE_EQ(12.5,this->ComputationTest());
}
}
