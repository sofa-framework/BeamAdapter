#include <string>
using std::string ;

#include <SofaTest/Sofa_test.h>
#include <sofa/helper/BackTrace.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <vector>

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

using std::vector;
using sofa::core::ExecParams;


#include "WireRestShape.h"

namespace sofa
{





class WireRestShapeTest : public Sofa_test<>,
        public::testing::WithParamInterface< std::vector<std::string> >
{
public:
    void simpleSceneTest(const std::vector< std::string >& v){
        EXPECT_EQ (v.size(),2);

        std::stringstream scene ;
        scene<<
        "<?xml version='1.0'?>"
        "<Node name='topoLines_cath'>"
        "       <WireRestShape template='Rigid' printLog='0' name='catheterRestShape' length='1000.0' straightLength='600' spireDiameter='4000.0' spireHeight='0.0' densityOfBeams='40 10' numEdges='200' numEdgesCollis='40 20'  youngModulus='10000' youngModulusExtremity='10000' edge2QuadMappingName='"<< v[0] << "' /> "
        "       <EdgeSetTopologyContainer name='meshLinesCath' />"
        "       <EdgeSetTopologyModifier   name='Modifier' /> "
        "       <EdgeSetTopologyAlgorithms name='TopoAlgo'   template='Rigid' />"
        "       <EdgeSetGeometryAlgorithms name='GeomAlgo'   template='Rigid' />"
        "       <MechanicalObject template='Rigid' name='dofTopo1' />"
        "       <Node name='mappedTopo'>"
        "           <MechanicalObject name='Quads'/>"
        "           <QuadSetTopologyContainer  name='ContainerCoils' />"
        "           <QuadSetTopologyModifier   name='Modifier' /> "
        "           <QuadSetTopologyAlgorithms name='TopoAlgo'  template='Vec3d' />"
        "           <QuadSetGeometryAlgorithms name='GeomAlgo'  template='Vec3d' />"
        "           <Edge2QuadTopologicalMapping name='topoMap' nbPointsOnEachCircle='10' radius='0.3'  input='@../meshLinesCath' output='@ContainerCoils' />    "
        "       </Node> "
        " </Node>	" ;

        Node::SPtr root = SceneLoaderXML::loadFromMemory ( "test1", scene.str().c_str(), scene.str().size());

        ASSERT_NE(root.get(), nullptr);

        MechanicalObject<Rigid3>* MO = nullptr;

        root->getTreeObject(MO);

        ASSERT_NE(MO, nullptr);

        EXPECT_TRUE(MO->getName() == "dofTopo1") ;
        component::engine::WireRestShape<Rigid3>* WRS  = nullptr;
        root->getTreeObject(WRS);

        //EXPECT_TRUE( WRS->getMessageLog(Error).size() ==0 )

        if(v[1]=="t"){
            EXPECT_MSG_NOEMIT(Error);
            root->init(ExecParams::defaultInstance()) ;
            root->bwdInit() ;
        }else{
            EXPECT_MSG_EMIT(Error);
            root->init(ExecParams::defaultInstance()) ;
            root->bwdInit() ;
        }

        //EXPECT_TRUE( WRS->getMessageLog(Error).size() !=0 )

        ASSERT_NE(WRS, nullptr);
    }
};

std::vector< std::vector<std::string> > unintvalues = {
    {"", "t"}, {"toto", "f"}, {"./mappedTopo/topoMap", "t" }
};


TEST_P(WireRestShapeTest, SimpleSceneValid) {
    ASSERT_NO_THROW(this->simpleSceneTest(GetParam() ) );
}

INSTANTIATE_TEST_CASE_P(checkRestShapeInits,
                        WireRestShapeTest,
                        ::testing::ValuesIn(unintvalues));


}
