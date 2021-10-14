#include <sofa/testing/BaseSimulationTest.h>

#include <regex>
#include <vector>
#include <string>
using std::string ;

#include <sofa/helper/BackTrace.h>

#include <SofaBaseUtils/initSofaBaseUtils.h>

#include <SofaBaseMechanics/MechanicalObject.h>

using sofa::core::topology::BaseMeshTopology ;
using sofa::core::objectmodel::Data ;

#include <SofaBaseTopology/TetrahedronSetTopologyContainer.h>
using sofa::component::topology::TetrahedronSetTopologyContainer ;

using sofa::helper::WriteAccessor ;
using sofa::defaulttype::Rigid3dTypes ;
using sofa::core::ExecParams ;

#include <SofaSimulationCommon/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML ;

#include <SofaSimulationGraph/DAGSimulation.h>
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;
using sofa::simulation::setSimulation ;
using sofa::core::objectmodel::New ;
using sofa::core::objectmodel::BaseData ;
using sofa::component::container::MechanicalObject ;

#include "../component/forcefield/AdaptiveBeamForceFieldAndMass.h"

namespace sofa
{

struct WireBeamInterpolationTest : public  sofa::testing::BaseSimulationTest,
        public ::testing::WithParamInterface<std::vector<std::string>>
{
    void simpleScene(const std::vector<std::string>& lines)
    {
        assert(lines.size()==3);
        sofa::component::initSofaBaseUtils();

        string scene =
                "<?xml version='1.0'?>"
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'>"
                "               <RequiredPlugin name='SofaBaseLinearSolver' />"
                "               <RequiredPlugin name='SofaImplicitOdeSolver' />"
                "   			<EulerImplicitSolver rayleighStiffness='0.08' rayleighMass='0.08' printLog='false' />"
                "               <CGLinearSolver iterations='100' threshold='1e-10' tolerance='1e-15' />"
                "               <EdgeSetTopologyContainer/> "
                "               <EdgeSetTopologyModifier /> "
                "               <EdgeSetGeometryAlgorithms template='Rigid3d' /> "
                "               $line1"
                "               <WireBeamInterpolation template='Rigid3d' name='Interpol' WireRestShape='@restShape' radius='0.1'/>"
                "               $line2"
                "</Node> " ;

        scene = std::regex_replace(scene, std::regex("\\$line1"), lines[0]) ;
        scene = std::regex_replace(scene, std::regex("\\$line2"), lines[1]) ;

        if(lines[2]=="T")
        {
            EXPECT_MSG_NOEMIT(Error, Warning) ;

            Node::SPtr root = SceneLoaderXML::loadFromMemory ( "test1", scene.c_str(), scene.size());
            ASSERT_NE(root.get(), nullptr);

            root->init(ExecParams::defaultInstance());
            root->reinit(ExecParams::defaultInstance()) ;
        }else if(lines[2]=="W")
        {
            EXPECT_MSG_EMIT(Error) ;
            EXPECT_MSG_NOEMIT(Warning) ;

            Node::SPtr root = SceneLoaderXML::loadFromMemory ( "test1", scene.c_str(), scene.size());
            ASSERT_NE(root.get(), nullptr);

            root->init(ExecParams::defaultInstance());
            root->reinit(ExecParams::defaultInstance()) ;
        }
    }
};

static std::vector<std::vector<std::string>> teststrings ={
    {
        "<Mesh name='meshSuture' edges='0 1' />"
        "<WireRestShape template='Rigid3d' name='restShape' length='600.0'  straightLength='400' spireDiameter='4000.0' spireHeight='0.0' densityOfBeams='15 0' numEdges='200' numEdgesCollis='40 40'  youngModulus='1000' youngModulusExtremity='1000'/>"
        "<MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ""
        ,""
        , "T"
    },
    {
        "<MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ,"<Mesh name='meshSuture' edges='0 1' />"
        , "T"
    },
    {
        "<Mesh name='meshSuture' edges='0 1' />"
        ,"<MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        , "W"
    },
    {
        "<MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ,"<AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='1.0'/>"
        , "W"
    },
    {
        "<Mesh name='meshSuture' edges='0 1' />"
        ,"<AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='1.0'/>"
        , "W"
    }
};

TEST_P(WireBeamInterpolationTest, checkMinimalScene) {
    ASSERT_NO_THROW(this->simpleScene(GetParam())) ;
}

INSTANTIATE_TEST_CASE_P(checkMinimalScene,
                        WireBeamInterpolationTest, ::testing::ValuesIn(teststrings) ) ;

}
