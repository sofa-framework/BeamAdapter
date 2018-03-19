/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
#include <regex>
#include <vector>
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

struct BeamAdapterTest : public Sofa_test<>,
        public ::testing::WithParamInterface<std::vector<std::string>>
{
    void simpleScene(const std::vector<std::string>& lines)
    {
        assert(lines.size()==3);
        string scene =
                "<?xml version='1.0'?>"
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'>"
                "   			<EulerImplicit rayleighStiffness='0.08' rayleighMass='0.08' printLog='false' />"
                "               <CGLinearSolver iterations='100' threshold='1e-10' tolerance='1e-15' />"
                "               $line1"
                "               <BeamInterpolation template='Rigid' name='Interpol' radius='0.1'/>"
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

std::vector<std::vector<std::string>> teststrings ={
    {
        "<Mesh name='meshSuture' edges='0 1' />"
        "<MechanicalObject template='Rigid' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ,""
        , "T"
    },
    {
        "<MechanicalObject template='Rigid' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ,"<Mesh name='meshSuture' edges='0 1' />"
        , "T"
    },
    {
        "<Mesh name='meshSuture' edges='0 1' />"
        ,"<MechanicalObject template='Rigid' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        , "W"
    },
    {
        "<MechanicalObject template='Rigid' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ,"<AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='1.0'/>"
        , "W"
    },
    {
        "<Mesh name='meshSuture' edges='0 1' />"
        ,"<AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='1.0'/>"
        , "W"
    }
};

TEST_P(BeamAdapterTest, checkMinimalScene) {
    ASSERT_NO_THROW(this->simpleScene(GetParam())) ;
}

INSTANTIATE_TEST_CASE_P(checkMinimalScene,
                        BeamAdapterTest, ::testing::ValuesIn(teststrings) ) ;

}
