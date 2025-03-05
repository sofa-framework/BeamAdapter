/******************************************************************************
*                              BeamAdapter plugin                             *
*                  (c) 2006 Inria, University of Lille, CNRS                  *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: see Authors.md                                                     *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/testing/BaseSimulationTest.h>

#include <sofa/helper/BackTrace.h>
#include <sofa/core/ExecParams.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/core/ExecParams.h>

using sofa::core::topology::BaseMeshTopology ;
using sofa::core::objectmodel::Data ;

#include <sofa/component/topology/container/dynamic/TetrahedronSetTopologyContainer.h>
using sofa::component::topology::container::dynamic::TetrahedronSetTopologyContainer ;

using sofa::helper::WriteAccessor ;
using sofa::defaulttype::Rigid3dTypes ;
using sofa::core::ExecParams ;

#include <sofa/simulation/common/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML;

#include <sofa/simulation/graph/DAGSimulation.h>
using sofa::simulation::graph::DAGSimulation;
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;
using sofa::core::objectmodel::New ;
using sofa::core::objectmodel::BaseData ;
using sofa::component::statecontainer::MechanicalObject ;

#include <sofa/simpleapi/SimpleApi.h>

#include <regex>
#include <vector>
#include <string>
using std::string;

#include <BeamAdapter/component/forcefield/AdaptiveBeamForceFieldAndMass.h>
#include <BeamAdapter/component/BeamInterpolation.h>
using sofa::component::fem::_beaminterpolation_::BeamInterpolation;

namespace sofa
{

struct BeamInterpolationTest : public  sofa::testing::BaseSimulationTest,
        public ::testing::WithParamInterface<std::vector<std::string>>
{
    void doSetUp() override
    {
        sofa::simpleapi::importPlugin(Sofa.Component.ODESolver.Backward);
        sofa::simpleapi::importPlugin(Sofa.Component.LinearSolver.Iterative);
        sofa::simpleapi::importPlugin(Sofa.Component.StateContainer);
        sofa::simpleapi::importPlugin(Sofa.Component.Topology.Container.Constant);
        sofa::simpleapi::importPlugin("BeamAdapter");
    }
    
    Node::SPtr createSingleBeam()
    {
        string scene =
            "<?xml version='1.0'?>"
            "<Node 	name='Root' gravity='0 -9.81 0' dt='0.01'>"
            "    <RequiredPlugin name='Sofa.Component.ODESolver.Backward' />"
            "    <RequiredPlugin name='Sofa.Component.LinearSolver.Direct' />"
            "    <RequiredPlugin name='Sofa.Component.Constraint.Projective' />"
            "    <RequiredPlugin name='Sofa.Component.StateContainer' />"
            "    <RequiredPlugin name='Sofa.Component.Topology.Container.Constant' />"
            "    <RequiredPlugin name='Sofa.Component.Topology.Container.Grid' />"
            "    <RequiredPlugin name='BeamAdapter' />"
            "    <DefaultAnimationLoop />"
            "    <DefaultVisualManagerLoop />"
            "    <Node name='BeamModel'/>"
            "        <EulerImplicitSolver rayleighStiffness='0.0' rayleighMass='0.0' />"
            "        <BTDLinearSolver />"
            "        <RegularGridTopology name='MeshLines' nx='200' ny='1' nz='1' xmax='100' xmin='0' ymin='0' ymax='0' zmax='0' zmin='0'/>"
            "        <MechanicalObject template='Rigid3d' name='DOFs' />"
            "        <FixedConstraint indices='0' />"
            "        <BeamInterpolation name='Interpol' radius='0.1'/>"
            "    </Node> "
            "</Node> ";

        Node::SPtr root = SceneLoaderXML::loadFromMemory("singleBeam", scene.c_str());
        sofa::simulation::node::initRoot(root.get());

        return root;
    }

    void checkDataInitialization(const Data<type::vector<SReal>>& data,
                                 const SReal& defaultValue,
                                 const sofa::Size& nbBeam,
                                 const SReal& value)
    {
        const auto& vector = helper::getReadAccessor(data);
        ASSERT_EQ(vector.size(), nbBeam);
        for (const auto& v: vector)
            ASSERT_FLOAT_EQ(v, value);
        ASSERT_FLOAT_EQ(defaultValue, value);
    }

    void checkCreation()
    {
        Node::SPtr root = createSingleBeam();

        // Search for Beam FF
        BeamInterpolation<Rigid3dTypes>* beamInterpolation = nullptr;
        root->getTreeObject(beamInterpolation);
        ASSERT_NE(beamInterpolation, nullptr);

        sofa::Size nbBeam = 199;

        // Check component state and Data default values
        checkDataInitialization(beamInterpolation->d_radius, beamInterpolation->m_defaultRadius, nbBeam, 0.1);
        checkDataInitialization(beamInterpolation->d_innerRadius, beamInterpolation->m_defaultInnerRadius, nbBeam, 0.0);
        checkDataInitialization(beamInterpolation->d_lengthY, beamInterpolation->m_defaultLengthY, nbBeam, 1.0);
        checkDataInitialization(beamInterpolation->d_lengthZ, beamInterpolation->m_defaultLengthZ, nbBeam, 1.0);
        checkDataInitialization(beamInterpolation->d_defaultYoungModulus, beamInterpolation->m_defaultYoungModulus, nbBeam, 1e5);
        checkDataInitialization(beamInterpolation->d_poissonRatio, beamInterpolation->m_defaultPoissonRatio, nbBeam, 0.4);
    }

    void simpleScene(const std::vector<std::string>& lines)
    {
        assert(lines.size()==3);
        string scene =
                "<?xml version='1.0'?>"
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'>"
                "   		    <EulerImplicitSolver rayleighStiffness='0.08' rayleighMass='0.08' printLog='false' />"
                "               <CGLinearSolver iterations='100' threshold='1e-10' tolerance='1e-15' />"
                "               $line1"
                "               <BeamInterpolation template='Rigid3d' name='Interpol' radius='0.1'/>"
                "               $line2"
                "</Node> " ;

        scene = std::regex_replace(scene, std::regex("\\$line1"), lines[0]) ;
        scene = std::regex_replace(scene, std::regex("\\$line2"), lines[1]) ;

        if(lines[2]=="T")
        {
            EXPECT_MSG_NOEMIT(Error, Warning) ;

            Node::SPtr root = SceneLoaderXML::loadFromMemory ( "test1", scene.c_str());
            ASSERT_NE(root.get(), nullptr);

            root->init(core::ExecParams::defaultInstance());
            root->reinit(core::ExecParams::defaultInstance()) ;
        }
        else if(lines[2]=="W")
        {
            EXPECT_MSG_EMIT(Error) ;
            EXPECT_MSG_NOEMIT(Warning) ;

            Node::SPtr root = SceneLoaderXML::loadFromMemory ( "test1", scene.c_str());
            ASSERT_NE(root.get(), nullptr);

            root->init(core::ExecParams::defaultInstance());
            root->reinit(core::ExecParams::defaultInstance()) ;
        }
    }
};


static std::vector<std::vector<std::string>> teststrings ={
    {
        "<MeshTopology name='meshSuture' edges='0 1' />"
        "<MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ,""
        , "T"
    },
    {
        "<MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ,"<MeshTopology name='meshSuture' edges='0 1' />"
        , "T"
    },
    {
        "<MeshTopology name='meshSuture' edges='0 1' />"
        ,"<MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        , "W"
    },
    {
        "<MechanicalObject template='Rigid3d' name='DOFs' showIndices='0' position='0 0 0 0 0 0 1   1 0 0 0 0 0 1'/>"
        ,"<AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='1.0'/>"
        , "W"
    },
    {
        "<MeshTopology name='meshSuture' edges='0 1' />"
        ,"<AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='1.0'/>"
        , "W"
    }
};

TEST_P(BeamInterpolationTest, checkMinimalScene) {
    ASSERT_NO_THROW(this->simpleScene(GetParam())) ;
}

TEST_F(BeamInterpolationTest, checkCreation) {
    ASSERT_NO_THROW(this->checkCreation());
}

INSTANTIATE_TEST_SUITE_P(checkMinimalScene,
                        BeamInterpolationTest, ::testing::ValuesIn(teststrings) ) ;

}
