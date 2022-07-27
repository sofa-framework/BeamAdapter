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

#include <sofa/component/statecontainer/MechanicalObject.h>
using sofa::core::topology::BaseMeshTopology ;
using sofa::core::objectmodel::Data ;

#include <sofa/component/topology/container/dynamic/TetrahedronSetTopologyContainer.h>
using sofa::component::topology::container::dynamic::TetrahedronSetTopologyContainer;

using sofa::helper::WriteAccessor ;
using sofa::defaulttype::Rigid3dTypes ;

#include <sofa/simulation/common/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML;

#include <sofa/simulation/graph/DAGSimulation.h>
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;
using sofa::simulation::setSimulation ;
using sofa::core::objectmodel::New ;
using sofa::core::objectmodel::BaseData ;
using sofa::component::statecontainer::MechanicalObject ;

#include <BeamAdapter/component/forcefield/AdaptiveBeamForceFieldAndMass.h>
using sofa::component::forcefield::AdaptiveBeamForceFieldAndMass;

#include <string>
using std::string;

namespace sofa
{

struct AdaptiveBeamForceFieldAndMassTest : public sofa::testing::BaseSimulationTest
{
    void simpleSceneTest(){
        string scene =
                "<?xml version='1.0'?>"
                "<Node 	name='Root' gravity='0 0 0' time='0' animate='0'>"
                "               <RequiredPlugin name='Sofa.Component.ODESolver.Backward' />"
                "               <RequiredPlugin name='Sofa.Component.LinearSolver.Iterative' />"
                "               <RequiredPlugin name='Sofa.Component.StateContainer' />"
                "               <RequiredPlugin name='Sofa.Component.Constraint.Projective' />"
                "               <RequiredPlugin name='Sofa.Component.Topology.Container.Constant' />"
                "               <RequiredPlugin name='BeamAdapter' />"
                "   			<EulerImplicitSolver rayleighStiffness='0.08' rayleighMass='0.08' printLog='false' />"
                "               <CGLinearSolver iterations='100' threshold='1e-10' tolerance='1e-15' />"
                "               <MeshTopology name='meshSuture' edges='0 1' />"
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
