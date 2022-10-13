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

#include <BeamAdapter/component/forcefield/AdaptiveBeamForceFieldAndMass.h>
#include <sofa/component/statecontainer/MechanicalObject.h>

#include <sofa/simulation/common/SceneLoaderXML.h>
#include <sofa/helper/system/thread/CTime.h>

#include <string>
using std::string;

namespace sofa
{
using sofa::component::forcefield::AdaptiveBeamForceFieldAndMass;
using sofa::component::statecontainer::MechanicalObject;
using sofa::helper::system::thread::ctime_t;

using sofa::simulation::Node;
using sofa::simulation::SceneLoaderXML;

struct AdaptiveBeamForceFieldAndMassTest : public sofa::testing::BaseSimulationTest
{
    using Rigid3dTypes = defaulttype::Rigid3dTypes;
    using VecCoord = MechanicalObject<Rigid3dTypes>::VecCoord;

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
        MechanicalObject<Rigid3dTypes>* mechanicalObject = nullptr;
        root->getTreeObject(mechanicalObject);

        ASSERT_NE(mechanicalObject, nullptr);
        EXPECT_TRUE(mechanicalObject->getName() == "DOFs") ;

        AdaptiveBeamForceFieldAndMass<Rigid3dTypes>* beamForceFieldMass  = nullptr;

        root->getTreeObject(beamForceFieldMass);

        ASSERT_NE(beamForceFieldMass, nullptr);
        ASSERT_NO_THROW(beamForceFieldMass->init());
        ASSERT_NO_THROW(beamForceFieldMass->reinit());
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
            "        <AdaptiveBeamForceFieldAndMass name='ForceField' interpolation='@Interpol' massDensity='10.0'/>"
            "    </Node> "
            "</Node> ";            

        Node::SPtr root = SceneLoaderXML::loadFromMemory("singleBeam", scene.c_str(), (unsigned int)(scene.size()));
        sofa::simulation::getSimulation()->init(root.get());

        return root;
    }


    void checkCreation()
    {
        Node::SPtr root = createSingleBeam();

        // Search for Beam FF
        AdaptiveBeamForceFieldAndMass<Rigid3dTypes>* beamForceFieldMass = nullptr;
        root->getTreeObject(beamForceFieldMass);
        ASSERT_NE(beamForceFieldMass, nullptr);

        // Check component state and Data default values
        ASSERT_EQ(beamForceFieldMass->d_componentState.getValue(), sofa::core::objectmodel::ComponentState::Valid);
        ASSERT_EQ(beamForceFieldMass->d_computeMass.getValue(), true);
        ASSERT_FLOAT_EQ(beamForceFieldMass->d_massDensity.getValue(), 10.0);
        ASSERT_FLOAT_EQ(beamForceFieldMass->rayleighMass.getValue(), 0.0);
        ASSERT_FLOAT_EQ(beamForceFieldMass->rayleighStiffness.getValue(), 0.0);
    }



    void checkValues()
    {
        int nbrStep = 500;
        int nbrGrid = 200;
        Node::SPtr root = createSingleBeam();

        // Search for Beam FF
        AdaptiveBeamForceFieldAndMass<Rigid3dTypes>* beamForceFieldMass = nullptr;
        root->getTreeObject(beamForceFieldMass);
        ASSERT_NE(beamForceFieldMass, nullptr);

        // Access mstate
        MechanicalObject<Rigid3dTypes>* dofs = nullptr;
        root->getTreeObject(dofs);

        // Access dofs
        const VecCoord& positions = dofs->x.getValue();
        ASSERT_EQ(positions.size(), nbrGrid);

        // Check position at init
        auto id = nbrGrid - 1;
        EXPECT_NEAR(positions[id][0], 100, 1e-4);
        EXPECT_NEAR(positions[id][1], 0, 1e-4);
        EXPECT_NEAR(positions[id][2], 0, 1e-4);

        // run some simulation steps
        auto simulation = sofa::simulation::getSimulation();
        for (int i = 0; i < nbrStep; i++)
        {
            simulation->animate(root.get(), 0.01);
        }

        // Check position after simulation
        EXPECT_NEAR(positions[id][0], -1.34228, 1e-4);
        EXPECT_NEAR(positions[id][1], -110.274221, 1e-4);
        EXPECT_NEAR(positions[id][2], 0, 1e-4);
    }


    void testPerformances()
    {
        Node::SPtr root = createSingleBeam();

        int nbrStep = 1000;
        int nbrTest = 10;

        double diffTimeMs = 0;
        double timeMin = std::numeric_limits<double>::max();
        double timeMax = std::numeric_limits<double>::min();

        auto simulation = sofa::simulation::getSimulation();
        for (int i = 0; i < nbrTest; ++i)
        {
            ctime_t startTime = sofa::helper::system::thread::CTime::getRefTime();
            for (int i = 0; i < nbrStep; i++)
            {
                simulation->animate(root.get(), 0.01);
            }

            ctime_t diffTime = sofa::helper::system::thread::CTime::getRefTime() - startTime;
            double diffTimed = sofa::helper::system::thread::CTime::toSecond(diffTime);

            if (timeMin > diffTimed)
                timeMin = diffTimed;
            if (timeMax < diffTimed)
                timeMax = diffTimed;

            diffTimeMs += diffTimed;
            simulation->reset(root.get());
        }

        //std::cout << "timeMean: " << diffTimeMs/nbrTest << std::endl;
        //std::cout << "timeMin: " << timeMin << std::endl;
        //std::cout << "timeMax: " << timeMax << std::endl;

        // Some logs (2022-08-08):
        //timeMean: 1.00409
        //timeMin : 0.977432
        //timeMax : 1.02568

    }
};

TEST_F(AdaptiveBeamForceFieldAndMassTest, SimpleScene) {
    ASSERT_NO_THROW(this->simpleSceneTest()) ;
}

TEST_F(AdaptiveBeamForceFieldAndMassTest, checkCreation) {
    ASSERT_NO_THROW(this->checkCreation());
}

TEST_F(AdaptiveBeamForceFieldAndMassTest, checkValues) {
    ASSERT_NO_THROW(this->checkValues());
}

TEST_F(AdaptiveBeamForceFieldAndMassTest, DISABLED_testPerformances) {
    ASSERT_NO_THROW(this->testPerformances());
}

}
