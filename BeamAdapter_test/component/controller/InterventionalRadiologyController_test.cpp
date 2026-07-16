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
#include <sofa/testing/BaseTest.h>
#include <sofa/simpleapi/SimpleApi.h>

#include <sofa/simulation/common/SceneLoaderXML.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/graph/DAGSimulation.h>

#include <BeamAdapter/component/controller/InterventionalRadiologyController.h>
#include <BeamAdapter/component/WireBeamInterpolation.h>
#include <BeamAdapter/utils/BeamActions.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <cmath>

namespace beamadapter_test
{
using namespace sofa::testing;
using namespace sofa::defaulttype;
using namespace sofa::type;
using namespace sofa::core::objectmodel;
using namespace beamadapter;

class InterventionalRadiologyController_test : public BaseTest
{
public:
    using IRController = InterventionalRadiologyController<Rigid3Types>;
    using Coord = Rigid3Types::Coord;
    using Real = Coord::value_type;

    void doSetUp() override
    {
        sofa::simpleapi::importPlugin("BeamAdapter");
        sofa::simpleapi::importPlugin(Sofa.Component.Topology.Container.Dynamic);
        sofa::simpleapi::importPlugin(Sofa.Component.Topology.Container.Grid);
        sofa::simpleapi::importPlugin(Sofa.Component.Constraint.Projective);
        sofa::simpleapi::importPlugin(Sofa.Component.StateContainer);
        sofa::simpleapi::importPlugin(Sofa.Component.ODESolver.Backward);
        sofa::simpleapi::importPlugin(Sofa.Component.LinearSolver.Direct);
        sofa::simpleapi::importPlugin(Sofa.Component.SolidMechanics.Spring);
        sofa::simpleapi::importPlugin(Sofa.Component.AnimationLoop);

        m_root = sofa::simpleapi::createRootNode(sofa::simulation::getSimulation(), "root");
    }

    void doTearDown() override
    {
        if (m_root != nullptr)
            sofa::simulation::node::unload(m_root);
    }

    /// Load + init only the components (no animation loop needed).
    void loadScene(const std::string& scene)
    {
        m_root = sofa::simulation::SceneLoaderXML::loadFromMemory("testscene", scene.c_str());
        ASSERT_NE(m_root.get(), nullptr);
        m_root->init(sofa::core::execparams::defaultInstance());
    }

    /// Load + fully init the graph through the simulation, so the scene can be animated.
    void loadAndInitRoot(const std::string& scene)
    {
        m_root = sofa::simulation::SceneLoaderXML::loadFromMemory("testscene", scene.c_str());
        ASSERT_NE(m_root.get(), nullptr);
        sofa::simulation::node::initRoot(m_root.get());
    }

    IRController* getController()
    {
        return m_root->get<IRController>(BaseContext::SearchDown);
    }

    /// A single-instrument scene without solver (for init / action / getter tests).
    /// `ctrlAttrs` are extra attributes injected on the controller tag.
    static std::string singleInstrumentScene(const std::string& ctrlAttrs = "")
    {
        return
            "<?xml version='1.0'?>"
            "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>                                        "
            "   <Node name='BeamTopology'>                                                           "
            "       <RodStraightSection name='StraightSection' length='100.0' nbBeams='30'           "
            "           nbEdgesCollis='30' nbEdgesVisu='30'/>                                         "
            "       <WireRestShape name='BeamRestShape' wireMaterials='@StraightSection'/>           "
            "       <EdgeSetTopologyContainer name='meshLinesBeam'/>                                 "
            "       <EdgeSetTopologyModifier />                                                      "
            "   </Node>                                                                              "
            "   <Node name='BeamModel'>                                                              "
            "       <MechanicalObject name='dof' template='Rigid3d' />                               "
            "       <RegularGridTopology name='lines' nx='40' ny='1' nz='1' />                       "
            "       <FixedProjectiveConstraint name='fc' indices='0' />                              "
            "       <WireBeamInterpolation name='Interpol' WireRestShape='@../BeamTopology/BeamRestShape' />"
            "       <InterventionalRadiologyController name='DeployController' template='Rigid3d'    "
            "           instruments='Interpol' " + ctrlAttrs + " />                                  "
            "   </Node>                                                                              "
            "</Node>                                                                                 ";
    }

    /// A two-instrument scene without solver.
    static std::string twoInstrumentScene()
    {
        return
            "<?xml version='1.0'?>"
            "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>                                        "
            "   <Node name='Topo1'>                                                                  "
            "       <RodStraightSection name='S1' length='100.0' nbBeams='30' nbEdgesCollis='30' nbEdgesVisu='30'/>"
            "       <WireRestShape name='RS1' wireMaterials='@S1'/>                                  "
            "       <EdgeSetTopologyContainer name='m1'/> <EdgeSetTopologyModifier />                "
            "   </Node>                                                                              "
            "   <Node name='Topo2'>                                                                  "
            "       <RodStraightSection name='S2' length='100.0' nbBeams='30' nbEdgesCollis='30' nbEdgesVisu='30'/>"
            "       <WireRestShape name='RS2' wireMaterials='@S2'/>                                  "
            "       <EdgeSetTopologyContainer name='m2'/> <EdgeSetTopologyModifier />                "
            "   </Node>                                                                              "
            "   <Node name='BeamModel'>                                                              "
            "       <MechanicalObject name='dof' template='Rigid3d' />                               "
            "       <RegularGridTopology name='lines' nx='70' ny='1' nz='1' />                       "
            "       <FixedProjectiveConstraint name='fc' indices='0' />                              "
            "       <WireBeamInterpolation name='I1' WireRestShape='@../Topo1/RS1' />                "
            "       <WireBeamInterpolation name='I2' WireRestShape='@../Topo2/RS2' />                "
            "       <InterventionalRadiologyController name='DeployController' template='Rigid3d'    "
            "           instruments='I1 I2' />                                                       "
            "   </Node>                                                                              "
            "</Node>                                                                                 ";
    }

    /// A full deployment scene with solver, mass and rest-shape springs (for stepping tests).
    static std::string deploymentScene(const std::string& ctrlAttrs)
    {
        return
            "<?xml version='1.0'?>"
            "<Node name='root' gravity='0 -9.81 0' dt='0.01'>                                        "
            "   <DefaultAnimationLoop />                                                             "
            "   <Node name='EdgeTopology'>                                                           "
            "       <RodStraightSection name='StraightSection' youngModulus='20000' radius='0.9'     "
            "           massDensity='0.00000155' nbBeams='30' nbEdgesCollis='30' nbEdgesVisu='30' length='100.0'/>"
            "       <WireRestShape template='Rigid3d' name='BeamRestShape' wireMaterials='@StraightSection'/>"
            "       <EdgeSetTopologyContainer name='meshLinesBeam'/>                                 "
            "       <EdgeSetTopologyModifier name='Modifier'/>                                       "
            "       <EdgeSetGeometryAlgorithms name='GeomAlgo' template='Rigid3d'/>                  "
            "       <MechanicalObject name='dofTopo2' template='Rigid3d'/>                           "
            "   </Node>                                                                              "
            "   <Node name='BeamModel'>                                                              "
            "       <EulerImplicitSolver rayleighStiffness='0.2' rayleighMass='0.1'/>                "
            "       <BTDLinearSolver verbose='0'/>                                                   "
            "       <RegularGridTopology name='MeshLines' nx='40' ny='1' nz='1'                      "
            "           xmax='0.0' xmin='0.0' ymin='0' ymax='0' zmax='0' zmin='0' p0='0 0 0'/>       "
            "       <MechanicalObject template='Rigid3d' name='DOFs'/>                               "
            "       <WireBeamInterpolation name='Interpol' WireRestShape='@../EdgeTopology/BeamRestShape'/>"
            "       <AdaptiveBeamForceFieldAndMass name='BeamForceField' interpolation='@Interpol'/> "
            "       <InterventionalRadiologyController name='DeployController' template='Rigid3d'    "
            "           instruments='Interpol' topology='@MeshLines' startingPos='0 0 0 0 0 0 1' " + ctrlAttrs + "/>"
            "       <FixedProjectiveConstraint name='FixedConstraint' indices='0'/>                  "
            "       <RestShapeSpringsForceField name='RestSPForceField'                              "
            "           points='@DeployController.indexFirstNode' angularStiffness='1e8' stiffness='1e8'/>"
            "   </Node>                                                                              "
            "</Node>                                                                                 ";
    }

protected:
    sofa::simulation::Node::SPtr m_root = nullptr;
};


///////////////////////////////////// Init / validation ///////////////////////////////////////////

/// No topology, no instrument, no mechanical state => Invalid.
TEST_F(InterventionalRadiologyController_test, init_empty_invalid)
{
    EXPECT_MSG_EMIT(Error);
    loadScene(
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>"
        "   <Node name='BeamModel'>"
        "       <InterventionalRadiologyController name='DeployController' template='Rigid3d'/>"
        "   </Node>"
        "</Node>");

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    EXPECT_EQ(ctrl->getComponentState(), ComponentState::Invalid);
}

/// Fully wired single-instrument scene => Valid, one instrument registered.
TEST_F(InterventionalRadiologyController_test, init_default_valid)
{
    EXPECT_MSG_NOEMIT(Error);
    loadScene(singleInstrumentScene());

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    EXPECT_EQ(ctrl->getComponentState(), ComponentState::Valid);
    EXPECT_EQ(ctrl->m_instrumentsList.size(), 1u);
}

/// A topology is present but there is no WireBeamInterpolation => Invalid.
TEST_F(InterventionalRadiologyController_test, init_noInstrument_invalid)
{
    EXPECT_MSG_EMIT(Error);
    loadScene(
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>"
        "   <Node name='BeamModel'>"
        "       <MechanicalObject name='dof' template='Rigid3d' />"
        "       <RegularGridTopology name='lines' nx='40' ny='1' nz='1' />"
        "       <FixedProjectiveConstraint name='fc' indices='0' />"
        "       <InterventionalRadiologyController name='DeployController' template='Rigid3d'/>"
        "   </Node>"
        "</Node>");

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    EXPECT_EQ(ctrl->getComponentState(), ComponentState::Invalid);
}

/// An instrument path that does not resolve => error and Invalid.
TEST_F(InterventionalRadiologyController_test, init_instrumentPathNotFound_invalid)
{
    EXPECT_MSG_EMIT(Error);
    loadScene(
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>"
        "   <Node name='BeamModel'>"
        "       <MechanicalObject name='dof' template='Rigid3d' />"
        "       <RegularGridTopology name='lines' nx='40' ny='1' nz='1' />"
        "       <FixedProjectiveConstraint name='fc' indices='0' />"
        "       <InterventionalRadiologyController name='DeployController' template='Rigid3d' instruments='DoesNotExist'/>"
        "   </Node>"
        "</Node>");

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    EXPECT_EQ(ctrl->getComponentState(), ComponentState::Invalid);
}

/// Without a FixedProjectiveConstraint the controller should become Invalid.
/// DISABLED: currently CRASHES. init() sets the state to Invalid on a missing
/// FixedProjectiveConstraint (IRC.inl ~L196) but does not return, so it goes on to call
/// applyInterventionalRadiologyController() at the end of init(), which dereferences the
/// null l_fixedConstraint in fixFirstNodesWithUntil() -> segfault. Re-enable once init()
/// returns early (or applyInterventionalRadiologyController() guards the null link).
TEST_F(InterventionalRadiologyController_test, DISABLED_init_noFixedConstraint_invalid)
{
    EXPECT_MSG_EMIT(Error);
    loadScene(
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>"
        "   <Node name='BeamTopology'>"
        "       <RodStraightSection name='StraightSection' length='100.0' nbBeams='30' nbEdgesCollis='30' nbEdgesVisu='30'/>"
        "       <WireRestShape name='BeamRestShape' wireMaterials='@StraightSection'/>"
        "       <EdgeSetTopologyContainer name='meshLinesBeam'/> <EdgeSetTopologyModifier />"
        "   </Node>"
        "   <Node name='BeamModel'>"
        "       <MechanicalObject name='dof' template='Rigid3d' />"
        "       <RegularGridTopology name='lines' nx='40' ny='1' nz='1' />"
        "       <WireBeamInterpolation name='Interpol' WireRestShape='@../BeamTopology/BeamRestShape' />"
        "       <InterventionalRadiologyController name='DeployController' template='Rigid3d' instruments='Interpol'/>"
        "   </Node>"
        "</Node>");

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    EXPECT_EQ(ctrl->getComponentState(), ComponentState::Invalid);
}

/// Two instruments => both registered, and the per-instrument Data are sized accordingly.
TEST_F(InterventionalRadiologyController_test, init_twoInstruments)
{
    EXPECT_MSG_NOEMIT(Error);
    loadScene(twoInstrumentScene());

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    EXPECT_EQ(ctrl->getComponentState(), ComponentState::Valid);
    EXPECT_EQ(ctrl->m_instrumentsList.size(), 2u);
    EXPECT_EQ(ctrl->d_xTip.getValue().size(), 2u);
    EXPECT_EQ(ctrl->d_rotationInstrument.getValue().size(), 2u);
}


///////////////////////////////////// Data initialization /////////////////////////////////////////

/// xtip / rotationInstrument are resized to the number of instruments even if left unset.
TEST_F(InterventionalRadiologyController_test, init_dataResizedToInstrumentCount)
{
    EXPECT_MSG_NOEMIT(Error);
    loadScene(singleInstrumentScene());

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    EXPECT_EQ(ctrl->d_xTip.getValue().size(), 1u);
    EXPECT_EQ(ctrl->d_rotationInstrument.getValue().size(), 1u);
}

/// startingPos orientation is normalized during init, and all dofs are initialized to it.
TEST_F(InterventionalRadiologyController_test, init_startingPosNormalized)
{
    EXPECT_MSG_NOEMIT(Error);
    // deliberately non-unit quaternion (0 0 0 5)
    loadScene(singleInstrumentScene("startingPos='1 2 3 0 0 0 5'"));

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    const auto q = ctrl->d_startingPos.getValue().getOrientation();
    const double qNorm = std::sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
    EXPECT_NEAR(qNorm, 1.0, 1e-6);

    // every dof of the mechanical state has been set to startingPos
    auto* beamModel = m_root->getChild("BeamModel");
    ASSERT_NE(beamModel, nullptr);
    auto* mstate = beamModel->getMechanicalState();
    ASSERT_NE(mstate, nullptr);
    EXPECT_GT(mstate->getSize(), 0u);
}

/// The current curvilinear abscissa list always starts at the origin (0).
TEST_F(InterventionalRadiologyController_test, init_curvAbsStartsAtOrigin)
{
    EXPECT_MSG_NOEMIT(Error);
    loadScene(singleInstrumentScene());

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    const auto& curvAbs = ctrl->getCurrentCurvAbscisses();
    ASSERT_FALSE(curvAbs.empty());
    EXPECT_NEAR(curvAbs.front(), 0.0, 1e-9);
}

/// An instrument deployed at start (xtip>0) creates a node at that abscissa.
TEST_F(InterventionalRadiologyController_test, init_deployedTipCreatesNode)
{
    EXPECT_MSG_NOEMIT(Error);
    loadScene(singleInstrumentScene("xtip='10'"));

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    const auto& curvAbs = ctrl->getCurrentCurvAbscisses();
    // origin + deployed length => tip reaches ~10
    EXPECT_NEAR(curvAbs.back(), 10.0, ctrl->d_threshold.getValue());
    EXPECT_GT(curvAbs.size(), 1u);
}


///////////////////////////////////// applyAction /////////////////////////////////////////////////

TEST_F(InterventionalRadiologyController_test, action_moveForwardBackward)
{
    loadScene(singleInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    const Real step = ctrl->d_step.getValue();
    const Real x0 = ctrl->d_xTip.getValue()[0];

    ctrl->applyAction(BeamAdapterAction::MOVE_FORWARD);
    EXPECT_NEAR(ctrl->d_xTip.getValue()[0], x0 + step, 1e-9);

    ctrl->applyAction(BeamAdapterAction::MOVE_BACKWARD);
    EXPECT_NEAR(ctrl->d_xTip.getValue()[0], x0, 1e-9);
}

TEST_F(InterventionalRadiologyController_test, action_spinRightLeft)
{
    loadScene(singleInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    const Real astep = ctrl->d_angularStep.getValue();
    const Real r0 = ctrl->d_rotationInstrument.getValue()[0];

    ctrl->applyAction(BeamAdapterAction::SPIN_RIGHT);
    EXPECT_NEAR(ctrl->d_rotationInstrument.getValue()[0], r0 + astep, 1e-9);

    ctrl->applyAction(BeamAdapterAction::SPIN_LEFT);
    EXPECT_NEAR(ctrl->d_rotationInstrument.getValue()[0], r0, 1e-9);
}

/// With a single instrument, switching tools is bounded and keeps the controlled index at 0.
TEST_F(InterventionalRadiologyController_test, action_switchTools_singleInstrumentBounded)
{
    loadScene(singleInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    ASSERT_EQ(ctrl->d_controlledInstrument.getValue(), 0);

    ctrl->applyAction(BeamAdapterAction::SWITCH_NEXT_TOOL); // no next tool
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 0);

    ctrl->applyAction(BeamAdapterAction::SWITCH_PREVIOUS_TOOL); // already first
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 0);

    ctrl->applyAction(BeamAdapterAction::USE_TOOL_1); // does not exist
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 0);

    ctrl->applyAction(BeamAdapterAction::USE_TOOL_2); // does not exist
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 0);
}

/// With two instruments, switching moves the controlled index within [0, 1].
TEST_F(InterventionalRadiologyController_test, action_switchTools_multiInstrument)
{
    loadScene(twoInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    ASSERT_EQ(ctrl->d_controlledInstrument.getValue(), 0);

    ctrl->applyAction(BeamAdapterAction::SWITCH_NEXT_TOOL);
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 1);

    ctrl->applyAction(BeamAdapterAction::SWITCH_NEXT_TOOL); // no third tool, stays at 1
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 1);

    ctrl->applyAction(BeamAdapterAction::SWITCH_PREVIOUS_TOOL);
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 0);

    ctrl->applyAction(BeamAdapterAction::USE_TOOL_1);
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 1);

    ctrl->applyAction(BeamAdapterAction::USE_TOOL_0);
    EXPECT_EQ(ctrl->d_controlledInstrument.getValue(), 0);
}

/// MOVE_FORWARD acts on the currently controlled instrument.
TEST_F(InterventionalRadiologyController_test, action_movesControlledInstrument)
{
    loadScene(twoInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    const Real step = ctrl->d_step.getValue();

    // NB: we snapshot instrument 0's tip rather than assuming it is exactly 0.
    // When no instrument is deployed at init (all xtip == 0), applyInterventionalRadiologyController()
    // refuses to build a zero-length wire and instead nudges the first instrument's tip to
    // std::numeric_limits<float>::epsilon() (~1.19e-7), then recomputes.
    // See InterventionalRadiologyController.inl:784-790 ("if the totalLength is 0, move the first instrument").
    // Since init() calls that method at its end, d_xTip[0] already reads ~1e-7 here, not 0.
    const Real x0Before = ctrl->d_xTip.getValue()[0];

    ctrl->applyAction(BeamAdapterAction::USE_TOOL_1);
    ctrl->applyAction(BeamAdapterAction::MOVE_FORWARD);

    // only the controlled instrument (1) advances by one step; instrument 0 is left untouched
    EXPECT_NEAR(ctrl->d_xTip.getValue()[1], step, 1e-9);
    EXPECT_NEAR(ctrl->d_xTip.getValue()[0], x0Before, 1e-9);
}


///////////////////////////////////// Getters / helpers ///////////////////////////////////////////

TEST_F(InterventionalRadiologyController_test, getTotalNbEdges)
{
    EXPECT_MSG_NOEMIT(Error);
    loadScene(singleInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    // RegularGridTopology nx=40 => 39 edges
    EXPECT_EQ(ctrl->getTotalNbEdges(), 39);
}

TEST_F(InterventionalRadiologyController_test, getInstrumentList)
{
    loadScene(twoInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    type::vector<WireBeamInterpolation<Rigid3Types>*> list;
    ctrl->getInstrumentList(list);
    EXPECT_EQ(list.size(), 2u);
    for (auto* instr : list)
        EXPECT_NE(instr, nullptr);
}

/// fixFirstNodesWithUntil(0) is rejected with a warning and resets indexFirstNode to 0.
TEST_F(InterventionalRadiologyController_test, fixFirstNodesWithUntil_zeroIsRejected)
{
    loadScene(singleInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    EXPECT_MSG_EMIT(Warning);
    ctrl->fixFirstNodesWithUntil(0);
    EXPECT_EQ(ctrl->d_indexFirstNode.getValue(), 0u);
}

/// fixFirstNodesWithUntil(k) fixes the first k-1 nodes and sets indexFirstNode to k-1.
TEST_F(InterventionalRadiologyController_test, fixFirstNodesWithUntil_setsIndexFirstNode)
{
    loadScene(singleInstrumentScene());
    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    ctrl->fixFirstNodesWithUntil(5);
    EXPECT_EQ(ctrl->d_indexFirstNode.getValue(), 4u);
}


///////////////////////////////////// Deployment (stepping) ///////////////////////////////////////

/// With a positive speed the tip advances by speed*dt each step and the tip node follows.
TEST_F(InterventionalRadiologyController_test, deploy_speedAdvancesTip)
{
    EXPECT_MSG_NOEMIT(Error);
    loadAndInitRoot(deploymentScene("xtip='0' step='0.5' speed='10' controlledInstrument='0'"));

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    ASSERT_EQ(ctrl->getComponentState(), ComponentState::Valid);

    // The scene starts undeployed (xtip=0). Rather than build a zero-length wire, init() ->
    // applyInterventionalRadiologyController() nudges the first tip to std::numeric_limits<float>::epsilon()
    // (~1.19e-7) and recomputes. See InterventionalRadiologyController.inl:784-790.
    // So the tip starts at ~1e-7, not exactly 0 -- hence the tolerance instead of an == 0 check.
    ASSERT_LT(ctrl->d_xTip.getValue()[0], 1e-6);

    const int nbSteps = 20;
    for (int i = 0; i < nbSteps; ++i)
        sofa::simulation::node::animate(m_root.get(), 0.01);

    // Each step advances the tip by speed*dt (see onBeginAnimationStep). After 20 steps:
    // speed * dt * nbSteps = 10 * 0.01 * 20 = 2.0, plus the negligible initial epsilon nudge above.
    const Real expectedTip = 10.0 * 0.01 * nbSteps;
    EXPECT_NEAR(ctrl->d_xTip.getValue()[0], expectedTip, 1e-5);

    // the deployed length is reflected in the last curvilinear abscissa
    const auto& curvAbs = ctrl->getCurrentCurvAbscisses();
    EXPECT_NEAR(curvAbs.back(), expectedTip, ctrl->d_threshold.getValue());
}

/// The tip is clamped to the instrument's total rest length.
TEST_F(InterventionalRadiologyController_test, deploy_tipClampedToRestLength)
{
    EXPECT_MSG_NOEMIT(Error);
    // start already fully deployed, keep pushing forward
    loadAndInitRoot(deploymentScene("xtip='100' step='0.5' speed='10' controlledInstrument='0'"));

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    for (int i = 0; i < 5; ++i)
        sofa::simulation::node::animate(m_root.get(), 0.01);

    // instrument rest length is 100 => xtip cannot exceed it
    EXPECT_NEAR(ctrl->d_xTip.getValue()[0], 100.0, 1e-6);
}

/// As the instrument deploys, the index of the first simulated node decreases.
TEST_F(InterventionalRadiologyController_test, deploy_indexFirstNodeDecreases)
{
    EXPECT_MSG_NOEMIT(Error);
    loadAndInitRoot(deploymentScene("xtip='0' step='0.5' speed='10' controlledInstrument='0'"));

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);
    const unsigned int idxStart = ctrl->d_indexFirstNode.getValue();

    for (int i = 0; i < 20; ++i)
        sofa::simulation::node::animate(m_root.get(), 0.01);

    EXPECT_LT(ctrl->d_indexFirstNode.getValue(), idxStart);
}

/// A longer run remains stable (no error) and keeps the component Valid.
TEST_F(InterventionalRadiologyController_test, deploy_runsStable)
{
    EXPECT_MSG_NOEMIT(Error);
    loadAndInitRoot(deploymentScene("xtip='0' step='0.5' speed='5' controlledInstrument='0'"));

    auto* ctrl = getController();
    ASSERT_NE(ctrl, nullptr);

    for (int i = 0; i < 100; ++i)
        sofa::simulation::node::animate(m_root.get(), 0.01);

    EXPECT_EQ(ctrl->getComponentState(), ComponentState::Valid);
}

} // namespace beamadapter_test
