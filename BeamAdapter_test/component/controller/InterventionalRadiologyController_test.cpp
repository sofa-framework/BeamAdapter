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
#include <sofa/simulation/graph/SimpleApi.h>

#include <sofa/simulation/common/SceneLoaderXML.h>
#include <sofa/simulation/Node.h>
#include <sofa/simulation/graph/DAGSimulation.h>

#include <BeamAdapter/component/controller/InterventionalRadiologyController.h>
#include <sofa/core/topology/BaseMeshTopology.h>

namespace beamadapter_test
{
using namespace sofa::testing;
using namespace sofa::defaulttype;
using namespace sofa::type;
using namespace sofa::core::objectmodel;
using namespace sofa::component::controller::_interventionalradiologycontroller_;

class InterventionalRadiologyController_test : public BaseTest
{
public:
    using InterventionalRadioCtrlRig3 = InterventionalRadiologyController<Rigid3Types>;
    typedef typename Rigid3Types::Coord    Coord;
    typedef typename Rigid3Types::Deriv    Deriv;
    typedef typename Coord::value_type   Real;
    //typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;

    void onSetUp() override
    {
        sofa::simpleapi::importPlugin("BeamAdapter");
        sofa::simpleapi::importPlugin("Sofa.Component.Topology.Container.Dynamic");
        sofa::simpleapi::importPlugin("Sofa.Component.Topology.Container.Grid");
        sofa::simpleapi::importPlugin("Sofa.Component.Constraint.Projective");

        m_root = sofa::simpleapi::createRootNode(sofa::simulation::getSimulation(), "root");
    }

    /// Unload the scene
    void onTearDown() override
    {
        if (m_root != nullptr) {
            sofa::simulation::node::unload(m_root);
        }
    }

    void loadScene(const std::string& scene)
    {
        m_root = sofa::simulation::SceneLoaderXML::loadFromMemory("testscene", scene.c_str());

        EXPECT_NE(m_root.get(), nullptr);

        m_root->init(sofa::core::execparams::defaultInstance());
    }


    /// Test creation of InterventionalRadiologyController in an empty scene without parameters
    void testEmptyInit();

    /// Test creation of InterventionalRadiologyController in a default scene with needed components
    void testDefaultInit();

    ///// Test creation of WireRestShape in a default scene and check parameters cohesion 
    void testParameterInit();

    /** TODO epernod 2023-05-03: 
     - Test creation of WireRestShape in a default scene and check created topology 
     - Test WireRestShape topology init from MeshLoader
     - Test WireRestShape transform methods 
     */   

private:
    /// Pointer to root Node
    sofa::simulation::Node::SPtr m_root = nullptr;

};



void InterventionalRadiologyController_test::testEmptyInit()
{
    std::string scene =
        "<?xml version='1.0'?>                                                                  "
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>                                       "
        "   <Node name='BeamModel'>                                                             "
        "       <InterventionalRadiologyController name='DeployController' template='Rigid3d'/> "
        "   </Node>                                                                             "
        "</Node>                                                                                ";

    EXPECT_MSG_EMIT(Error);
    loadScene(scene);

    InterventionalRadioCtrlRig3::SPtr IRCtrlPtr = m_root->get< InterventionalRadioCtrlRig3 >(sofa::core::objectmodel::BaseContext::SearchDown);
    EXPECT_NE(IRCtrlPtr.get(), nullptr);
    EXPECT_EQ(IRCtrlPtr->getComponentState(), ComponentState::Invalid);
}


void InterventionalRadiologyController_test::testDefaultInit()
{
    std::string scene =
        "<?xml version='1.0'?>                                                                  "
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>                                       "
        "   <Node name='BeamTopology'>                                                          "
        "       <RodStraightSection name='StraightSection'/>                                    "
        "       <WireRestShape name='BeamRestShape' wireMaterials='@StraightSection'/>          "
        "       <EdgeSetTopologyContainer name='meshLinesBeam'/>                                "
        "       <EdgeSetTopologyModifier />                                                     "
        "   </Node>                                                                             "
        "   <Node name='BeamModel'>                                                             "
        "       <MechanicalObject name='dof' template='Rigid3d' />                              "
        "       <RegularGridTopology name='lines' nx='60' ny='1' nz='1' />                      "
        "       <FixedConstraint name='fc' indices='0' />                                       "
        "       <WireBeamInterpolation name='BeamInterpolation' WireRestShape='@../BeamTopology/BeamRestShape' />"
        "       <InterventionalRadiologyController name='DeployController' template='Rigid3d' instruments='BeamInterpolation' /> "
        "   </Node>                                                                             "
        "</Node>                                                      ";

    EXPECT_MSG_NOEMIT(Error);
    loadScene(scene);

    InterventionalRadioCtrlRig3::SPtr IRCtrlPtr = m_root->get< InterventionalRadioCtrlRig3 >(sofa::core::objectmodel::BaseContext::SearchDown);
    EXPECT_NE(IRCtrlPtr.get(), nullptr);
    EXPECT_EQ(IRCtrlPtr->getComponentState(), ComponentState::Valid);

    EXPECT_EQ(IRCtrlPtr.get()->m_instrumentsList.size(), 1);
}



TEST_F(InterventionalRadiologyController_test, test_init_empty) {
    testEmptyInit();
}

TEST_F(InterventionalRadiologyController_test, test_init_default) {
    testDefaultInit();
}

} // namespace beamadapter_test
