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

#include <BeamAdapter/component/engine/WireRestShape.h>

namespace beamadapter_test
{
using namespace sofa::testing;
using namespace sofa::defaulttype;
using namespace sofa::core::objectmodel;
using namespace sofa::component::engine::_wirerestshape_;

class WireRestShape_test : public BaseTest
{
public:
    void onSetUp() override
    {
        sofa::simpleapi::importPlugin("BeamAdapter");
        sofa::simpleapi::importPlugin("Sofa.Component.Topology.Container.Dynamic");
        //sofa::simpleapi::importPlugin("Sofa.Component.Engine.Select");



            //<RequiredPlugin name = "Sofa.Component.LinearSolver.Direct" / > <!--Needed to use components[BTDLinearSolver] -->
            //<RequiredPlugin name = "Sofa.Component.ODESolver.Backward" / > <!--Needed to use components[EulerImplicitSolver] -->
            //<RequiredPlugin name = "Sofa.Component.SolidMechanics.Spring" / > <!--Needed to use components[RestShapeSpringsForceField] -->
            //<RequiredPlugin name = "Sofa.Component.StateContainer" / > <!--Needed to use components[MechanicalObject] -->
            //<RequiredPlugin name = "Sofa.Component.Topology.Container.Dynamic" / > <!--Needed to use components[EdgeSetGeometryAlgorithms EdgeSetTopologyContainer EdgeSetTopologyModifier] -->
            //<RequiredPlugin name = "Sofa.Component.Topology.Container.Grid" / > <!--Needed to use components[RegularGridTopology] -->
            //<RequiredPlugin name = "Sofa.Component.Visual" / > <!--Needed to use components[VisualStyle] -->

        m_simu = sofa::simpleapi::createSimulation("DAG");
        m_root = sofa::simpleapi::createRootNode(m_simu, "root");
    }

    /// Unload the scene
    void onTearDown() override
    {
        if (m_simu != nullptr && m_root != nullptr) {
            m_simu->unload(m_root);
        }
    }

    void loadScene(const std::string& scene)
    {
        m_root = sofa::simulation::SceneLoaderXML::loadFromMemory("testscene",
            scene.c_str(),
            scene.size());

        EXPECT_NE(m_root.get(), nullptr);

        m_root->init(sofa::core::execparams::defaultInstance());
    }


    /// Test creation of WireRestShape in an empty scene without parameters
    void testEmptyInit();

    /// Test creation of WireRestShape in an empty scene without parameters
    void testDefaultInit();

private:
    /// Pointer to SOFA simulation
    sofa::simulation::Simulation::SPtr m_simu = nullptr;
    /// Pointer to root Node
    sofa::simulation::Node::SPtr m_root = nullptr;

};



void WireRestShape_test::testEmptyInit()
{
    std::string scene =
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>             "
        "   <Node name='BeamTopology'>                                "
        "       <WireRestShape name='BeamRestShape'/>                 "
        "   </Node>                                                   "
        "</Node>                                                      ";

    EXPECT_MSG_EMIT(Error);
    loadScene(scene);

    WireRestShape<Rigid3Types>::SPtr wireRShape = m_root->get< WireRestShape<Rigid3Types> >(sofa::core::objectmodel::BaseContext::SearchDown);
    EXPECT_NE(wireRShape.get(), nullptr);
    EXPECT_EQ(wireRShape->getComponentState(), ComponentState::Invalid);
}


void WireRestShape_test::testDefaultInit()
{
    std::string scene =
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>             "
        "   <Node name='BeamTopology'>                                "
        "       <WireRestShape name='BeamRestShape'/>                 "
        "       <EdgeSetTopologyContainer name='meshLinesBeam'/>      "
        "       <EdgeSetTopologyModifier />                           "
        "   </Node>                                                   "
        "</Node>                                                      ";

    loadScene(scene);


    WireRestShape<Rigid3Types>::SPtr wireRShape = this->m_root->get< WireRestShape<Rigid3Types> >(sofa::core::objectmodel::BaseContext::SearchDown);
    EXPECT_NE(wireRShape.get(), nullptr);
    EXPECT_EQ(wireRShape->getComponentState(), ComponentState::Valid);
}










TEST_F(WireRestShape_test, test_init_empty) {
    testEmptyInit();
}

TEST_F(WireRestShape_test, test_init_default) {
    testDefaultInit();
}

} // namespace beamadapter_test