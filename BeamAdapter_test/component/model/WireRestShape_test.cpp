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
#include <sofa/core/topology/BaseMeshTopology.h>

namespace beamadapter_test
{
using namespace sofa::testing;
using namespace sofa::defaulttype;
using namespace sofa::type;
using namespace sofa::core::objectmodel;
using namespace sofa::component::engine::_wirerestshape_;

class WireRestShape_test : public BaseTest
{
public:
    using WireRestShapeRig3 = WireRestShape<Rigid3Types>;
    typedef typename Rigid3Types::Coord    Coord;
    typedef typename Rigid3Types::Deriv    Deriv;
    typedef typename Coord::value_type   Real;
    typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;

    void onSetUp() override
    {
        m_simu = sofa::simpleapi::createSimulation("DAG");
        sofa::simulation::setSimulation(m_simu.get());
        m_root = sofa::simpleapi::createRootNode(m_simu, "root");

        sofa::simpleapi::importPlugin("BeamAdapter");
        sofa::simpleapi::importPlugin("Sofa.Component.Topology.Container.Dynamic");
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
        m_root = sofa::simulation::SceneLoaderXML::loadFromMemory("testscene", scene.c_str());

        EXPECT_NE(m_root.get(), nullptr);

        m_root->init(sofa::core::execparams::defaultInstance());
    }


    /// Test creation of WireRestShape in an empty scene without parameters
    void testEmptyInit();

    /// Test creation of WireRestShape in a default scene with needed Topology components
    void testDefaultInit();

    /// Test creation of WireRestShape in a default scene and check parameters cohesion 
    void testParameterInit();

    /// Test creation of WireRestShape in a default scene and check created topology 
    void testTopologyInit();

    // TODO test on topology init from MeshLoader

    /// Test WireRestShape transform methods 
    void testTransformMethods();

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

    WireRestShapeRig3::SPtr wireRShape = m_root->get< WireRestShapeRig3 >(sofa::core::objectmodel::BaseContext::SearchDown);
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


    WireRestShapeRig3::SPtr wireRShape = this->m_root->get< WireRestShapeRig3 >(sofa::core::objectmodel::BaseContext::SearchDown);
    EXPECT_NE(wireRShape.get(), nullptr);
    //TODO: change component to match behavior expected by tests
    //EXPECT_EQ(wireRShape->getComponentState(), ComponentState::Valid);
}


void WireRestShape_test::testParameterInit()
{
    std::string scene =
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>             "
        "   <Node name='BeamTopology'>                                "
        "       <WireRestShape name='BeamRestShape' length='100.0'    "
        "        straightLength='95.0'/>                              "
        "       <EdgeSetTopologyContainer name='meshLinesBeam'/>      "
        "       <EdgeSetTopologyModifier />                           "
        "   </Node>                                                   "
        "</Node>                                                      ";

    loadScene(scene);

    WireRestShapeRig3::SPtr wireRShape = this->m_root->get< WireRestShapeRig3 >(sofa::core::objectmodel::BaseContext::SearchDown);
    WireRestShapeRig3* wire = wireRShape.get();
    EXPECT_NE(wire, nullptr);

    Real fullLength = wire->getLength();
    EXPECT_EQ(fullLength, 100.0);

    Real straightLength = wire->getReleaseCurvAbs();
    EXPECT_EQ(straightLength, 95.0);
    
    vector<Real> keysPoints, keysPoints_ref = { 0, straightLength, fullLength };
    Real ratio = straightLength / fullLength;
    vector<int> nbP_density, nbP_density_ref = { int(floor(5.0 * ratio)), int(floor(20.0 * (1 - ratio))) };
    
    wire->getSamplingParameters(keysPoints, nbP_density);
    EXPECT_EQ(keysPoints.size(), 3);
    EXPECT_EQ(keysPoints, keysPoints_ref);
   
    EXPECT_EQ(nbP_density.size(), 2);
    EXPECT_EQ(nbP_density, nbP_density_ref);

    Real dx1, dx2, dx3, nbEdgesCol_ref = 20;
    wire->getCollisionSampling(dx1, 0.0);
    wire->getCollisionSampling(dx2, fullLength);
    wire->getCollisionSampling(dx3, 90.0);
    EXPECT_EQ(dx1, straightLength / nbEdgesCol_ref);
    EXPECT_EQ(dx2, (fullLength - straightLength) / nbEdgesCol_ref);
    EXPECT_EQ(dx3, straightLength / nbEdgesCol_ref);
}


void WireRestShape_test::testTopologyInit()
{
    std::string scene =
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>             "
        "   <Node name='BeamTopology'>                                "
        "       <WireRestShape name='BeamRestShape' length='100.0'    "
        "        straightLength='95.0' numEdges='30'/>        "
        "       <EdgeSetTopologyContainer name='meshLinesBeam'/>      "
        "       <EdgeSetTopologyModifier />                           "
        "   </Node>                                                   "
        "</Node>                                                      ";

    loadScene(scene);

    // get node of the mesh
    sofa::simulation::Node* beam = m_root->getChild("BeamTopology");
    EXPECT_NE(beam, nullptr);

    // getting topology
    sofa::core::topology::BaseMeshTopology* topo = beam->getMeshTopology();
    EXPECT_NE(topo, nullptr);

    // checking topo created by WireRestShape
    int numbEdgesVisu = 30;
    EXPECT_EQ(topo->getNbPoints(), numbEdgesVisu + 1);
    EXPECT_EQ(topo->getNbEdges(), numbEdgesVisu);

    Real dx = 100.0 / Real(numbEdgesVisu);

    EXPECT_EQ(topo->getPX(0), 0.0);
    EXPECT_EQ(topo->getPX(1), dx);
    EXPECT_EQ(topo->getPX(10), 10*dx);
    EXPECT_EQ(topo->getPY(10), 0.0);

    EXPECT_EQ(topo->getEdge(10)[0], 10);
    EXPECT_EQ(topo->getEdge(10)[1], 11);
}


void WireRestShape_test::testTransformMethods()
{
    std::string scene =
        "<?xml version='1.0'?>"
        "<Node name='Root' gravity='0 -9.81 0' dt='0.01'>             "
        "   <Node name='BeamTopology'>                                "
        "       <WireRestShape name='BeamRestShape' length='100.0'    "
        "        straightLength='95.0' numEdges='30'/>        "
        "       <EdgeSetTopologyContainer name='meshLinesBeam'/>      "
        "       <EdgeSetTopologyModifier />                           "
        "   </Node>                                                   "
        "</Node>                                                      ";

    loadScene(scene);

    WireRestShapeRig3::SPtr wireRShape = this->m_root->get< WireRestShapeRig3 >(sofa::core::objectmodel::BaseContext::SearchDown);
    WireRestShapeRig3* wire = wireRShape.get();
    EXPECT_NE(wire, nullptr);

    Real fullLength = wire->getLength();
    Real straightLength = wire->getReleaseCurvAbs();
    Real middHook = (fullLength + straightLength) / 2;

    Transform transfo_0, transfo_1, transfo_2, transfo_3;
    Transform transfo_0_ref, transfo_1_ref, transfo_2_ref, transfo_3_ref;
    transfo_0_ref.set(Vec3(0.0, 0.0, 0.0), Quat<Real>());
    transfo_1_ref.set(Vec3(straightLength / 2, 0.0, 0.0), Quat<Real>());
    transfo_2_ref.set(Vec3(94.972534, 0.15907438, -0.0082204072), Quat<Real>(-0.00456144, -0.28667623, 0.01524152, 0.9578954));
    transfo_3_ref.set(Vec3(94.9857, 0.07953719, -0.0020991669), Quat<Real>(-0.0023050411, -0.14486668, 0.015741583, 0.989323));

    wire->getRestTransformOnX(transfo_0, 0.0);
    wire->getRestTransformOnX(transfo_1, straightLength / 2);
    wire->getRestTransformOnX(transfo_2, fullLength);
    wire->getRestTransformOnX(transfo_3, middHook);
    
    for (int i = 0; i < 3; i++)
    {
        EXPECT_NEAR(transfo_0.getOrigin()[i], transfo_0_ref.getOrigin()[i], 0.0001);
        EXPECT_NEAR(transfo_1.getOrigin()[i], transfo_1_ref.getOrigin()[i], 0.0001);
        EXPECT_NEAR(transfo_2.getOrigin()[i], transfo_2_ref.getOrigin()[i], 0.0001);
        EXPECT_NEAR(transfo_3.getOrigin()[i], transfo_3_ref.getOrigin()[i], 0.0001);
    }

    for (int i = 0; i < 4; i++)
    {
        EXPECT_NEAR(transfo_0.getOrientation()[i], transfo_0_ref.getOrientation()[i], 0.0001);
        EXPECT_NEAR(transfo_1.getOrientation()[i], transfo_1_ref.getOrientation()[i], 0.0001);
        EXPECT_NEAR(transfo_2.getOrientation()[i], transfo_2_ref.getOrientation()[i], 0.0001);
        EXPECT_NEAR(transfo_3.getOrientation()[i], transfo_3_ref.getOrientation()[i], 0.0001);
    }

}



TEST_F(WireRestShape_test, test_init_empty) {
    testEmptyInit();
}

TEST_F(WireRestShape_test, test_init_default) {
    testDefaultInit();
}

TEST_F(WireRestShape_test, test_init_parameters) {
    testParameterInit();
}

TEST_F(WireRestShape_test, test_init_topology) {
    testTopologyInit();
}

TEST_F(WireRestShape_test, test_Transform_methods) {
    testTransformMethods();
}

} // namespace beamadapter_test
