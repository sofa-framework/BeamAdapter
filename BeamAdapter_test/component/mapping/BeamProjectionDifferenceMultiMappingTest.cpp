/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <sofa/testing/BaseSimulationTest.h>
using sofa::testing::BaseSimulationTest;

#include <SceneCreator/SceneCreator.h>

#include <sofa/simulation/graph/DAGSimulation.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/component/topology/container/constant/MeshTopology.h>
#include <sofa/component/mapping/testing/Multi2MappingTestCreation.h>

#include <BeamAdapter/component/mapping/BeamProjectionDifferenceMultiMapping.h>
#include <BeamAdapter/component/BeamInterpolation.h>

namespace beamadapter::test {

using namespace sofa::core;
using namespace sofa::component;
using sofa::type::Vec;
using sofa::type::Mat;
using sofa::core::objectmodel::New;


/**  Test suite for BeamProjectionDifferenceMultiMapping.
The test cases are defined in the #Test_Cases member group.
  */
template <typename _BeamProjectionDifferenceMultiMapping>
struct BeamProjectionDifferenceMultiMappingTest : public sofa::Multi2Mapping_test<_BeamProjectionDifferenceMultiMapping>
{

    typedef _BeamProjectionDifferenceMultiMapping BeamProjectionDifferenceMultiMapping;
    typedef sofa::Multi2Mapping_test<BeamProjectionDifferenceMultiMapping> Inherit;

    typedef typename BeamProjectionDifferenceMultiMapping::In1 In1DataTypes;
    typedef typename BeamProjectionDifferenceMultiMapping::In2 In2DataTypes;
    typedef statecontainer::MechanicalObject<In1DataTypes> In1MechanicalObject;
    typedef statecontainer::MechanicalObject<In2DataTypes> In2MechanicalObject;

    typedef typename BeamProjectionDifferenceMultiMapping::Out OutDataTypes;
    typedef statecontainer::MechanicalObject<OutDataTypes> OutMechanicalObject;

    typedef typename In1DataTypes::VecCoord In1VecCoord;
    typedef typename In2DataTypes::VecCoord In2VecCoord;
    typedef typename OutDataTypes::VecCoord OutVecCoord;

    BeamProjectionDifferenceMultiMapping* m_mapping;
    sofa::component::fem::BeamInterpolation<sofa::defaulttype::Rigid3Types>* m_interpolation;
    sofa::component::topology::container::constant::MeshTopology* m_topology;


    BeamProjectionDifferenceMultiMappingTest()
    {
        this->setupScene();
        this->errorMax = 10;

        m_mapping = static_cast<BeamProjectionDifferenceMultiMapping*>( this->mapping );
    }


    /** @name Test_Cases
      For each of these cases, we can test if the mapping work
      */
    bool test_oneBeam_twoParticles()
    {
        const int Nin1=1, Nin2=2;
        this->in1Dofs.resize(Nin1);
        this->in2Dofs.resize(Nin2);
        this->outDofs->resize(Nin1);

        // parent position
        In1VecCoord xin1(Nin1);
        In1DataTypes::set( xin1[0], 0.5, 0., 0.);

        // the beam
        In2VecCoord xin2(Nin2);
        In2DataTypes::set( xin2[0], 0., 0., 0. );
        In2DataTypes::set( xin2[1], 1., 0., 0. );

        // expected mapped values
        OutVecCoord expectedChildCoords(Nin1);
        OutDataTypes::set( expectedChildCoords[0], 0., 0., 0.);

        m_topology = sofa::modeling::addNew<sofa::component::topology::container::constant::MeshTopology>(this->parentsIn2).get();
        m_topology->addEdge(0,1);
        m_mapping->l_in2Topology.set(m_topology);
        m_interpolation = sofa::modeling::addNew<sofa::component::fem::BeamInterpolation<sofa::defaulttype::Rigid3Types>>(this->parentsIn2).get();
        m_mapping->l_interpolation.set(m_interpolation);
        m_mapping->d_indices.setValue(sofa::vector<sofa::Index>{0});
        sofa::type::vector<bool> directions{0, 1, 1, 0, 0, 0, 0};
        m_mapping->d_directions.setValue(directions);

        return this->runTest(sofa::vector<In1VecCoord>{xin1}, sofa::vector<In2VecCoord>{xin2}, expectedChildCoords);
    }
};



// Define the list of types to instanciate. We do not necessarily need to test all combinations.
using ::testing::Types;
typedef Types<
beamadapter::mapping::BeamProjectionDifferenceMultiMapping<sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types>
> DataTypes; // the types to instanciate.

// Test suite for all the instanciations
TYPED_TEST_SUITE(BeamProjectionDifferenceMultiMappingTest, DataTypes);

TYPED_TEST( BeamProjectionDifferenceMultiMappingTest, oneBeam_twoParticles )
{
    // child coordinates given directly in parent frame
    ASSERT_TRUE(this->test_oneBeam_twoParticles());
}

} // namespace
