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

#include <string>
using std::string;

using sofa::core::topology::BaseMeshTopology ;
using sofa::core::objectmodel::Data ;

using sofa::helper::WriteAccessor ;
using sofa::defaulttype::Rigid3dTypes ;

#include <sofa/simulation/common/SceneLoaderXML.h>
using sofa::simulation::SceneLoaderXML ;

#include <sofa/simulation/graph/DAGSimulation.h>
using sofa::simulation::graph::DAGSimulation;
using sofa::simulation::Simulation ;
using sofa::simulation::Node ;
using sofa::simulation::setSimulation ;
using sofa::core::objectmodel::New ;
using sofa::core::objectmodel::BaseData ;
using sofa::component::statecontainer::MechanicalObject ;

#include <BeamAdapter/component/constraint/AdaptiveBeamSlidingConstraint.h>
#include <BeamAdapter/component/WireBeamInterpolation.h>
using sofa::component::constraintset::AdaptiveBeamSlidingConstraint ;
using sofa::component::fem::WireBeamInterpolation ;


namespace sofa
{

template <typename DataTypes>
struct AdaptiveBeamSlidingConstraintTest : public sofa::testing::BaseSimulationTest, AdaptiveBeamSlidingConstraint<DataTypes>
{
    void normalBehavior(){
        Simulation* simu;
        setSimulation(simu = new DAGSimulation());

        typename AdaptiveBeamSlidingConstraint<DataTypes>::SPtr thisObject = New<AdaptiveBeamSlidingConstraint<DataTypes>>();
        thisObject->setName("myname");
        EXPECT_TRUE(thisObject->getName() == "myname");
        EXPECT_TRUE(thisObject->findLink("interpolation") != nullptr );

        Node::SPtr node = simu->createNewGraph("root");
        typename MechanicalObject<DataTypes>::SPtr mecaobject = New<MechanicalObject<DataTypes> >();
        typename WireBeamInterpolation<DataTypes>::SPtr interpolation = New<WireBeamInterpolation<DataTypes> >();

        interpolation->findData("name")->read("wireInterpolation");
        node->addObject(mecaobject);
        node->addObject(interpolation);
        node->addObject(thisObject);
        thisObject->findLink("interpolation")->read("@./wireInterpolation");

        EXPECT_NO_THROW( thisObject->init() );
        EXPECT_NO_THROW( thisObject->reset() );
    }

};

using ::testing::Types;
typedef Types<Rigid3dTypes> DataTypes;

TYPED_TEST_SUITE(AdaptiveBeamSlidingConstraintTest, DataTypes);

TYPED_TEST(AdaptiveBeamSlidingConstraintTest, NormalBehavior) {
    this->normalBehavior() ;
}

}
