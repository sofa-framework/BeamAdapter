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
#include <gtest/gtest.h>
#include <sofa/simulation/graph/DAGNode.h>
#include <sofa/simpleapi/SimpleApi.h>
#include <sofa/core/ObjectFactory.h>

TEST(AdaptiveBeamController, target)
{
    sofa::simpleapi::importPlugin("BeamAdapter");
    
    const auto node = sofa::simpleapi::createNode("node");
    const auto controller = sofa::simpleapi::createObject(node, "AdaptiveBeamController");

    const auto& creators = sofa::core::ObjectFactory::getInstance()->getEntry("AdaptiveBeamController").creatorMap;

    const auto it = creators.find(sofa::defaulttype::Rigid3Types::Name());
    EXPECT_NE(it, creators.end());

    EXPECT_EQ(std::string(it->second->getTarget()), std::string("BeamAdapter"));
}
