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
#pragma once

#include <BeamAdapter/component/model/RodStraightSection.h>
#include <BeamAdapter/component/model/BaseRodSectionMaterial.inl>
#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa::beamadapter
{

template <class DataTypes>
RodStraightSection<DataTypes>::RodStraightSection()
    : BaseRodSectionMaterial<DataTypes>()
{

}


template <class DataTypes>
void RodStraightSection<DataTypes>::initSection()
{

}


template <class DataTypes>
void RodStraightSection<DataTypes>::getRestTransformOnX(Transform& global_H_local, const Real& x_used, const Real& x_start)
{
    global_H_local.set(Vec3(x_start + x_used, 0.0, 0.0), Quat());
}


} // namespace sofa::beamadapter
