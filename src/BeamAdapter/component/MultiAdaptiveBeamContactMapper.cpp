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
#define SOFA_COMPONENT_COLLISION_MULTIADAPTIVEBEAMCONTACTMAPPER_CPP


#include "MultiAdaptiveBeamContactMapper.inl"
#include <sofa/helper/Factory.inl>

namespace sofa
{

namespace component
{

namespace collision
{

using namespace defaulttype;



template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamContactMapper<BSplineModel<2>,Vec3Types>;


ContactMapperCreator< ContactMapper<BSplineModel<2> > > BSplineMultiContactMapperClass("MultiAdaptiveBeamContactMapper",true);
template class  ContactMapper<BSplineModel<2> >;



} // namespace collision

} // namespace component

} // namespace sofa

