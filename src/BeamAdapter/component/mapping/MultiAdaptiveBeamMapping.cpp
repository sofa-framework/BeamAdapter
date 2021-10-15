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
//
// C++ Implementation : UnifiedMultiMultiAdaptiveBeamMapping
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include <BeamAdapter/component/mapping/MultiAdaptiveBeamMapping.inl>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace mapping
{
using namespace defaulttype;
using namespace core;
using namespace core::behavior;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////

// Register in the Factory
static int MultiAdaptiveBeamMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs")
.add< MultiAdaptiveBeamMapping< Rigid3Types, Vec3Types > >()

;


    template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamMapping< Rigid3Types, Vec3Types >;



} // namespace mapping

} // namespace component

} // namespace sofa

