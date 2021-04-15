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
//
// C++ Implementation : AdaptivePlasticBeamForceField
//
// Description: Implementation of plasticity over the AdaptiveBeamForceFieldAndMass
// force field interface.
//
//
// Author: Camille Krewcuns, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

#define SOFA_PLUGIN_BEAMADAPTER_ADAPTVEPLASTICBEAMFORCEFIELD_CPP
#include "AdaptivePlasticBeamForceField.inl"

namespace sofa::plugin::beamadapter::component::forcefield
{

namespace _adaptiveplasticbeamforcefield_
{

using sofa::core::RegisterObject ;
using sofa::defaulttype::Rigid3fTypes;
using sofa::defaulttype::Rigid3dTypes;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////

static int AdaptivePlasticBeamForceFieldClass = RegisterObject("Beam element with plastic mechanical behavior and dynamically changing topology")
.add<AdaptivePlasticBeamForceField<Rigid3Types> >();

template class SOFA_BEAMADAPTER_API AdaptivePlasticBeamForceField<Rigid3Types>;

////////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace _adaptiveplasticbeamforcefield_
} // namespace sofa::plugin::beamadapter::component::forcefield