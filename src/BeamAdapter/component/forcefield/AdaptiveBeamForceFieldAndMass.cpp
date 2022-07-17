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
// C++ Implementation : WireBeamInterpolation / AdaptiveBeamForceFieldAndMass
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

#define SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEBEAMFORCEFIELD_CPP
#include <BeamAdapter/component/forcefield/AdaptiveBeamForceFieldAndMass.inl>

namespace sofa::component::forcefield
{

namespace _adaptivebeamforcefieldandmass_
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

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
static int AdaptiveBeamForceFieldAndMassClass = RegisterObject("Adaptive Beam finite elements")
.add< AdaptiveBeamForceFieldAndMass<Rigid3Types> >()
;

template class SOFA_BEAMADAPTER_API AdaptiveBeamForceFieldAndMass<Rigid3Types>;

////////////////////////////////////////////////////////////////////////////////////////////////////

} /// _adaptivebeamforcefiedlandmass_

} /// namespace sofa::component::forcefield
