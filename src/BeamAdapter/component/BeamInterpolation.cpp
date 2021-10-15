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
// C++ Implementation : BeamInterpolation / AdaptiveBeamForceFieldAndMass
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

/// This define is here to prevent the declaration of the template instances as "extern".
/// Have a look a the end of BeamInterpolation.h
#define SOFA_BEAMINTERPOLATION_CPP
#include <BeamAdapter/component/BeamInterpolation.inl>

namespace sofa
{

namespace component
{

namespace fem
{

namespace _beaminterpolation_
{

using namespace sofa::defaulttype;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
static int BeamInterpolationClass = core::RegisterObject("Adaptive Beam Interpolation")
.add< BeamInterpolation<Rigid3Types> >()
;



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Explicit template instanciation of extern template.
////////////////////////////////////////////////////////////////////////////////////////////////////
template class SOFA_BEAMADAPTER_API BeamInterpolation<Rigid3Types>;

} /// namespace _beaminterpolation_

} /// namespace fem

} /// namespace component

} /// namespace sofa

