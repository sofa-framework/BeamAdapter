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
#define SOFA_PLUGIN_BEAMADAPTER_BEAMINTERPOLATION_CPP

#include <sofa/defaulttype/RigidTypes.h>
#include <BeamAdapter/config.h>
#include <sofa/core/ObjectFactory.h>

/// This define is here to prevent the declaration of the template instances as "extern".
/// Have a look a the end of BeamInterpolation.h
#include <BeamAdapter/component/BeamInterpolation.inl>


namespace sofa::component::fem::_beaminterpolation_
{

template class SOFA_BEAMADAPTER_API BeamInterpolation<sofa::defaulttype::Rigid3Types>;

} // namespace sofa::component::fem::_beaminterpolation_

namespace beamadapter
{

void registerBeamInterpolation(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Adaptive Beam Interpolation")
                             .add< sofa::component::fem::_beaminterpolation_::BeamInterpolation<sofa::defaulttype::Rigid3Types> >());
}

} // namespace beamadapter
