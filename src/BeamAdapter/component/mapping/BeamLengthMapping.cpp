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
// C++ Implementation : BeamLengthMapping
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
#define BEAMADAPTER_BEAMLENGTHMAPPING_CPP
#include <sofa/core/behavior/MechanicalState.h>
#include <BeamAdapter/config.h>
#include <sofa/core/ObjectFactory.h>

#include <BeamAdapter/component/mapping/BeamLengthMapping.inl>


namespace sofa::component::mapping
{

namespace _beamlengthmapping_
{
    template class SOFA_BEAMADAPTER_API BeamLengthMapping<sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec1Types>;
}

} // namespace sofa::component::mapping

namespace beamadapter
{

void registerBeamLengthMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Compute the lengths of the beams.")
                             .add< sofa::component::mapping::BeamLengthMapping<sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec1Types> >());
}

} // namespace beamadapter
