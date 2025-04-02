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
// Author: Eulalie Coevoet
//
// Copyright: See COPYING file that comes with this distribution

#define BEAMADAPTER_MAPPING_BEAMPROJECTIONDIFFERENCEMULTIMAPPING_CPP

#include <BeamAdapter/component/mapping/BeamProjectionDifferenceMultiMapping.inl>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace beamadapter
{

void registerBeamProjectionDifferenceMultiMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Computes the difference between given points and their projection on a beam.")
                             .add< BeamProjectionDifferenceMultiMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types > >());
}

template class SOFA_BEAMADAPTER_API BeamProjectionDifferenceMultiMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;

} // namespace beamadapter
