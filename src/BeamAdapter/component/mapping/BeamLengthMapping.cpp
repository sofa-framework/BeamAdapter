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

//using namespace defaulttype;
using namespace core;
using namespace core::behavior;
using namespace sofa::defaulttype;


/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
//SOFA_DECL_CLASS(BeamLengthMapping)

// Register in the Factory
int BeamLengthMappingClass = core::RegisterObject("computes the lengths of the beams")
        .add< BeamLengthMapping<Rigid3Types, Vec1dTypes   > >(true) //default template
        //.add< BeamLengthMapping<Rigid3Types, Rigid3Types > >()
;

namespace _beamlengthmapping_
{
    template class SOFA_BEAMADAPTER_API BeamLengthMapping<Rigid3dTypes, Vec1dTypes   >;
}

} // namespace sofa::component::mapping




