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
// C++ Implementation : AdaptiveBeamController
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

#define SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEBEAMCONTROLLER_CPP

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <BeamAdapter/config.h>
#include <BeamAdapter/component/controller/AdaptiveBeamController.inl>


namespace sofa::component::controller::_adaptivebeamcontroller_
{

using sofa::defaulttype::Rigid3Types;
using sofa::defaulttype::Rigid3Types;
using core::RegisterObject;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////

//TODO(dmarchal 2017-06-01): Il faut remplacer les descriptions dans RegisterObject par un vrai description
static int AdaptiveBeamControllerClass = RegisterObject("")
.add< AdaptiveBeamController<Rigid3Types> >()
;

template class SOFA_BEAMADAPTER_API AdaptiveBeamController<Rigid3Types>;


} // namespace sofa::component::controller::_adaptivebeamcontroller_


