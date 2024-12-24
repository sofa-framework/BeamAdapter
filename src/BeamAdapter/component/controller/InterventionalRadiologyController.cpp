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
// C++ Implementation : InterventionalRadiologyController
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#define SOFA_PLUGIN_BEAMADAPTER_INTERVENTIONALRADIOCONTROLLER_CPP

#include <BeamAdapter/config.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <BeamAdapter/component/controller/InterventionalRadiologyController.inl>


namespace sofa::component::controller::_interventionalradiologycontroller_
{

template class SOFA_BEAMADAPTER_API InterventionalRadiologyController<sofa::defaulttype::Rigid3Types>;

} // namespace sofa::component::controller::_interventionalradiologycontroller_

namespace beamadapter
{

void registerInterventionalRadiologyController(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Provides a Mouse & Keyboard user control on an EdgeSet Topology.")
                             .add< sofa::component::controller::_interventionalradiologycontroller_::InterventionalRadiologyController<sofa::defaulttype::Rigid3Types> >());
}

} // namespace beamadapter
