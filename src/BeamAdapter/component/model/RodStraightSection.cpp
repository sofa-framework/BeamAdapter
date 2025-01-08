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
#define SOFA_PLUGIN_BEAMADAPTER_RODSTRAIGHTSECTION_CPP

#include <BeamAdapter/component/model/RodStraightSection.inl>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <BeamAdapter/component/model/BaseRodSectionMaterial.inl>

namespace sofa::beamadapter
{

template class SOFA_BEAMADAPTER_API RodStraightSection<sofa::defaulttype::Rigid3Types>;

}// namespace sofa::beamadapter

namespace beamadapter
{

void registerRodStraightSection(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Class defining a rod straight section Material, defining material and geometry parameters.")
                             .add< sofa::beamadapter::RodStraightSection<sofa::defaulttype::Rigid3Types> >());
}

} // namespace beamadapter
