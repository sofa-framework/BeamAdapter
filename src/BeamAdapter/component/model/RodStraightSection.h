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
#pragma once

#include <BeamAdapter/config.h>
#include <BeamAdapter/component/model/BaseRodSectionMaterial.h>

namespace sofa::beamadapter
{

/**
 * \class RodStraightSection
 * \brief Describe a rod straight section
 *  
 * 
 */
template <class DataTypes>
class RodStraightSection : public sofa::beamadapter::BaseRodSectionMaterial<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RodStraightSection, DataTypes), SOFA_TEMPLATE(BaseRodSectionMaterial, DataTypes));

    /// Default Constructor
    RodStraightSection();

    void getRestTransformOnX(Transform& global_H_local, const Real& x_used, const Real& x_start) override;

protected:
    void initSection() override;
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_RODSTRAIGHTSECTION_CPP)
extern template class SOFA_BEAMADAPTER_API RodStraightSection<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace sofa::beamadapter
