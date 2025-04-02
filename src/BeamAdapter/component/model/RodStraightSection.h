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

namespace beamadapter
{

/**
 * \class RodStraightSection
 * \brief Specialization class of @sa BaseRodSectionMaterial describing a rod straight section.
 *  
 * This class will describe a rod straight section which will parametrized only using the @sa BaseRodSectionMaterial Data
 * Method @sa getRestTransformOnX will return: Vec3(current_x, 0 0)
 */
template <class DataTypes>
class RodStraightSection : public BaseRodSectionMaterial<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RodStraightSection, DataTypes), SOFA_TEMPLATE(BaseRodSectionMaterial, DataTypes));

    using Real = typename DataTypes::Real;
    using Transform = typename sofa::defaulttype::SolidTypes<Real>::Transform;
    using Quat = sofa::type::Quat<Real>;

    /// Default Constructor
    RodStraightSection();

    /// Override method to get the rest position of the beam. In this implementation, it will basically returns Vec3(x_start + x_used, 0 0)
    void getRestTransformOnX(Transform& global_H_local, const Real& x_used, const Real& x_start) override;

protected:
    /// Internal method to init the section. Called by @sa BaseRodSectionMaterial::init() method
    bool initSection() override;
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_RODSTRAIGHTSECTION_CPP)
extern template class SOFA_BEAMADAPTER_API RodStraightSection<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace beamadapter
