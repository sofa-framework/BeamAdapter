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

#include <BeamAdapter/component/model/BaseRodSectionMaterial.h>

namespace beamadapter
{

template <class DataTypes>
BaseRodSectionMaterial<DataTypes>::BaseRodSectionMaterial()
    : d_poissonRatio(initData(&d_poissonRatio, (Real)0.49, "poissonRatio", "Poisson Ratio of this section"))
    , d_youngModulus(initData(&d_youngModulus, (Real)5000, "youngModulus", "Young Modulus of this section"))
    , d_massDensity(initData(&d_massDensity, (Real)1.0, "massDensity", "Density of the mass (usually in kg/m^3)"))
    , d_radius(initData(&d_radius, (Real)1.0, "radius", "Full radius of this section"))
    , d_innerRadius(initData(&d_innerRadius, (Real)0.0, "innerRadius", "Inner radius of this section if hollow"))   
    , d_length(initData(&d_length, (Real)1.0, "length", "Total length of this section"))
    , d_nbBeams(initData(&d_nbBeams, (Size)5, "nbBeams", "Number of Beams for the mechanical model"))
    , d_nbEdgesVisu(initData(&d_nbEdgesVisu, (Size)10, "nbEdgesVisu", "Number of Edges for the visual model"))
    , d_nbEdgesCollis(initData(&d_nbEdgesCollis, (Size)20, "nbEdgesCollis", "Number of Edges for the collision model"))
{

}


template <class DataTypes>
void BaseRodSectionMaterial<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Loading);

    if(!d_nbBeams.isSet())
    {
        msg_deprecated() << "nbBeams is now required but it was not set. Its value will be copied from nbEdgesCollis as a temporary compatibility solution.";
        d_nbBeams.setValue(d_nbEdgesCollis.getValue());
    }
    
    // Prepare beam sections
    double r = this->d_radius.getValue();
    double rInner = this->d_innerRadius.getValue();
    this->m_beamSection._r = r;
    this->m_beamSection._rInner = rInner;
    this->m_beamSection._Iz = M_PI * (r * r * r * r - rInner * rInner * rInner * rInner) / 4.0;
    this->m_beamSection._Iy = this->m_beamSection._Iz;
    this->m_beamSection._J = this->m_beamSection._Iz + this->m_beamSection._Iy;
    this->m_beamSection._A = M_PI * (r * r - rInner * rInner);
    this->m_beamSection._Asy = 0.0;
    this->m_beamSection._Asz = 0.0;

    // call delegate method to init the section
    bool res = initSection();

    if (res)
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    else
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
}


template <class DataTypes>
void BaseRodSectionMaterial<DataTypes>::getInterpolationParameters(Real& _A, Real& _Iy, Real& _Iz, Real& _Asy, Real& _Asz, Real& _J) const
{
    _A = m_beamSection._A; 
    _Iy = m_beamSection._Iy;
    _Iz = m_beamSection._Iz;
    _Asy = m_beamSection._Asy;
    _Asz = m_beamSection._Asz;
    _J = m_beamSection._J;
}


template <class DataTypes>
void BaseRodSectionMaterial<DataTypes>::getMechanicalParameters(Real& youngModulus, Real& cPoisson, Real& massDensity) const
{
    youngModulus = this->d_youngModulus.getValue();
    cPoisson = this->d_poissonRatio.getValue();
    massDensity = this->d_massDensity.getValue();
}

} // namespace beamadapter
