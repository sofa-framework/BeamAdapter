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

namespace sofa::beamadapter
{

template <class DataTypes>
BaseRodSectionMaterial<DataTypes>::BaseRodSectionMaterial()
    : d_poissonRatio(initData(&d_poissonRatio, (Real)0.49, "poissonRatio", "Poisson Ratio of this section"))
    , d_youngModulus(initData(&d_youngModulus, (Real)5000, "youngModulus", "Young Modulus of this section"))
    , d_massDensity(initData(&d_massDensity, (Real)1.0, "massDensity", "Density of the mass (usually in kg/m^3)"))
    , d_radius(initData(&d_radius, (Real)1.0, "radius", "Full radius of this section"))
    , d_innerRadius(initData(&d_innerRadius, (Real)0.0, "innerRadius", "Inner radius of this section if hollow"))   
    , d_length(initData(&d_length, (Real)1.0, "length", "Total length of this section"))
    , d_nbEdgesVisu(initData(&d_nbEdgesVisu, (Size)10, "nbEdgesVisu", "number of Edges for the visual model"))
    , d_nbEdgesCollis(initData(&d_nbEdgesCollis, (Size)20, "nbEdgesCollis", "number of Edges for the collision model"))
{

}


template <class DataTypes>
void BaseRodSectionMaterial<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Loading);

    // Prepare beam sections
    double r = this->d_radius.getValue();
    double rInner = this->d_innerRadius.getValue();
    this->beamSection._r = r;
    this->beamSection._rInner = rInner;
    this->beamSection._Iz = M_PI * (r * r * r * r - rInner * rInner * rInner * rInner) / 4.0;
    this->beamSection._Iy = this->beamSection._Iz;
    this->beamSection._J = this->beamSection._Iz + this->beamSection._Iy;
    this->beamSection._A = M_PI * (r * r - rInner * rInner);
    this->beamSection._Asy = 0.0;
    this->beamSection._Asz = 0.0;

    // call delegate method to init the section
    bool res = initSection();

    if (res)
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    else
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
}


template <class DataTypes>
void BaseRodSectionMaterial<DataTypes>::getYoungModulusAtX(Real& youngModulus, Real& cPoisson) const
{
    youngModulus = this->d_youngModulus.getValue();
    cPoisson = this->d_poissonRatio.getValue();
}


template <class DataTypes>
void BaseRodSectionMaterial<DataTypes>::getInterpolationParam(Real& _rho, Real& _A, Real& _Iy, Real& _Iz, Real& _Asy, Real& _Asz, Real& _J) const
{
    if (d_massDensity.isSet())
        _rho = d_massDensity.getValue();

    if (d_radius.isSet())
    {
        _A = beamSection._A;
        _Iy = beamSection._Iy;
        _Iz = beamSection._Iz;
        _Asy = beamSection._Asy;
        _Asz = beamSection._Asz;
        _J = beamSection._J;
    }
}

} // namespace sofa::beamadapter
