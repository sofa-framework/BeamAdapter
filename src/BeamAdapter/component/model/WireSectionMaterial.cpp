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
#define SOFA_PLUGIN_BEAMADAPTER_WIRERESTSHAPE_CPP

#include <BeamAdapter/component/model/WireSectionMaterial.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

namespace sofa::beamadapter
{

using namespace sofa::defaulttype;

int WireSectionMaterialClass = core::RegisterObject("Wire Section Material.")
.add< WireSectionMaterial >();


WireSectionMaterial::WireSectionMaterial()
    : d_poissonRatio(initData(&d_poissonRatio, (SReal)0.49, "poissonRatio", "Poisson Ratio"))
    , d_youngModulus(initData(&d_youngModulus, (SReal)5000, "youngModulus", "Young Modulus"))
    , d_radius(initData(&d_radius, (SReal)1.0f, "radius", "radius"))
    , d_innerRadius(initData(&d_innerRadius, (SReal)0.0f, "innerRadius", "inner radius if it applies"))
    , d_massDensity(initData(&d_massDensity, (SReal)1.0, "massDensity", "Density of the mass (usually in kg/m^3)"))
    , d_nbEdgesVisu(initData(&d_nbEdgesVisu, 10, "nbEdgesVisu", "number of Edges for the visual model"))
    , d_nbEdgesCollis(initData(&d_nbEdgesCollis, 20, "nbEdgesCollis", "number of Edges for the collision model"))
{

}


void WireSectionMaterial::init()
{
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

    if (int nbrEdgesVisu = d_nbEdgesVisu.getValue() <= 0)
    {
        msg_warning() << "Number of visual edges has been set to an invalid value: " << nbrEdgesVisu << ". Value should be a positive integer. Setting to default value: 10";
        d_nbEdgesVisu.setValue(10);
    }


    if (int nbEdgesCollis = d_nbEdgesCollis.getValue() <= 0)
    {
        msg_warning() << "Number of collision edges has been set to an invalid value: " << nbEdgesCollis << ". Value should be a positive integer. Setting to default value: 20";
        d_nbEdgesVisu.setValue(10);
    }
}


/// This function gives the Young modulus and Poisson's coefficient of the beam depending on the beam position
void WireSectionMaterial::getYoungModulusAtX(SReal& youngModulus, SReal& cPoisson)
{
    youngModulus = this->d_youngModulus.getValue();
    cPoisson = this->d_poissonRatio.getValue();
}


/// This function gives the mass density and the BeamSection data depending on the beam position
void WireSectionMaterial::getInterpolationParam(SReal& _rho, SReal& _A, SReal& _Iy, SReal& _Iz, SReal& _Asy, SReal& _Asz, SReal& _J)
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


}// namespace sofa::beamadapter
