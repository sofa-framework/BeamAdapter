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
#include <BeamAdapter/utils/BeamSection.h>
#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa::beamadapter
{

/**
 * \class WireRestShape
 * \brief Describe the shape functions on multiple segments
 *  
 *  Describe the full shape of a Wire with a given length and radius. The wire is discretized by a set of beams (given by the keyPoints and the relatives Beam density)
 *  This component compute the beam discretization and the shape functions on multiple segments using curvilinear abscissa.
 */
class SOFA_BEAMADAPTER_API WireSectionMaterial : public core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(WireSectionMaterial, core::objectmodel::BaseObject);

    /// Default Constructor
    WireSectionMaterial();


    void init();

    /// This function gives the Young modulus and Poisson's coefficient of the beam depending on the beam position
    void getYoungModulusAtX(SReal& youngModulus, SReal& cPoisson);

    /// This function gives the mass density and the BeamSection data depending on the beam position
    void getInterpolationParam(SReal& _rho, SReal& _A, SReal& _Iy, SReal& _Iz, SReal& _Asy, SReal& _Asz, SReal& _J);

   
    /// User Data about the Young modulus
    Data<SReal> d_poissonRatio;
    Data<SReal> d_youngModulus;

    /// Radius
    Data<SReal> d_radius;
    Data<SReal> d_innerRadius;
    Data<SReal> d_massDensity;
protected:
    BeamSection beamSection;
};

} // namespace sofa::beamadapter
