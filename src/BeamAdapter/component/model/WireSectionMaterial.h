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
    void getYoungModulusAtX(float& youngModulus, float& cPoisson);
    void getYoungModulusAtX(double& youngModulus, double& cPoisson);

    /// This function gives the mass density and the BeamSection data depending on the beam position
    void getInterpolationParam(float& _rho, float& _A, float& _Iy, float& _Iz, float& _Asy, float& _Asz, float& _J);
    void getInterpolationParam(double& _rho, double& _A, double& _Iy, double& _Iz, double& _Asy, double& _Asz, double& _J);

    [[nodiscard]] int getNbVisualEdges() const { return d_nbEdgesVisu.getValue(); }

    [[nodiscard]] int getNbCollisionEdges() const { return d_nbEdgesCollis.getValue(); }
     

   
    /// User Data about the Young modulus
    Data<SReal> d_poissonRatio;
    Data<SReal> d_youngModulus;

    /// Radius
    Data<SReal> d_radius;
    Data<SReal> d_innerRadius;
    Data<SReal> d_massDensity;

    Data< int > d_nbEdgesVisu;
    Data< int > d_nbEdgesCollis;
protected:
    BeamSection beamSection;
};

} // namespace sofa::beamadapter
