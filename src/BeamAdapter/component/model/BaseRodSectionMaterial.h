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
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/objectmodel/BaseObject.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyModifier.h>
#include <sofa/core/loader/MeshLoader.h>

namespace sofa::beamadapter
{

using sofa::core::loader::MeshLoader;

/**
 * \class WireRestShape
 * \brief Describe the shape functions on multiple segments
 *  
 *  Describe the full shape of a Wire with a given length and radius. The wire is discretized by a set of beams (given by the keyPoints and the relatives Beam density)
 *  This component compute the beam discretization and the shape functions on multiple segments using curvilinear abscissa.
 */
template <class DataTypes>
class BaseRodSectionMaterial : public core::objectmodel::BaseObject
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BaseRodSectionMaterial, DataTypes), core::objectmodel::BaseObject);

    using Coord = typename DataTypes::Coord;
    using Real = typename Coord::value_type;
    using Transform = typename sofa::defaulttype::SolidTypes<Real>::Transform;
    using Vec3 = sofa::type::Vec<3, Real>;
    using Quat = sofa::type::Quat<Real>;

    /// Default Constructor
    BaseRodSectionMaterial();

    void init() override;

    /// This function gives the Young modulus and Poisson's coefficient of the beam depending on the beam position
    void getYoungModulusAtX(Real& youngModulus, Real& cPoisson) const;

    /// This function gives the mass density and the BeamSection data depending on the beam position
    void getInterpolationParam(Real& _rho, Real& _A, Real& _Iy, Real& _Iz, Real& _Asy, Real& _Asz, Real& _J) const;

    /// This function is called by the force field to evaluate the rest position of each beam
    virtual void getRestTransformOnX(Transform& global_H_local, const Real& x_used, const Real& x_start)
    {
        SOFA_UNUSED(global_H_local);
        SOFA_UNUSED(x_used);
        SOFA_UNUSED(x_start);
    }

    [[nodiscard]] int getNbVisualEdges() const { return d_nbEdgesVisu.getValue(); }

    [[nodiscard]] int getNbCollisionEdges() const { return d_nbEdgesCollis.getValue(); }

    [[nodiscard]] Real getLength() const { return d_length.getValue(); }
     
protected:
    virtual void initSection() {}

public:
    /// User Data about the Young modulus
    Data<Real> d_poissonRatio;
    Data<Real> d_youngModulus;

    /// Radius
    Data<Real> d_radius;
    Data<Real> d_innerRadius;
    Data<Real> d_massDensity;

    Data<Real> d_length;

    Data< int > d_nbEdgesVisu;
    Data< int > d_nbEdgesCollis;

private:
    BeamSection beamSection;
};

} // namespace sofa::beamadapter
