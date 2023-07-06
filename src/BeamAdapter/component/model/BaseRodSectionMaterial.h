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
 * \class BaseRodSectionMaterial
 * \brief Base class describing a Rod section which will define a set of beam elements.
 *  
 * This class provide an api to define a rod/wire section using physical and geometry parameters.
 * The section will then be modellized by a set of beam elements. Inheriting class should provide the geometry structure:
 * @sa RodMeshSection to define a rod using a mesh file, @sa RodSpireSection or @sa RodStraightSection to define procedural shapes.
 * Method @sa initSection and @sa getRestTransformOnX should be overriden to provide the correct creation and interpolation.
 * 
 * The rod section is described by:
 * - Topology parameters: vertices and edges @sa d_nbEdgesVisu and @sa d_nbEdgesCollis
 * - Geometry parameters: radius @sa d_radius, @sa d_innerRadius and length @sa d_length
 * - Mechanical parameters: @sa d_poissonRatio and @sa d_youngModulus
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
    using Size = sofa::Size;

    /////////////////////////// Inherited from BaseObject //////////////////////////////////////////

    /// Default Constructor
    BaseRodSectionMaterial();

    /// init method from BaseObject API. Will call internal @see initSection to be overriden by children
    void init() override;


    /////////////////////////// Geometry and physics Getter //////////////////////////////////////////

    /// Returns the number of visual edges of this section. To be set or computed by child.
    [[nodiscard]] int getNbVisualEdges() const { return d_nbEdgesVisu.getValue(); }

    /// Returns the number of collision edges of this section. To be set or computed by child.
    [[nodiscard]] int getNbCollisionEdges() const { return d_nbEdgesCollis.getValue(); }

    /// Returns the total length of this section. To be set or computed by child.
    [[nodiscard]] Real getLength() const { return d_length.getValue(); }

    /// Returns the Young modulus and Poisson's coefficient of this section
    void getYoungModulusAtX(Real& youngModulus, Real& cPoisson) const;

    /// Returns the mass density and the BeamSection of this section
    void getInterpolationParam(Real& _rho, Real& _A, Real& _Iy, Real& _Iz, Real& _Asy, Real& _Asz, Real& _J) const;

    /// This function is called to get the rest position of the beam depending on the current curved abscisse given in parameter 
    virtual void getRestTransformOnX(Transform& global_H_local, const Real& x_used, const Real& x_start)
    {
        SOFA_UNUSED(global_H_local);
        SOFA_UNUSED(x_used);
        SOFA_UNUSED(x_start);
    }
 
protected:
    /// Internal method to init the section. to be overidden by child.
    virtual bool initSection() { return false; }

public:
    Data<Real> d_poissonRatio; ///< Data defining the mehcanical Poisson ratio of this section
    Data<Real> d_youngModulus; ///< Data defining the mehcanical Young Modulus of this section
    Data<Real> d_massDensity; ///< Data defining the mehcanical mass density of this section
    
    Data<Real> d_radius; ///< Data defining the geometry radius of this section
    Data<Real> d_innerRadius; ///< Data defining the geometry internal radius of this section is hollow 
    Data<Real> d_length; ///< Data defining the geometry length of this section

    Data<Size> d_nbEdgesVisu; ///< Data defining the number of visual edges composing this section
    Data<Size> d_nbEdgesCollis; ///< Data defining the number of collision edges composing this section

private:
    /// Internal structure to store physical parameter of the a beam section
    BeamSection beamSection;
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_BASERODSECTIONMATERIAL_CPP)
extern template class SOFA_BEAMADAPTER_API BaseRodSectionMaterial<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace sofa::beamadapter
