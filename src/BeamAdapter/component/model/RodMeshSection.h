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
#include <sofa/core/loader/MeshLoader.h>

namespace sofa::beamadapter
{

using sofa::core::loader::MeshLoader;

/**
 * \class RodMeshSection
 * \brief Specialization class of @sa BaseRodSectionMaterial describing a rod section created using a Mesh file
 *  
 * This class will describe a rod section defined by a mesh file structure using the link @sa l_loader
 * Method @sa initFromLoader and @sa initRestConfig will define the beam structure using the geometry of the given mesh
 * as well as the Length. Mechanical parameters are set using the @sa BaseRodSectionMaterial Data
 * Method @sa getRestTransformOnX will return the current position of the curviline abscisse along the mesh structure.
 */
template <class DataTypes>
class RodMeshSection : public sofa::beamadapter::BaseRodSectionMaterial<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RodMeshSection, DataTypes), SOFA_TEMPLATE(BaseRodSectionMaterial, DataTypes));

    using Real = typename DataTypes::Real;
    using Transform = typename sofa::defaulttype::SolidTypes<Real>::Transform;
    using Coord = typename DataTypes::Coord;
    using Quat = sofa::type::Quat<Real>;

    /// Default Constructor
    RodMeshSection();

    /// Override method to get the rest position of the beam. In this implementation, it will interpolate along the loaded mesh geometry
    void getRestTransformOnX(Transform& global_H_local, const Real& x_used, const Real& x_start) override;
      
protected:
    /// Internal method to init the section. Called by @sa BaseRodSectionMaterial::init() method
    bool initSection() override;

    /// Internal method called by initSection to init from a linked MeshLoader @sa l_loader
    bool initFromLoader();
    /// Internal method called by initFromLoader to compute @sa m_localRestPositions and @sa m_localRestTransforms given the mesh structure
    void initRestConfig();
    /// Method to check if the given loader has a edge set structure
    bool checkLoaderTopology();

    /// Tool method to rotate the input frame @param input given an axis @param x. Result is set in @param output
    void rotateFrameForAlignX(const type::Quat<Real>& input, type::Vec3& x, type::Quat<Real>& output);

public:
    /// Link to be set to the topology container in the component graph.
    SingleLink<RodMeshSection<DataTypes>, MeshLoader, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_loader;

private:
    /// Pointer to the MeshLoader, should be set using @sa l_loader, otherwise will search for one in current Node.
    MeshLoader* p_loader{ nullptr };

    type::vector<type::Vec3> 		m_localRestPositions; ///< rest position of the key points interpolated on the mesh geometry
    type::vector<Transform> m_localRestTransforms; ///< rest transform of the key points interpolated on the mesh geometry
    type::vector<Real>      m_curvAbs; ///< set of absciss curviline points
    Real 					m_absOfGeometry{ 0 }; ///< max curv absciss of this mesh structure
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_RODMESHSECTION_CPP)
extern template class SOFA_BEAMADAPTER_API RodMeshSection<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace sofa::beamadapter
