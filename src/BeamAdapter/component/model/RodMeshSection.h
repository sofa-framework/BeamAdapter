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
 * \class WireRestShape
 * \brief Describe the shape functions on multiple segments
 *  
 *  Describe the full shape of a Wire with a given length and radius. The wire is discretized by a set of beams (given by the keyPoints and the relatives Beam density)
 *  This component compute the beam discretization and the shape functions on multiple segments using curvilinear abscissa.
 */
template <class DataTypes>
class RodMeshSection : public sofa::beamadapter::BaseRodSectionMaterial<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(RodMeshSection, DataTypes), SOFA_TEMPLATE(BaseRodSectionMaterial, DataTypes));

    using Coord = typename DataTypes::Coord;
    using Real = typename Coord::value_type;
    using Vec3 = sofa::type::Vec<3, Real>;

    /// Default Constructor
    RodMeshSection();
      
protected:
    void initSection() override;
    void initFromLoader();
    bool checkLoaderTopology();

public:
    /// Link to be set to the topology container in the component graph.
    SingleLink<RodMeshSection<DataTypes>, MeshLoader, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_loader;

private:
    /// Pointer to the MeshLoader, should be set using @sa l_loader, otherwise will search for one in current Node.
    MeshLoader* loader{ nullptr };

    /// Data required for the File loading
    type::vector<Vec3> m_localRestPositions;
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_RODMESHSECTION_CPP)
extern template class SOFA_BEAMADAPTER_API RodMeshSection<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace sofa::beamadapter
