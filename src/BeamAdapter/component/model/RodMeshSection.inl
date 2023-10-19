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

#include <BeamAdapter/component/model/RodMeshSection.h>
#include <BeamAdapter/component/model/BaseRodSectionMaterial.inl>
#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa::beamadapter
{

template <class DataTypes>
RodMeshSection<DataTypes>::RodMeshSection()
    : BaseRodSectionMaterial<DataTypes>()
    , l_loader(initLink("loader", "link to the MeshLoader"))
{

}


template <class DataTypes>
void RodMeshSection<DataTypes>::initSection()
{
    
}


template <class DataTypes>
void RodMeshSection<DataTypes>::initFromLoader()
{
    if (!checkLoaderTopology())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    type::vector<Vec3> vertices;
    sofa::core::topology::BaseMeshTopology::SeqEdges edges;

    //get the topology position
    auto topoVertices = sofa::helper::getReadAccessor(loader->d_positions);

    //copy the topology edges in a local vector
    auto topoEdges = sofa::helper::getReadAccessor(loader->d_edges);
    edges = topoEdges.ref();

    /** renumber the vertices  **/
    type::vector<unsigned int> verticesConnexion; //gives the number of edges connected to a vertex
    for (unsigned int i = 0; i < topoVertices.size(); i++)
        verticesConnexion.push_back(2);

    for (const auto& ed : edges)
    {
        verticesConnexion[ed[0]]--;
        verticesConnexion[ed[1]]--;
    }

    msg_info() << "Successfully compute the vertex connexion";

    // check for the first corner of the edge
    unsigned int firstIndex = 0;
    bool found = false;
    while ((firstIndex < verticesConnexion.size()) && !found)
    {
        if (verticesConnexion[firstIndex] == 1)
            found = true;
        else
            firstIndex++;
    }

    if (firstIndex == verticesConnexion.size())
    {
        msg_error() << "The first vertex of the beam structure is not found, probably because of a closed structure";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    vertices.push_back(topoVertices[firstIndex]);

    while (edges.size() > 0)
    {
        auto it = edges.begin();
        auto end = edges.end();

        bool notFound = true;
        while (notFound && (it != end))
        {
            const auto& ed = (*it);
            auto toDel = it;
            it++;
            if (ed[0] == firstIndex)
            {
                vertices.push_back(topoVertices[ed[1]]);
                firstIndex = ed[1];
                edges.erase(toDel);
                notFound = false;

            }
            else if (ed[1] == firstIndex)
            {
                vertices.push_back(topoVertices[ed[0]]);
                firstIndex = ed[0];
                edges.erase(toDel);
                notFound = false;
            }
        }
    }

    msg_info() << "Successfully computed the topology";

    m_localRestPositions = vertices;

    //for (unsigned int i = 0; i < m_localRestPositions.size() - 1; i++)
    //    m_localRestPositions[i] *= d_nonProceduralScale.getValue();

    //TODO on the WireRestShape
    //initRestConfig();
    // TODO epernod 2022-08-05: Init from loader seems quite buggy, need to check if this is still needed and working
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}


template <class DataTypes>
bool RodMeshSection<DataTypes>::checkLoaderTopology()
{
    if (!loader->d_edges.getValue().size())
    {
        msg_error() << "There is no edges in the topology loaded by " << loader->getName();
        return false;
    }

    if (loader->d_triangles.getValue().size())
    {
        msg_error() << "There are triangles in the topology loaded by " << loader->getName();
        return false;
    }

    if (loader->d_quads.getValue().size())
    {
        msg_error() << "There are quads in the topology loaded by " << loader->getName();
        return false;
    }

    if (loader->d_polygons.getValue().size())
    {
        msg_error() << "There are polygons in the topology loaded by " << loader->getName();
        return false;
    }

    //TODO(dmarchal 2017-05-17) when writing a TODO please specify:
    // who will do that
    // when it will be done
    /// \todo check if the topology is like a wire


    return true;
}



} // namespace sofa::beamadapter
