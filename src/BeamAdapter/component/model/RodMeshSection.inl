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

#define EPSILON 0.0001

namespace sofa::beamadapter
{

template <class DataTypes>
RodMeshSection<DataTypes>::RodMeshSection()
    : BaseRodSectionMaterial<DataTypes>()
    , l_loader(initLink("loader", "link to the MeshLoader"))
{

}


template <class DataTypes>
bool RodMeshSection<DataTypes>::initSection()
{
    // Get meshLoader, check first if loader has been set using link. Otherwise will search in current context.
    p_loader = l_loader.get();

    if (!p_loader)
        this->getContext()->get(p_loader);

    if (!p_loader) {
        msg_error() << "Cannot find a mesh loader. Please insert a MeshObjLoader in the same node or use l_loader to specify the path in the scene graph.";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }
    else
    {
        msg_info() << "Found a mesh with " << p_loader->d_edges.getValue().size() << " edges";
        return initFromLoader();
    }
}


template <class DataTypes>
void RodMeshSection<DataTypes>::getRestTransformOnX(Transform& global_H_local, const Real& x_used, const Real& x_start)
{
    Real abs_curr = x_used - x_start;
    abs_curr = abs_curr /(this->d_length.getValue()) * m_absOfGeometry;

    Coord p;

    /*** find the range which includes the "requested" abs ***/
    Real startingAbs = 0;
    unsigned int index = 0;
    while ((startingAbs < abs_curr) && (index < m_localRestPositions.size()))
    {
        index++;
        startingAbs = m_curvAbs[index];
    }

    /*** OOB ***/
    if (abs_curr > startingAbs)
    {
        msg_error() << "Out of bound position requested= " << abs_curr << " with startingAbs = " << startingAbs;
        return;
    }
    else /*** Expected case ***/
    {
        const Real alpha = (abs_curr - m_curvAbs[index - 1]) / (m_curvAbs[index] - m_curvAbs[index - 1]);
        const Real one_minus_alpha = 1 - alpha;
        const type::Vec3 result = m_localRestTransforms[index - 1].getOrigin() * one_minus_alpha + m_localRestTransforms[index].getOrigin() * alpha;
        Quat slerp;
        slerp.slerp(m_localRestTransforms[index - 1].getOrientation(), m_localRestTransforms[index].getOrientation(), alpha, true);
        slerp.normalize();

        p.getCenter() = result;
        p.getOrientation() = slerp;
    }

    type::Vec3 PosEndCurve = p.getCenter();
    Quat ExtremityQuat = p.getOrientation();
    type::Vec3 ExtremityPos = PosEndCurve + type::Vec3(x_start, 0, 0);

    global_H_local.set(ExtremityPos, ExtremityQuat);
}


template <class DataTypes>
bool RodMeshSection<DataTypes>::initFromLoader()
{
    if (!checkLoaderTopology())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return false;
    }

    type::vector<type::Vec3> vertices;
    sofa::core::topology::BaseMeshTopology::SeqEdges edges;

    //get the topology position
    auto topoVertices = sofa::helper::getReadAccessor(p_loader->d_positions);

    //copy the topology edges in a local vector
    auto topoEdges = sofa::helper::getReadAccessor(p_loader->d_edges);
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
        return false;
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

    initRestConfig();

    return true;
}


template <class DataTypes>
void RodMeshSection<DataTypes>::initRestConfig()
{
    m_curvAbs.clear();
    Real tot = 0;
    m_curvAbs.push_back(0);
    Quat input, output;
    input.identity();
    m_localRestTransforms.resize(m_localRestPositions.size());
    m_localRestTransforms[0].setOrigin(type::Vec3(0, 0, 0));
    m_localRestTransforms[0].setOrientation(input);

    for (unsigned int i = 0; i < m_localRestPositions.size() - 1; i++)
    {
        type::Vec3 vec = m_localRestPositions[i + 1] - m_localRestPositions[i];
        Real norm = vec.norm();
        tot += norm;

        this->rotateFrameForAlignX(input, vec, output);

        input = output;

        m_localRestTransforms[i + 1].setOrientation(output);

        type::Vec3 localPos = m_localRestPositions[i + 1] - m_localRestPositions[0];

        m_localRestTransforms[i + 1].setOrigin(localPos);

        m_curvAbs.push_back(tot);
    }
    m_absOfGeometry = tot;

    this->d_length.setValue(m_absOfGeometry);

    msg_info() << "Length of the loaded shape = " << m_absOfGeometry;
}


template <class DataTypes>
bool RodMeshSection<DataTypes>::checkLoaderTopology()
{
    if (!p_loader->d_edges.getValue().size())
    {
        msg_error() << "There is no edges in the topology loaded by " << p_loader->getName();
        return false;
    }

    if (p_loader->d_triangles.getValue().size())
    {
        msg_error() << "There are triangles in the topology loaded by " << p_loader->getName();
        return false;
    }

    if (p_loader->d_quads.getValue().size())
    {
        msg_error() << "There are quads in the topology loaded by " << p_loader->getName();
        return false;
    }

    if (p_loader->d_polygons.getValue().size())
    {
        msg_error() << "There are polygons in the topology loaded by " << p_loader->getName();
        return false;
    }

    return true;
}


template <class DataTypes>
void RodMeshSection<DataTypes>::rotateFrameForAlignX(const Quat& input, type::Vec3& x, Quat& output)
{
    x.normalize();
    type::Vec3 x0 = input.inverseRotate(x);

    Real cTheta = x0[0];
    Real theta;
    if (cTheta > (1 - EPSILON))
    {
        output = input;
    }
    else
    {
        theta = acos(cTheta);
        // axis of rotation
        type::Vec3 dw(0, -x0[2], x0[1]);
        dw.normalize();

        // computation of the rotation
        Quat inputRoutput;
        inputRoutput.axisToQuat(dw, theta);

        output = input * inputRoutput;
    }
}




} // namespace sofa::beamadapter
