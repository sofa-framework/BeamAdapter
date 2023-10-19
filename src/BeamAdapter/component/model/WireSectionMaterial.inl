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

#include <BeamAdapter/component/model/WireSectionMaterial.h>
#include <sofa/core/objectmodel/BaseObject.h>

namespace sofa::beamadapter
{

template <class DataTypes>
WireSectionMaterial<DataTypes>::WireSectionMaterial()
    : d_poissonRatio(initData(&d_poissonRatio, (Real)0.49, "poissonRatio", "Poisson Ratio"))
    , d_youngModulus(initData(&d_youngModulus, (Real)5000, "youngModulus", "Young Modulus"))
    , d_radius(initData(&d_radius, (Real)1.0, "radius", "radius"))
    , d_innerRadius(initData(&d_innerRadius, (Real)0.0, "innerRadius", "inner radius if it applies"))
    , d_massDensity(initData(&d_massDensity, (Real)1.0, "massDensity", "Density of the mass (usually in kg/m^3)"))
    , d_length(initData(&d_length, (Real)1.0, "length", "total length of the wire instrument"))
    , d_density(initData(&d_density, 10, "densityOfBeams", "density of beams between key points"))
    , d_nbEdgesVisu(initData(&d_nbEdgesVisu, 10, "nbEdgesVisu", "number of Edges for the visual model"))
    , d_nbEdgesCollis(initData(&d_nbEdgesCollis, 20, "nbEdgesCollis", "number of Edges for the collision model"))
    , l_loader(initLink("loader", "link to the MeshLoader"))
{

}


template <class DataTypes>
void WireSectionMaterial<DataTypes>::init()
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

    if (!l_loader.empty())
    {
        // Get meshLoader, check first if loader has been set using link. Otherwise will search in current context.
        loader = l_loader.get();
        initFromLoader();
    }
    else
    {
        if (int nbrEdgesVisu = d_nbEdgesVisu.getValue() <= 0)
        {
            msg_warning() << "Number of visual edges has been set to an invalid value: " << nbrEdgesVisu << ". Value should be a positive integer. Setting to default value: 10";
            d_nbEdgesVisu.setValue(10);
        }


        if (int nbEdgesCollis = d_nbEdgesCollis.getValue() <= 0)
        {
            msg_warning() << "Number of collision edges has been set to an invalid value: " << nbEdgesCollis << ". Value should be a positive integer. Setting to default value: 20";
            d_nbEdgesCollis.setValue(10);
        }
    }

    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}


template <class DataTypes>
void WireSectionMaterial<DataTypes>::getYoungModulusAtX(Real& youngModulus, Real& cPoisson) const
{
    youngModulus = this->d_youngModulus.getValue();
    cPoisson = this->d_poissonRatio.getValue();
}


template <class DataTypes>
void WireSectionMaterial<DataTypes>::getInterpolationParam(Real& _rho, Real& _A, Real& _Iy, Real& _Iz, Real& _Asy, Real& _Asz, Real& _J) const
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


template <class DataTypes>
void WireSectionMaterial<DataTypes>::initFromLoader()
{
    if (!checkTopology())
    {
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    // this is... dirty but d_straightLength itself is not relevant for loader-created wire
    // but the code uses it actively and expects it to not be zero.
    if (d_straightLength.getValue() <= 0.0)
    {
        msg_warning() << "straightLength cannot be 0 (or negative...). Setting a minimum value.";
        d_straightLength.setValue(0.0001);
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

    for (unsigned int i = 0; i < m_localRestPositions.size() - 1; i++)
        m_localRestPositions[i] *= d_nonProceduralScale.getValue();

    initRestConfig();
    // TODO epernod 2022-08-05: Init from loader seems quite buggy, need to check if this is still needed and working
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}


} // namespace sofa::beamadapter
