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
//
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#pragma once

#include <BeamAdapter/component/engine/WireRestShape.h>

#include <sofa/core/behavior/MechanicalState.h>

#include <sofa/simulation/TopologyChangeVisitor.h>
#include <sofa/core/visual/VisualParams.h>

#define EPSILON 0.0000000001
#define VERIF 1

namespace beamadapter
{

using sofa::type::vector ;
using sofa::core::objectmodel::TagSet ;
using sofa::core::objectmodel::BaseContext ;

/*!
 * @brief Default Constructor.
 */
template <class DataTypes>
WireRestShape<DataTypes>::WireRestShape()
    : d_density(initData(&d_density, "densityOfBeams", "density of beams between key points"))
    , d_keyPoints(initData(&d_keyPoints,"keyPoints","key points of the shape (curv absc)"))
    , l_sectionMaterials(initLink("wireMaterials", "link to Wire Section Materials (to be ordered according to the instrument, from handle to tip)"))
    , l_topology(initLink("topology", "link to the topology container"))
{
    d_density.setReadOnly(true); // density is supposed to be filled using the section materials
    d_keyPoints.setReadOnly(true); // key points are supposed to be filled using the section materials
}


template<class DataTypes>
void WireRestShape<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Loading);

    //////////////////////////////////////////////
    ////////// get and fill local topology ///////
    //////////////////////////////////////////////

    if (!l_topology)
    {
        l_topology.set(this->getContext()->getMeshTopologyLink());
    }

    if (l_topology)
    {
        msg_info() << "Found topology named "<< l_topology->getName() ;
    }
    else
    {
        msg_error() << "Cannot find topology container. Please specify the link to the topology or insert one in the same node.";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    if (l_sectionMaterials.empty())
    {
        msg_error() << "No BaseRodSectionMaterial set. At least one material should be set and link using wireMaterials.";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }
    
    // Workaround around the fact that d_density and d_keyPoints are supposed to be output only
    if(d_density.isSet())
    {
        msg_warning() << "The density field will be ignored (output only Data).";
    }
    if(d_keyPoints.isSet())
    {
        msg_warning() << "The keyPoints field will be ignored (output only Data).";
    }


    ////////////////////////////////////////////////////////
    ////////// keyPoint list and Density Assignement ///////
    ////////////////////////////////////////////////////////

    initLengths();

  
    initTopology();

    
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    msg_info() << "WireRestShape end init";    
}


template <class DataTypes>
void WireRestShape<DataTypes>::initLengths()
{
    auto keyPointList = sofa::helper::getWriteOnlyAccessor(d_keyPoints);
    auto densityList = sofa::helper::getWriteOnlyAccessor(d_density);
    keyPointList.resize(l_sectionMaterials.size() + 1);
    keyPointList[0] = Real(0.0);
    densityList.resize(l_sectionMaterials.size());
    
    for (unsigned int i = 0; i < l_sectionMaterials.size(); ++i)
    {
        auto rodSection = l_sectionMaterials.get(i);
        keyPointList[i+1] = keyPointList[i] + rodSection->getLength();
        densityList[i] = rodSection->getNbBeams();
    }
}


template <class DataTypes>
bool WireRestShape<DataTypes>::initTopology()
{
    /// fill topology :
    l_topology->clear();
    l_topology->cleanup();

    const type::vector<Real>& keyPts = d_keyPoints.getValue();
    if (l_sectionMaterials.size() != keyPts.size() - 1)
    {
        msg_error() << "Wrong number of inputs. Component can't be init. Number of input materials: " << l_sectionMaterials.size() << ", should be equal to keyPointList.size()-1. keyPointList.size() is equal to: " << keyPts.size();
        return false;
    }
    
    Real prev_length = 0.0;
    int prev_edges = 0;
    int startPtId = 0; 
    for (sofa::Size i = 0; i < l_sectionMaterials.size(); ++i)
    {
        // Add topology of the material 
        int nbrVisuEdges = l_sectionMaterials.get(i)->getNbVisualEdges();
        Real length = fabs(keyPts[i + 1] - keyPts[i]);
        Real dx = length / nbrVisuEdges;

        // add points from the material
        for (int i = startPtId; i < nbrVisuEdges + 1; i++) {
            l_topology->addPoint(prev_length + i * dx, 0, 0);
        }

        // add segments from the material
        for (int i = prev_edges; i < prev_edges + nbrVisuEdges; i++) {
            l_topology->addEdge(i, i + 1);
        }

        prev_length = length;
        prev_edges = nbrVisuEdges;
        startPtId = 1; // Assume the last point of mat[n] == first point of mat[n+1]
    }

    return true;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getSamplingParameters(type::vector<Real>& xP_noticeable,
                                                     type::vector<sofa::Size>& nbP_density) const
{
    xP_noticeable = d_keyPoints.getValue();
    nbP_density = d_density.getValue();
}


template <class DataTypes>
void WireRestShape<DataTypes>::getMechanicalSampling(Real &dx, const Real x_curv)
{
    unsigned int numLines = 0;
    Real x_used = x_curv - EPSILON;

    const Real totalLength = this->getLength();
    x_used = std::clamp(x_used, 0.0, totalLength);

    const type::vector<Real>& keyPts = d_keyPoints.getValue();
    
    // verify that size of number of materials == size of keyPoints-1
    if (l_sectionMaterials.size() != keyPts.size() - 1)
    {
        msg_error() << "Problem size of number of materials: " << l_sectionMaterials.size()
                    << " !=  size of keyPoints-1 " << keyPts.size()-1
                    << ". Returning default values.";
        numLines = 20;
        dx = totalLength / numLines;
        return;
    }
    
    // Check in which section x_used belongs to and get access to this section material
    for (sofa::Size i = 1; i< keyPts.size(); ++i)
    {
        if (x_used <= keyPts[i])
        {
            numLines = l_sectionMaterials.get(i-1)->getNbBeams();

            Real length = fabs(keyPts[i] - keyPts[i-1]);
            dx = length / numLines;
            return;
        }
    }

    // If x_used is out of bounds. Warn user and returns default value.
    numLines = 20;
    dx = totalLength / numLines;
    msg_error() << " problem in getMechanicalSampling : x_curv " << x_used << " is not between keyPoints" << d_keyPoints.getValue();
}

template <class DataTypes>
void WireRestShape<DataTypes>::getCollisionSampling(Real &dx, const Real x_curv)
{
    unsigned int numLines = 0;
    Real x_used = x_curv - EPSILON;

    const Real totalLength = this->getLength();
    if (x_used > totalLength) {
        x_used = totalLength;
    }
    else if (x_used < 0.0) {
        x_used = 0.0;
    }

    const type::vector<Real>& keyPts = d_keyPoints.getValue();
    
    // verify that size of number of materials == size of keyPoints-1
    if (l_sectionMaterials.size() != keyPts.size() - 1)
    {
        msg_error() << "Problem size of number of materials: " << l_sectionMaterials.size()
                    << " !=  size of keyPoints-1 " << keyPts.size()-1 
                    << ". Returning default values.";
        numLines = 20;
        dx = totalLength / numLines;
        return;
    }
    
    // Check in which section x_used belongs to and get access to this section material
    for (sofa::Size i = 1; i< keyPts.size(); ++i)
    {
        if (x_used <= keyPts[i])
        {
            numLines = l_sectionMaterials.get(i-1)->getNbCollisionEdges();

            Real length = fabs(keyPts[i] - keyPts[i-1]);
            dx = length / numLines;
            return;
        }
    }

    // If x_used is out of bounds. Warn user and returns default value.
    numLines = 20;
    dx = totalLength / numLines;
    msg_error() << " problem in getCollisionSampling : x_curv " << x_used << " is not between keyPoints" << d_keyPoints.getValue();
}

template <class DataTypes>
void WireRestShape<DataTypes>::getNumberOfCollisionSegment(Real &dx, sofa::Size& numLines)
{
    numLines = 0;
    for (sofa::Size i = 0; i < l_sectionMaterials.size(); ++i)
    {
        numLines += l_sectionMaterials.get(i)->getNbCollisionEdges();
    }
    dx = getLength() / numLines;
}

template <class DataTypes>
sofa::Size WireRestShape<DataTypes>::getTotalNumberOfBeams()
{
    sofa::Size numBeams = 0;
    for (sofa::Size i = 0; i < l_sectionMaterials.size(); ++i)
    {
        numBeams += l_sectionMaterials.get(i)->getNbBeams();
    }
    
    return numBeams;
}

template <class DataTypes>
void WireRestShape<DataTypes>::getRestTransformOnX(Transform &global_H_local, const Real x)
{
    Real x_used = x - EPSILON;

    const Real totalLength = this->getLength();
    if (x_used > totalLength) {
        x_used = totalLength;
    }
    else if (x_used < 0.0) {
        x_used = 0.0;
    }
   
    const type::vector<Real>& keyPts = d_keyPoints.getValue();
    for (sofa::Size i = 1; i < keyPts.size(); ++i)
    {
        if (x_used <= keyPts[i])
        {
            return l_sectionMaterials.get(i - 1)->getRestTransformOnX(global_H_local, x_used, keyPts[i - 1]);
        }
    }

    msg_error() << " problem in getRestTransformOnX : x_curv " << x_used << " is not between keyPoints" << d_keyPoints.getValue();
}


template <class DataTypes>
const BeamSection& WireRestShape<DataTypes>::getBeamSectionAtX(const Real x_curv) const
{
    const Real x_used = x_curv - Real(EPSILON);
    const type::vector<Real>& keyPts = d_keyPoints.getValue();

    // Check in which section x_used belongs to and get access to this section material
    for (sofa::Size i = 1; i < keyPts.size(); ++i)
    {
        if (x_used <= keyPts[i])
        {
            return l_sectionMaterials.get(i - 1)->getBeamSection();
        }
    }

    msg_error() << " problem in getBeamSection : x_curv " << x_curv << " is not between keyPoints" << keyPts;
    static const BeamSection emptySection{};
    return emptySection;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getInterpolationParametersAtX(const Real x_curv, Real &_A, Real &_Iy , Real &_Iz, Real &_Asy, Real &_Asz, Real &_J) const
{
    const Real x_used = x_curv - Real(EPSILON);
    const type::vector<Real>& keyPts = d_keyPoints.getValue();

    // Check in which section x_used belongs to and get access to this section material
    for (sofa::Size i = 1; i < keyPts.size(); ++i)
    {
        if (x_used <= keyPts[i])
        {
            return l_sectionMaterials.get(i - 1)->getInterpolationParameters(_A, _Iy, _Iz, _Asy, _Asz, _J);
        }
    }

    msg_error() << " problem in getInterpolationParam : x_curv " << x_curv << " is not between keyPoints" << keyPts;
}


template <class DataTypes>
void WireRestShape<DataTypes>::getMechanicalParametersAtX(const Real x_curv, Real& youngModulus, Real& cPoisson, Real& massDensity) const
{
    const Real x_used = x_curv - Real(EPSILON);
    const type::vector<Real>& keyPts = d_keyPoints.getValue();

    // Check in which section x_used belongs to and get access to this section material
    for (std::size_t i = 1; i < keyPts.size(); ++i)
    {
        if (x_used <= keyPts[i])
        {
            return l_sectionMaterials.get(i - 1)->getMechanicalParameters(youngModulus, cPoisson, massDensity);
        }
    }

    msg_error() << " problem in getMechanicalParamAtX : x_curv " << x_curv << " is not between keyPoints" << keyPts;
}


template <class DataTypes>
typename WireRestShape<DataTypes>::Real WireRestShape<DataTypes>::getLength()
{
    return d_keyPoints.getValue().back();
}

template <class DataTypes>
void WireRestShape<DataTypes>::computeOrientation(const Vec3& AB, const Quat& Q, Quat &result)
{
    Vec3 PQ = AB;
    Quat quat = Q;

    Vec3 x = quat.rotate(Vec3(1,0,0));
    PQ.normalize();

    if (dot(x, PQ) > 0.9999999)
        result = Q;

    Vec3 y;
    double alpha;

    if (dot(x, PQ) < -0.9999999)
    {
        y = quat.rotate(Vec3(0,0,1));
        alpha = M_PI;
    }
    else
    {
        y = cross(x, PQ);
        y.normalize();
        alpha = acos(dot(x, PQ));
    }

    Quat qaux = Quat(y, alpha);
    result = qaux * quat;

}


} // namespace beamadapter
