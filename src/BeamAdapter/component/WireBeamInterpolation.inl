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
// C++ Implementation : WireBeamInterpolation / AdaptiveBeamForceFieldAndMass
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

#include <BeamAdapter/component/WireBeamInterpolation.h>
#include <BeamAdapter/component/BeamInterpolation.inl>


namespace sofa::component::fem::_wirebeaminterpolation_
{

using sofa::component::engine::WireRestShape ;


template <class DataTypes>
WireBeamInterpolation<DataTypes>::WireBeamInterpolation(sofa::component::engine::WireRestShape<DataTypes> *_restShape)
: Inherited()
, m_restShape(initLink("WireRestShape", "link to the component on the scene"), _restShape)
{


}
//////////// useful tool

template <class DataTypes>
void WireBeamInterpolation<DataTypes>::init()
{
    Inherited::init();

    if( m_restShape.get() == nullptr )
    {
        msg_error() << "Missing WireRestShape. The component is thus de-activated" ;
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    type::vector<Real> xP_noticeable;
    type::vector< int> nbP_density;

    m_restShape.get()->getSamplingParameters(xP_noticeable, nbP_density);

    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}


template <class DataTypes>
 void WireBeamInterpolation<DataTypes>::bwdInit()
{
    Inherited::bwdInit();

    if (this->isControlled()){
        msg_info() << "external controller for this ForceField is detected" ;
    }else{
        msg_error() << "not straightRestShape or this->Edge_List is assigned" ;
    }
}


template<class DataTypes>
void WireBeamInterpolation<DataTypes>::addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1,
                                               const Transform &DOF0_H_Node0, const Transform &DOF1_H_Node1)
{
    auto edgeList = sofa::helper::getWriteOnlyAccessor(this->d_edgeList);
    auto lengthList = sofa::helper::getWriteOnlyAccessor(this->d_lengthList);
    auto DOF0TransformNode0 = sofa::helper::getWriteOnlyAccessor(this->d_DOF0TransformNode0);
    auto DOF1TransformNode1 = sofa::helper::getWriteOnlyAccessor(this->d_DOF1TransformNode1);
    auto curvAbsList = sofa::helper::getWriteOnlyAccessor(this->d_curvAbsList);

    edgeList.push_back(eID);
    lengthList.push_back(length);

    curvAbsList.push_back(Vec2(x0, x1));

    // as an angle is set between DOFs and Beam, they are no more aligned
    this->d_dofsAndBeamsAligned.setValue(false);
    DOF0TransformNode0.push_back(DOF0_H_Node0);
    DOF1TransformNode1.push_back(DOF1_H_Node1);
}



template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getRestTransform(unsigned int edgeInList, Transform &local0_H_local1_rest)
{
    msg_warning() << "GetRestTransform not implemented for not straightRestShape" ;

    // the beam is straight: the transformation between local0 and local1 is provided by the length of the beam
    local0_H_local1_rest.set(Vec3(this->d_lengthList.getValue()[edgeInList], 0, 0), Quat());
}


template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest)
{
    if (this->isControlled() && this->m_restShape!=nullptr)
    {
        const Vec2 &curvAbs = this->d_curvAbsList.getValue()[edgeInList];
        auto restShape = this->m_restShape.get();

        Real x_middle = (curvAbs.x() + curvAbs.y()) / 2;
        Transform global_H_local_middle, global_H_local_0, global_H_local_1;

        restShape->getRestTransformOnX(global_H_local_middle, x_middle);
        restShape->getRestTransformOnX(global_H_local_0, curvAbs.x());
        restShape->getRestTransformOnX(global_H_local_1, curvAbs.y());

        local_H_local0_rest = global_H_local_middle.inversed() * global_H_local_0;
        local_H_local1_rest = global_H_local_middle.inversed() * global_H_local_1;

        return;
    }

    msg_warning() << "getRestTransform not implemented for not straightRestShape" ;


    /// the beam is straight: local is in the middle of local0 and local1
    /// the transformation between local0 and local1 is provided by the length of the beam
    double edgeMidLength = this->d_lengthList.getValue()[edgeInList] / 2.0;

    local_H_local0_rest.set(-Vec3(edgeMidLength,0,0), Quat());
    local_H_local1_rest.set(Vec3(edgeMidLength,0,0), Quat());
}


template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getCurvAbsAtBeam(const unsigned int &edgeInList_input, const Real& baryCoord_input, Real& x_output)
{
    ///TODO(dmarchal 2017-05-17): Please tell who and when it will be done.
    // TODO : version plus complete prenant en compte les coupures et autres particularites de ce modele ?
    x_output = 0;
    for(unsigned int i=0; i<edgeInList_input; i++)
        x_output += this->getLength(i);

    x_output += this->getLength(edgeInList_input) * baryCoord_input;
}

template<class DataTypes>
bool WireBeamInterpolation<DataTypes>::getApproximateCurvAbs(const Vec3& x_input, const VecCoord& x, Real& x_output)
{
    if(x.size() <= 1)
    {
        x_output = 0.0;
        return false;
    }

    // Initialize with the first vertex
    Transform globalHlocal0, globalHlocal1;
    this->computeTransform(0, globalHlocal0, globalHlocal1, x);
    Real closestDist = (x_input-globalHlocal0.getOrigin()).norm2();
    Real beamBary = 0.0;
    bool projected = false;
    unsigned int beamIndex = 0;

    // Just look for the closest point on the curve
    // Returns false if this point is not a projection on the curve
    unsigned int nb = this->getNumBeams();
    for(unsigned int i=0; i<nb; i++)	// Check each segment and each vertex
    {
        this->computeTransform(i, globalHlocal0, globalHlocal1, x);
        Vec3 A = globalHlocal0.getOrigin(), B = globalHlocal1.getOrigin();
        Real r = ((x_input-A) * (B-A)) / (B-A).norm2();

        if(r >= 0 && r <= 1)
        {
            Vec3 proj = A + (B-A) * r;
            double dist = (x_input-proj).norm2();
            if(dist < closestDist)
            {
                beamIndex = i;
                beamBary = r;
                projected = true;
                closestDist = dist;
            }
        }
        else if(i != nb-1) // Also check vertices between segments (not the last one)
        {
            double dist = (x_input-B).norm2();
            if(dist < closestDist)
            {
                beamIndex = i;
                beamBary = 1.0;
                projected = true;
                closestDist = dist;
            }
        }
    }

    // Also test the last vertex
    double dist = (x_input-globalHlocal1.getOrigin()).norm2();
    if(dist < closestDist)
    {
        beamIndex = nb - 1;
        beamBary = 1.0;
        projected = false;
    }

    // We know the beam the point can be projected to, translate that to an abscissa
    getCurvAbsAtBeam(beamIndex, beamBary, x_output);
    return projected;
}



template<class DataTypes>
template<class T>
typename T::SPtr  WireBeamInterpolation<DataTypes>::create(T* tObj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
{
    WireRestShape<DataTypes>* _restShape = nullptr;
    std::string _restShapePath;
    bool pathOK = false;

    if(arg)
    {
        if (arg->getAttribute("WireRestShape",nullptr) != nullptr)
        {
            _restShapePath = arg->getAttribute("WireRestShape");
            context->findLinkDest(_restShape, _restShapePath, nullptr);

            if(_restShape == nullptr)
              msg_warning(context) << " ("<< WireBeamInterpolation<DataTypes>::GetClass()->className <<") : WireRestShape attribute not set correctly, WireBeamInterpolation will be constructed with a default WireRestShape" ;
            else
                pathOK = true;
        }
        else
            msg_error(context) << " (" << WireBeamInterpolation<DataTypes>::GetClass()->className <<") : WireRestShape attribute not used, WireBeamInterpolation will be constructed with a default WireRestShape" ;


        if (!pathOK)
        {
            _restShapePath=" ";
            _restShape = new WireRestShape<DataTypes>();
        }
    }

    typename T::SPtr obj = sofa::core::objectmodel::New<T>(_restShape);
    obj->setPathToRestShape(_restShapePath);
    if (context) context->addObject(obj);
    if (arg) obj->parse(arg);
    return obj;
}

} // namespace sofa::component::fem::_wirebeaminterpolation_


