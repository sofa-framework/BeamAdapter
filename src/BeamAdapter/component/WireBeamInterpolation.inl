/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
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
#ifndef SOFA_COMPONENT_FEM_WIREBEAMINTERPOLATION_INL
#define SOFA_COMPONENT_FEM_WIREBEAMINTERPOLATION_INL

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/decompose.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/gl/Cylinder.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/gl/Axis.h>
#include <algorithm>

#include <BeamAdapter/component/engine/WireRestShape.h>
#include <BeamAdapter/component/WireBeamInterpolation.h>
#include <BeamAdapter/component/BeamInterpolation.inl>

namespace sofa
{

namespace component
{

namespace fem
{

namespace _wirebeaminterpolation_
{

using sofa::component::engine::WireRestShape ;


template <class DataTypes>
WireBeamInterpolation<DataTypes>::WireBeamInterpolation(sofa::component::engine::WireRestShape<DataTypes> *_restShape)
: Inherited()
, m_restShape(initLink("WireRestShape", "link to the component on the scene"), _restShape)
{


}

template <class DataTypes>
WireBeamInterpolation<DataTypes>::~WireBeamInterpolation()
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
        this->d_componentState = sofa::core::objectmodel::ComponentState::Invalid ;
        return;
    }

    type::vector<Real> xP_noticeable;
    type::vector< int> nbP_density;

    m_restShape.get()->getSamplingParameters(xP_noticeable, nbP_density);
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
    VecElementID &edgeList = *this->d_edgeList.beginEdit();
    vector< double > &lengthList = *this->d_lengthList.beginEdit();
    vector< Transform > &DOF0TransformNode0 = *this->d_DOF0TransformNode0.beginEdit();
    vector< Transform > &DOF1TransformNode1 = *this->d_DOF1TransformNode1.beginEdit();
    vector< Vec2 > &curvAbsList = *this->d_curvAbsList.beginEdit();

    edgeList.push_back(eID);
    lengthList.push_back(length);

    curvAbsList.push_back(Vec2(x0, x1));

    // as an angle is set between DOFs and Beam, they are no more aligned
    this->d_dofsAndBeamsAligned.setValue(false);
    DOF0TransformNode0.push_back(DOF0_H_Node0);
    DOF1TransformNode1.push_back(DOF1_H_Node1);

    this->d_edgeList.endEdit();
    this->d_lengthList.endEdit();
    this->d_DOF0TransformNode0.endEdit();
    this->d_DOF1TransformNode1.endEdit();
    this->d_curvAbsList.endEdit();
}



template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getRestTransform(unsigned int edgeInList, Transform &local0_H_local1_rest)
{
    msg_warning() << "GetRestTransform not implemented for not straightRestShape" ;

    // the beam is straight: the transformation between local0 and local1 is provided by the length of the beam
    local0_H_local1_rest.set(Vec3(this->d_lengthList.getValue()[edgeInList], 0, 0), Quat<Real>());
}


template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest)
{
    if (this->isControlled() && this->m_restShape!=NULL)
    {
        const Vec2 &curvAbs = this->d_curvAbsList.getValue()[edgeInList];

        Real x_middle = (curvAbs.x() + curvAbs.y()) / 2;
        Transform global_H_local_middle, global_H_local_0, global_H_local_1;

        this->m_restShape.get()->getRestTransformOnX(global_H_local_middle, x_middle);
        this->m_restShape.get()->getRestTransformOnX(global_H_local_0, curvAbs.x());
        this->m_restShape.get()->getRestTransformOnX(global_H_local_1, curvAbs.y());

        local_H_local0_rest = global_H_local_middle.inversed() * global_H_local_0;
        local_H_local1_rest = global_H_local_middle.inversed() * global_H_local_1;

        return;
    }

    msg_warning() << "getRestTransform not implemented for not straightRestShape" ;


    /// the beam is straight: local is in the middle of local0 and local1
    /// the transformation between local0 and local1 is provided by the length of the beam
    double edgeMidLength = this->d_lengthList.getValue()[edgeInList] / 2.0;

    local_H_local0_rest.set(-Vec3(edgeMidLength,0,0), Quat<Real>());
    local_H_local1_rest.set(Vec3(edgeMidLength,0,0), Quat<Real>());
}



template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getBeamAtCurvAbs(const Real& x_input, unsigned int &edgeInList_output,
                                                        Real& baryCoord_output, unsigned int start)
{
    if(this->m_brokenInTwo )
    {
        Real x_abs_broken = this->m_restShape.get()->getReleaseCurvAbs();

        ////////// case 1.a : broken part !!
        if (x_input > x_abs_broken)
        {

            /// x_i = curv_abs from the begining of the broken part
            Real x_i = x_input-x_abs_broken;
            Real x=0.0;

            for (unsigned int e=0; e<this->m_numBeamsNotUnderControl; e++)
            {
                x += this->getLength(e);
                if(x > x_i)
                {
                    edgeInList_output = e;
                    Real x0 = x - this->getLength(e);
                    baryCoord_output =(x_i-x0) / this->getLength(e);
                    return;
                }
            }

            edgeInList_output = this->m_numBeamsNotUnderControl-1;
            baryCoord_output = 1.0;
            return;

        }
        ////////// case 1.b : controlled part !!
        else
        {
            start = this->m_numBeamsNotUnderControl;
        }
    }

    Inherited::getBeamAtCurvAbs(x_input, edgeInList_output, baryCoord_output, start);
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
    this->computeTransform2(0, globalHlocal0, globalHlocal1, x);
    Real closestDist = (x_input-globalHlocal0.getOrigin()).norm2();
    Real beamBary = 0.0;
    bool projected = false;
    unsigned int beamIndex = 0;

    // Just look for the closest point on the curve
    // Returns false if this point is not a projection on the curve
    unsigned int nb = this->getNumBeams();
    for(unsigned int i=0; i<nb; i++)	// Check each segment and each vertex
    {
        this->computeTransform2(i, globalHlocal0, globalHlocal1, x);
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
bool WireBeamInterpolation<DataTypes>::getCurvAbsOfProjection(const Vec3& x_input, const VecCoord& vecX, Real& xcurv_output, const Real& tolerance)
{
    // We have put all the code in a new class, because it uses a lot of custom functions and data
    ProjectionSearch<DataTypes> ps(this, x_input, vecX, xcurv_output, tolerance) ;
    return ps.doSearch(xcurv_output) ;
}

template<class DataTypes>
bool WireBeamInterpolation<DataTypes>::breaksInTwo(const Real &x_min_out,  Real &x_break, int &numBeamsNotUnderControlled )
{
    const Real eps = 0.0000000001;

    if (this->m_brokenInTwo)
    {
        msg_error() << "Already broken" ;
        return false;
    }

    if (!this->isControlled() || this->m_restShape == NULL || x_min_out <= eps)
    {
        msg_error() << "Problem with function breaksInTwo ";
        return false;
    }

    // if the release point is not "out" (x_min_out> x_break), then the break is not possible
    x_break = m_restShape.get()->getReleaseCurvAbs();
    if (x_min_out > x_break)
        return false;

    // put the info of the "released" part of the beam in the beginning of the beams;
    this->m_numBeamsNotUnderControl=0;
    unsigned int duplicatePoint=0;

    // browse the curvilinear abscissa to find the point that needs to be duplicate
    // put the info of the second part of the wire at the begining
    unsigned int i=0;

    VecElementID &edgeList = *this->d_edgeList.beginEdit();
    vector< double > &lengthList = *this->d_lengthList.beginEdit();
    vector< Transform > &DOF0TransformNode0 = *this->d_DOF0TransformNode0.beginEdit();
    vector< Transform > &DOF1TransformNode1 = *this->d_DOF1TransformNode1.beginEdit();
    vector< Vec2 > &curvAbsList = *this->d_curvAbsList.beginEdit();

    const unsigned int curvAbsListSize = curvAbsList.size();

    for (unsigned int e = 1; e < curvAbsListSize; e++)
    {
        if (fabs(curvAbsList[e].x() - x_break) < eps)
        {
            duplicatePoint = e;
            this->m_numBeamsNotUnderControl = curvAbsListSize - e;
        }

        if (curvAbsList[e].x() > (x_break - eps))
        {
            edgeList[i] = edgeList[e];
            lengthList[i] = lengthList[e];
            curvAbsList[i] = curvAbsList[e];

            // When the instrument are rotated we apply a transformation between DOF and beam node
            // (should always be the case and dofsAndBeamsAligned should be false)
            if (!this->d_dofsAndBeamsAligned.getValue())
            {
                DOF0TransformNode0[i] = DOF0TransformNode0[e];
                DOF1TransformNode1[i] = DOF1TransformNode1[e];
            }

            i++;
        }
    }

    this->d_edgeList.endEdit();
    this->d_lengthList.endEdit();
    this->d_DOF0TransformNode0.endEdit();
    this->d_DOF1TransformNode1.endEdit();
    this->d_curvAbsList.endEdit();

    if (duplicatePoint == 0)
        msg_error() << " Problem no point were found at the x_break position ! getReleaseCurvAbs() should provide a <<notable>> point" ;

    m_restShape.get()->releaseWirePart();
    numBeamsNotUnderControlled = this->m_numBeamsNotUnderControl;

    this->m_brokenInTwo = true;

    return true;
}

template<class DataTypes>
template<class T>
typename T::SPtr  WireBeamInterpolation<DataTypes>::create(T* tObj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
{
    WireRestShape<DataTypes>* _restShape = NULL;
    std::string _restShapePath;
    bool pathOK = false;

    if(arg)
    {
        if (arg->getAttribute("WireRestShape",NULL) != NULL)
        {
            _restShapePath = arg->getAttribute("WireRestShape");
            context->findLinkDest(_restShape, _restShapePath, NULL);

            if(_restShape == NULL)
              msg_warning(context) << " ("<< tObj->getClassName() <<") : WireRestShape attribute not set correctly, WireBeamInterpolation will be constructed with a default WireRestShape" ;
            else
                pathOK = true;
        }
        else
            msg_error(context) << " (" << tObj->getClassName() <<") : WireRestShape attribute not used, WireBeamInterpolation will be constructed with a default WireRestShape" ;


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

template<class DataTypes>
bool ProjectionSearch<DataTypes>::doSearch(Real& result)
{
    initSearch(m_e);

    // Do one pass of the newton method
    newtonMethod();

    Real dist = computeDistAtCurvAbs(m_e);

    // If the new estimate is good, go with it
    if(m_found)
    {
        result = m_e;
        return testForProjection(result);
    }

    // Look if the new estimate is closer to the target than the previous estimate
    if(dist < m_de)
    {
        // If it is, continue with the newton method, it should converge fast
        while(m_totalIterations < m_interpolation->getNumBeams()+1)
        {
            m_de = dist;
            newtonMethod();
            m_totalIterations++;

            if(m_found)
            {
                // Go back to global curv abs
                result = m_beamStart + (m_beamEnd - m_beamStart) * m_le;
                return testForProjection(result);
            }

            dist = computeDistAtCurvAbs(m_e);
            if(dist > m_de)	// Diverging
                break;

            if(m_le < 0 || m_le > 1)	// We have to change beam, use the dichotomic method instead
                break;
        }
    }

    // If the estimate is outside the beam, or is further than the previous one, change beam
    //  and use a dichotomic search until we find a solution
    if(!m_found)
    {
        m_totalIterations = 0;

        while(m_totalIterations < m_interpolation->getNumBeams() + 10 && m_dichotomicIterations < 10)
        {
            m_totalIterations++;
            m_dichotomicIterations++;

            // We will compute 10 samples
            range = m_segEnd - m_segStart;
            rangeSampling = range / s_sampling;
            for(unsigned int i=0; i<=s_sampling; i++)
            {
                Real curvAbs = m_segStart + rangeSampling * (Real)i;
                m_distTab[i] = computeDistAtCurvAbs(curvAbs);
            }
            unsigned int minIndex = std::min_element(m_distTab, m_distTab + s_sampling + 1) - m_distTab;	// Where is the minimum
            m_e = m_segStart + rangeSampling * minIndex;
            if(testForProjection(m_e))
            {
                result = m_e;
                return true;
            }

            // If the minium is at one extremity, change beam
            if(minIndex == 0)
                changeCurrentBeam(m_beamIndex - 1);
            else if(minIndex == s_sampling)
                changeCurrentBeam(m_beamIndex + 1);
            else
            {	// We continue the search with a smaller interval (keeping only 3 points)
                m_segStart = m_e - rangeSampling;
                m_segEnd = m_e + rangeSampling;
            }
        }
    }

    result = m_e;
    return testForProjection(result);
}

template<class DataTypes>
void ProjectionSearch<DataTypes>::initSearch(Real curvAbs)
{
    m_e = curvAbs;
    if(m_e < 0)
    {
        m_beamIndex = 0;
        m_le = 0;
    }
    else if(m_e > m_interpolation->getRestTotalLength())
    {
        m_beamIndex = m_interpolation->getNumBeams()-1;
        m_le = 1;
    }
    else
        m_interpolation->getBeamAtCurvAbs(m_e, m_beamIndex, m_le);

    m_searchDirection = 0;
    m_interpolation->getAbsCurvXFromBeam(m_beamIndex, m_beamStart, m_beamEnd);
    m_segStart = m_beamStart;
    m_segEnd = m_beamEnd;
    m_interpolation->getSplinePoints(m_beamIndex, m_x, P0, P1, P2, P3);
    m_de = computeDistAtCurvAbs(m_e);
}

template<class DataTypes>
void ProjectionSearch<DataTypes>::newtonMethod()
{
    Real bx = m_le, bx2 = bx*bx, bx3 = bx2*bx;
    Real obx = 1-bx, obx2 = obx*obx, obx3 = obx2*obx;
    Vec3 P = P0*obx3 + P1*3*bx*obx2 + P2*3*bx2*obx + P3*bx3;
    Vec3 dP = P0*(-3*obx2) + P1*(3-12*bx+9*bx2) + P2*(6*bx-9*bx2) + P3*(3*bx2);
    Real f_x = dot( (m_target - P) , dP ) ;

    if (f_x==0.0 || fabs(f_x)/dP.norm() < m_tolerance) // reach convergence
    {
        m_interpolation->getCurvAbsAtBeam(m_beamIndex, m_le, m_e);
        m_found = true;
        return;
    }
    Vec3 d2P = P0*6*(1-bx) + P1*(-12 + 18*bx) + P2*(6-18*bx) + P3*6*bx;

    Real df_x = dot(-dP,dP) + dot((m_target-P), d2P);

    if (fabs(df_x) < 1e-5*m_tolerance)
    {
        return;
    }

    Real d_bx = -f_x/df_x;
    m_le += d_bx;

    m_e = m_beamStart + (m_beamEnd - m_beamStart) * m_le;

    // NOTE : bx+d_bx-1.0 ne donne pas une estimation correcte de la position dans l'autre beam, puisque sa longueur peut �tre diff�rente !
}

template<class DataTypes>
bool ProjectionSearch<DataTypes>::changeCurrentBeam(int index)
{
    // If at the end of the thread
    if(index < 0)
    {
        m_segEnd = m_segStart + rangeSampling;
        return false;
    }
    else if(index > static_cast<int>(m_interpolation->getNumBeams())-1)
    {
        m_segStart = m_segEnd - rangeSampling;
        return false;
    }

    int dir = index - m_beamIndex;
    if(m_searchDirection * dir < 0)
    {	// Changing the direction of search means we are looking for a point near an extremity
        // We know we are looking for a point inside the interval [0.9;1.1]
        // but the ranges for each beam can be different
        Real nStart, nEnd;
        m_interpolation->getAbsCurvXFromBeam(index, nStart, nEnd);
        Real nRangeSampling = (nEnd - nStart) / s_sampling;

        if(dir < 0)
        {
            m_segStart = m_beamStart - nRangeSampling;
            m_segEnd = m_beamStart + rangeSampling;
        }
        else
        {
            m_segStart = m_beamEnd - rangeSampling;
            m_segEnd = m_beamEnd + nRangeSampling;
        }
        return false;
    }


    if(dir < 0)
        m_searchDirection = -1;
    else if(dir > 0)
        m_searchDirection = 1;

    // Really changing beam
    m_beamIndex = index;
    m_interpolation->getAbsCurvXFromBeam(m_beamIndex, m_beamStart, m_beamEnd);
    m_segStart = m_beamStart;
    m_segEnd = m_beamEnd;
    m_interpolation->getSplinePoints(m_beamIndex, m_x, P0, P1, P2, P3);
    m_dichotomicIterations = 0;

    return true;
}

template<class DataTypes>
typename ProjectionSearch<DataTypes>::Real ProjectionSearch<DataTypes>::computeDistAtCurvAbs(Real curvAbs)
{
    if(curvAbs >= m_beamStart && curvAbs <= m_beamEnd)
    {	// We can use the control points we saved
        Real bx = (curvAbs - m_beamStart) / (m_beamEnd - m_beamStart), bx2 = bx*bx, bx3 = bx2*bx;
        Real obx = 1-bx, obx2 = obx*obx, obx3 = obx2*obx;
        Vec3 P = P0*obx3 + P1*3*bx*obx2 + P2*3*bx2*obx + P3*bx3;

        return (m_target-P).norm();
    }
    else
    {
        // TODO(dmarchal 2017-05-17) Please specify who/when this will be done
        // TODO : save all the control points so we don't have to compute them again
        Real bx;
        unsigned int index;
        Vec3 tP0, tP1, tP2, tP3;
        m_interpolation->getBeamAtCurvAbs(curvAbs, index, bx);
        m_interpolation->getSplinePoints(index, m_x, tP0, tP1, tP2, tP3);
        Real bx2 = bx*bx, bx3 = bx2*bx;
        Real obx = 1-bx, obx2 = obx*obx, obx3 = obx2*obx;
        Vec3 P = tP0*obx3 + tP1*3*bx*obx2 + tP2*3*bx2*obx + tP3*bx3;

        return (m_target-P).norm();
    }
}

template<class DataTypes>
bool ProjectionSearch<DataTypes>::testForProjection(Real curvAbs)
{
    Real bx;
    unsigned int index;
    Vec3 tP0, tP1, tP2, tP3;
    m_interpolation->getBeamAtCurvAbs(curvAbs, index, bx);
    m_interpolation->getSplinePoints(index, m_x, tP0, tP1, tP2, tP3);
    Real bx2 = bx*bx, bx3 = bx2*bx;
    Real obx = 1-bx, obx2 = obx*obx, obx3 = obx2*obx;
    Vec3 P = tP0*obx3 + tP1*3*bx*obx2 + tP2*3*bx2*obx + tP3*bx3;
    Vec3 dP = tP0*(-3*obx2) + tP1*(3-12*bx+9*bx2) + tP2*(6*bx-9*bx2) + tP3*(3*bx2);
    Real f_x = dot( (m_target - P) , dP ) ;

    if (fabs(f_x)/dP.norm() < m_tolerance)
        return true;

    return false;
}

} // namespace _wirebeaminterpolation_

} // namespace fem

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_FEM_WIREBEAMINTERPOLATION_INL */
