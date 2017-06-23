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
// C++ Implementation : BeamInterpolation / AdaptiveBeamForceFieldAndMass
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_FEM_BEAMINTERPOLATION_INL
#define SOFA_COMPONENT_FEM_BEAMINTERPOLATION_INL

#include "BeamInterpolation.h"

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/helper/decompose.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <sofa/helper/OptionsGroup.h>

#include <sofa/helper/gl/Cylinder.h>
#include <sofa/simulation/Simulation.h>
#include <sofa/helper/gl/Axis.h>

namespace sofa
{

namespace component
{

namespace fem
{

namespace _beaminterpolation_
{

using core::behavior::MechanicalState;

//////////// useful tool

template <class DataTypes>
void BeamInterpolation<DataTypes>::RotateFrameForAlignX(const Quat &input, Vec3 &x, Quat &output)
{
    x.normalize();
    Vec3 x0=input.inverseRotate(x);

    Real cTheta=x0[0];
     Real theta;
    if (cTheta>0.9999999999)
    {
        output = input;
    }
    else
    {
        theta=acos(cTheta);
        // axis of rotation
        Vec3 dw(0,-x0[2],x0[1]);
        dw.normalize();

        // computation of the rotation
        Quat inputRoutput;
        inputRoutput.axisToQuat(dw, theta);

        output=input*inputRoutput;
    }
}

template <class DataTypes>
BeamInterpolation<DataTypes>::BeamInterpolation() :
   crossSectionShape(initData(&crossSectionShape, sofa::helper::OptionsGroup(3,"circular","elliptic (not available)","rectangular"), "crossSectionShape", "shape of the cross-section. Can be: circular, elliptic, square, rectangular. Default is circular" ))
, radius(initData(&radius, (Real)1.0f, "radius", "radius of the beam (if circular cross-section is considered)"))
, innerRadius(initData(&innerRadius, (Real)0.0f, "innerRadius", "inner radius of the beam if it applies"))
, sideLength(initData(&sideLength, (Real)1.0f, "sideLength", "side length of the beam (if square cross-section is considered)"))
, smallRadius(initData(&smallRadius, (Real)1.0f, "smallRadius", "small radius of the beam (if elliptic cross-section is considered)"))
, largeRadius(initData(&largeRadius, (Real)1.0f, "largeRadius", "large radius of the beam (if elliptic cross-section is considered)"))
, lengthY(initData(&lengthY, (Real)1.0f, "lengthY", "length of the beam section along Y (if rectangular cross-section is considered)"))
, lengthZ(initData(&lengthZ, (Real)1.0f, "lengthZ", "length of the beam section along Z (if rectangular cross-section is considered)"))
, dofsAndBeamsAligned(initData(&dofsAndBeamsAligned, true, "dofsAndBeamsAligned", "if false, a transformation for each beam is computed between the DOF and the beam nodes"))
, defaultYoungModulus(initData(&defaultYoungModulus, (Real) 100000, "defaultYoungModulus", "value of the young modulus if not defined in an other component"))
, straight(initData(&straight,true,"straight","If true, will consider straight beams for the rest position"))
, mStateNodes(sofa::core::objectmodel::New< sofa::component::container::MechanicalObject<sofa::defaulttype::Vec3Types> >())
, m_edgeList(initData(&m_edgeList, "edgeList", "list of the edge in the topology that are concerned by the Interpolation"))
, m_lengthList(initData(&m_lengthList, "lengthList", "list of the length of each beam"))
, m_DOF0TransformNode0(initData(&m_DOF0TransformNode0, "DOF0TransformNode0", "Optional rigid transformation between the degree of Freedom and the first node of the beam"))
, m_DOF1TransformNode1(initData(&m_DOF1TransformNode1, "DOF1TransformNode1", "Optional rigid transformation between the degree of Freedom and the second node of the beam"))
, m_curvAbsList(initData(&m_curvAbsList, "curvAbsList", ""))
, m_beamCollision(initData(&m_beamCollision, "beamCollision", "list of beam (in edgeList) that needs to be considered for collision"))
, d_vecID(initData(&d_vecID, OptionsGroup(3,"current","free","rest" ), "vecID", "input pos and vel (current, free pos/vel, rest pos)" ))
, d_InterpolationInputs(initData(&d_InterpolationInputs, "InterpolationInputs", "vector containing (beamID, baryCoord)"))
, d_InterpolatedPos(initData(&d_InterpolatedPos, "InterpolatedPos", "output Interpolated Position"))
, d_InterpolatedVel(initData(&d_InterpolatedVel, "InterpolatedVel", "output Interpolated Velocity"))
, _topology(NULL)
, _mstate(NULL)
{
    this->brokenInTwo=false;
    _isControlled=false;
    this->_numBeamsNotUnderControl = 0;
    mStateNodes->setName("bezierNodes");
    this->addSlave(mStateNodes);
}


template <class DataTypes>
void BeamInterpolation<DataTypes>::computeCrossSectionInertiaMatrix()
{
    if ( this->crossSectionShape.getValue().getSelectedItem() == "elliptic")
    {
        /* code */
    }
    else if ( this->crossSectionShape.getValue().getSelectedItem() == "rectangular" )
    {
        Real Ly = this->lengthY.getValue();
        Real Lz = this->lengthZ.getValue();

        this->_constantSection._Iy=Ly*Lz*Lz*Lz/12.0;
        this->_constantSection._Iz=Lz*Ly*Ly*Ly/12.0;
        this->_constantSection._J=this->_constantSection._Iy + this->_constantSection._Iz;
        this->_constantSection._A = Ly*Lz;

        this->_constantSection._Asy = 0.0;
        this->_constantSection._Asz = 0.0;
    }
    else
    {
        std::cout << "CROSS SECTION SHAPE !!!! " << this->crossSectionShape.getValue().getSelectedItem() << std::endl;
        this->_constantSection._r = this->radius.getValue();
        this->_constantSection._rInner = this->innerRadius.getValue();
        double r = this->radius.getValue();
        if (r <= 0.0)
        {
            serr << "Radius must be positive" << sendl;
            if (r>0)
            {
                serr << "Radius must be positive. Use absolute value instead. r = " << r << sendl;
                r = -r;
            }
        }
        double rInner = this->innerRadius.getValue();
        this->_constantSection._Iz = M_PI*(r*r*r*r - rInner*rInner*rInner*rInner)/4.0;
        //_Iz = M_PI*(r*r*r*r)/4.0;
        this->_constantSection._Iy = this->_constantSection._Iz ;
        this->_constantSection._J = this->_constantSection._Iz + this->_constantSection._Iy;
        this->_constantSection._A = M_PI*(r*r - rInner*rInner);

        this->_constantSection._Asy = 0.0;
        this->_constantSection._Asz = 0.0;
    }
}


/* ************* ADAPTIVE INTERPOLATION ************** */
template <class DataTypes>
 void BeamInterpolation<DataTypes>::init()
{
    this->f_printLog.setValue(true);

    computeCrossSectionInertiaMatrix();

    sofa::core::objectmodel::BaseContext* context = this->getContext();

    // Init Adaptive Topology:
    this->_topology = context->getMeshTopology();
}


template <class DataTypes>
void BeamInterpolation<DataTypes>::bwdInit()
{
    if(!this->verifyTopology())
        serr<<"no topology found"<<sendl;

    sofa::core::objectmodel::BaseContext* context = this->getContext();

    _mstate = dynamic_cast< sofa::core::behavior::MechanicalState<DataTypes> *> (context->getMechanicalState());

    if (_mstate == NULL)
    {
        serr << "WARNING no MechanicalState found bwdInit failed" << sendl;
        return;
    }

    if (!interpolationIsAlreadyInitialized())
    {
        VecElementID &edgeList = *m_edgeList.beginEdit();

        edgeList.clear();

        for (int i=0; i<this->_topology->getNbEdges(); i++)
        {
            edgeList.push_back(i);
        }

        vector<Transform> &DOF0TransformNode0 = *m_DOF0TransformNode0.beginEdit();
        vector<Transform> &DOF1TransformNode1 = *m_DOF1TransformNode1.beginEdit();

        if (!dofsAndBeamsAligned.getValue())
        {
            DOF0TransformNode0.resize(edgeList.size());
            DOF1TransformNode1.resize(edgeList.size());
        }

        m_edgeList.endEdit();
        m_DOF0TransformNode0.endEdit();
        m_DOF1TransformNode1.endEdit();

        helper::ReadAccessor<Data<VecCoord> > statePos = _mstate->read(sofa::core::ConstVecCoordId::position()) ;

        vector< double > &lengthList = *m_lengthList.beginEdit();
        lengthList.clear();

        const unsigned int edgeListSize = m_edgeList.getValue().size();
        unsigned int nd0Id=0, nd1Id=0;

        for (unsigned int i = 0; i < edgeListSize; i++)
        {
            this->getNodeIndices(i, nd0Id, nd1Id);

            Vec3 beam_segment = statePos[nd1Id].getCenter() - statePos[nd0Id].getCenter();

            lengthList.push_back(beam_segment.norm());
        }

        m_lengthList.endEdit();

        if (!this->verifyTopology())
            serr << "WARNING bwdInit failed" << sendl;
    }
}


template<class DataTypes>
bool BeamInterpolation<DataTypes>::interpolationIsAlreadyInitialized()
{
    if (m_edgeList.getValue().size() == 0)
        return false;

    const unsigned int nbEdges = m_edgeList.getValue().size();

    if (m_DOF0TransformNode0.getValue().size() != nbEdges)
        return false;

    if (m_DOF1TransformNode1.getValue().size() != nbEdges)
        return false;

    if (m_lengthList.getValue().size() != nbEdges)
        return false;

    return true;
}


// verify that we got a non-null pointer on the topology
// verify that the this->m_edgeList do not contain non existing edges

template<class DataTypes>
bool BeamInterpolation<DataTypes>::verifyTopology()
{
    if (this->_topology==NULL)
    {
        serr << "ERROR(BeamInterpolation): object must have a BaseMeshTopology (i.e. EdgeSetTopology or MeshTopology)."<<sendl;
        return false;
    }
    else
    {
        if(this->_topology->getNbEdges()==0)
        {
            serr << "ERROR(BeamInterpolation): topology is empty."<<sendl;
            return false;
        }

        this->_topologyEdges = &this->_topology->getEdges();
        if (this->f_printLog.getValue())
            sout << "the vector _topologyEdges is now set with " << _topologyEdges->size() << " edges" << sendl;


        const VecElementID &edgeList = m_edgeList.getValue();

        for (unsigned int j = 0; j < edgeList.size(); j++)
        {
            if(edgeList[j] > this->_topologyEdges->size())
            {
                serr << "WARNING defined this->m_edgeList is not compatible with topology" << sendl;
                return false;
            }
        }

    }

    return true;
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::clear()
{
    VecElementID &edgeList = *m_edgeList.beginEdit();
    vector< double > &lengthList = *m_lengthList.beginEdit();
    vector< Transform > &DOF0TransformNode0 = *m_DOF0TransformNode0.beginEdit();
    vector< Transform > &DOF1TransformNode1 = *m_DOF1TransformNode1.beginEdit();
    vector< Vec2 > &curvAbsList = *m_curvAbsList.beginEdit();

    if(this->brokenInTwo)
    {
        edgeList.resize(_numBeamsNotUnderControl);
        lengthList.resize(_numBeamsNotUnderControl);
        DOF0TransformNode0.resize(_numBeamsNotUnderControl);
        DOF1TransformNode1.resize(_numBeamsNotUnderControl);
        curvAbsList.resize(_numBeamsNotUnderControl);
    }
    else
    {
        edgeList.clear();
        lengthList.clear();
        DOF0TransformNode0.clear();
        DOF1TransformNode1.clear();
        curvAbsList.clear();
    }

    m_edgeList.endEdit();
    m_lengthList.endEdit();
    m_DOF0TransformNode0.endEdit();
    m_DOF1TransformNode1.endEdit();
    m_curvAbsList.endEdit();
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1, const Real &angle)
{
    VecElementID &edgeList = *m_edgeList.beginEdit();
    vector< double > &lengthList = *m_lengthList.beginEdit();
    vector< Transform > &DOF0TransformNode0 = *m_DOF0TransformNode0.beginEdit();
    vector< Transform > &DOF1TransformNode1 = *m_DOF1TransformNode1.beginEdit();
    vector< Vec2 > &curvAbsList = *m_curvAbsList.beginEdit();

    curvAbsList.push_back(Vec2(x0, x1));

    edgeList.push_back(eID);
    lengthList.push_back(length);

    Quat QuatX ;
    QuatX.axisToQuat(Vec3(1,0,0), angle);
    QuatX.normalize();

    // as an angle is set between DOFs and Beam, they are no more aligned
    this->dofsAndBeamsAligned.setValue(false);
    DOF0TransformNode0.push_back(Transform(Vec3(0, 0, 0), QuatX ));
    DOF1TransformNode1.push_back(Transform(Vec3(0, 0, 0), QuatX ));

    m_edgeList.endEdit();
    m_lengthList.endEdit();
    m_DOF0TransformNode0.endEdit();
    m_DOF1TransformNode1.endEdit();
    m_curvAbsList.endEdit();
}


template<class DataTypes>
 void BeamInterpolation<DataTypes>::getBeamAtCurvAbs(const Real& x_input, unsigned int &edgeInList_output, Real& baryCoord_output, unsigned int start)
 {
     //lTotalRest = total length of the
     Real lTotalRest = this->getRestTotalLength();
     //LTotal =  length sum of the beams that are "out"
     Real LTotal=0.0;

     const unsigned int edgeListSize = m_edgeList.getValue().size();


     // we find the length of the beam that is "out"
     for (unsigned int e = start; e < edgeListSize; e++)
     {
         LTotal += this->getLength(e);
     }


     // x_i = abs_curv from the begining of the instrument
     Real  x_i = x_input + LTotal - lTotalRest;

     if( x_i < 0.0)
     {
         edgeInList_output = start;
         baryCoord_output = 0;
         return;
     }

     // we compute the x value of each node :the topology (stored in Edge_list) is supposed to be a regular seq of segment
     Real x = 0;

     for (unsigned int e = start; e < edgeListSize; e++)
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

     edgeInList_output = edgeListSize - 1;
     baryCoord_output = 1.0;
}


template <class DataTypes>
void BeamInterpolation<DataTypes>::getAbsCurvXFromBeam(int beam, Real& x_curv)
{
    x_curv = m_curvAbsList.getValue()[beam].y();
}

template <class DataTypes>
void BeamInterpolation<DataTypes>::getAbsCurvXFromBeam(int beam, Real& x_curv_start, Real& x_curv_end)
{
    Vec2 ca = m_curvAbsList.getValue()[beam];
    x_curv_start = ca.x();
    x_curv_end = ca.y();
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getDOFtoLocalTransform(unsigned int edgeInList, Transform &DOF0_H_local0, Transform &DOF1_H_local1)
{
    if (dofsAndBeamsAligned.getValue())
    {
        DOF0_H_local0.clear();
        DOF1_H_local1.clear();
        return;
    }

    DOF0_H_local0 = m_DOF0TransformNode0.getValue()[edgeInList];
    DOF1_H_local1 = m_DOF1TransformNode1.getValue()[edgeInList];
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getDOFtoLocalTransformInGlobalFrame(unsigned int edgeInList, Transform &DOF0Global_H_local0, Transform &DOF1Global_H_local1, const VecCoord &x)
{

    Transform DOF0_H_local0, DOF1_H_local1;
    getDOFtoLocalTransform(edgeInList, DOF0_H_local0, DOF1_H_local1);

    unsigned int node0Idx, node1Idx;
    this->getNodeIndices(edgeInList, node0Idx, node1Idx);

    // Computes the Rotation to global from DOF0 and DOF1
    Transform global_R_DOF0(Vec3(0,0,0), x[node0Idx].getOrientation());
    Transform global_R_DOF1(Vec3(0,0,0), x[node1Idx].getOrientation());

    // - rotation due to the optional transformation
    DOF0Global_H_local0 = global_R_DOF0*DOF0_H_local0;
    DOF1Global_H_local1 = global_R_DOF1*DOF1_H_local1;

}


template<class DataTypes>
int BeamInterpolation<DataTypes>::computeTransform(unsigned int edgeInList,
                                                   Transform &global_H0_local,
                                                   Transform &global_H1_local,
                                                   Transform &local0_H_local1,
                                                   Quat &local_R_local0,
                                                   const VecCoord &x)
{

    //<<"computeTransform for edge"<< edgeInList<<std::endl;

    // 1. Get the indices of element and nodes
    unsigned int node0Idx, node1Idx;
    if (getNodeIndices( edgeInList,  node0Idx, node1Idx ) == -1)
    {
        serr << "[computeTransform2] Error in getNodeIndices(). Aborting....." << sendl;
        return -1;
    }

    //2. Computes the optional rigid transformation of DOF0_Transform_node0 and DOF1_Transform_node1
    Transform DOF0_H_local0, DOF1_H_local1;
    getDOFtoLocalTransform(edgeInList, DOF0_H_local0,  DOF1_H_local1);


    // 3. Computes the transformation global To local for both nodes

    Transform global_H_DOF0(x[node0Idx].getCenter(), x[node0Idx].getOrientation());
    Transform global_H_DOF1(x[node1Idx].getCenter(), x[node1Idx].getOrientation());
        // - add a optional transformation
    Transform global_H_local0 = global_H_DOF0*DOF0_H_local0;
    Transform global_H_local1 = global_H_DOF1*DOF1_H_local1;


     // 4. Compute the local frame

    // SIMPLIFICATION: local = local0:
    local_R_local0.clear();

    global_H_DOF0.set(Vec3(0,0,0), x[node0Idx].getOrientation());
    global_H_DOF1.set(Vec3(0,0,0), x[node1Idx].getOrientation());
    // - rotation due to the optional transformation
    global_H_local0 = global_H_DOF0*DOF0_H_local0;
    global_H_local1 = global_H_DOF1*DOF1_H_local1;


    global_H0_local = global_H_local0;
    Quat local0_R_local1 = local0_H_local1.getOrientation();
    Transform local0_HR_local1(Vec3(0,0,0), local0_R_local1);

    global_H1_local = global_H_local1 * local0_HR_local1.inversed();

    return 1;

}


template<class DataTypes>
int BeamInterpolation<DataTypes>::computeTransform2(unsigned int edgeInList,
                                                    Transform &global_H_local0,
                                                    Transform &global_H_local1,
                                                    const VecCoord &x)
{
    // 1. Get the indices of element and nodes
    unsigned int node0Idx, node1Idx;
    if ( getNodeIndices( edgeInList,  node0Idx, node1Idx ) == -1)
    {
        serr << "[computeTransform2] Error in getNodeIndices(). Aborting....." << sendl;
        return -1;
    }

    // 2. Computes the optional rigid transformation of DOF0_Transform_node0 and DOF1_Transform_node1
    Transform DOF0_H_local0, DOF1_H_local1;
    getDOFtoLocalTransform(edgeInList, DOF0_H_local0,  DOF1_H_local1);

    // 3. Computes the transformation global To local for both nodes
    Transform global_H_DOF0(x[node0Idx].getCenter(),x[node0Idx].getOrientation());
    Transform global_H_DOF1(x[node1Idx].getCenter(),x[node1Idx].getOrientation());
    // - add a optional transformation
    global_H_local0 = global_H_DOF0*DOF0_H_local0;
    global_H_local1 = global_H_DOF1*DOF1_H_local1;

    return 1; //without error
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest)
{



    if(straight.getValue())
    {
        // the beam is straight: local is in the middle of local0 and local1
        // the transformation between local0 and local1 is provided by the length of the beam
        local_H_local0_rest.set(-Vec3(this->m_lengthList.getValue()[edgeInList]/2,0,0), Quat());
        local_H_local1_rest.set(Vec3(this->m_lengthList.getValue()[edgeInList]/2,0,0), Quat());
    }
    else
    {
        MechanicalState<DataTypes>* state = dynamic_cast<MechanicalState<DataTypes>*>(getContext()->getMechanicalState());

        if(state)
        {
            unsigned int node0Id, node1Id;
            getNodeIndices(edgeInList,node0Id,node1Id);

            Coord global_0 = state->read(core::VecCoordId::restPosition())->getValue()[node0Id];
            Coord global_1 = state->read(core::VecCoordId::restPosition())->getValue()[node1Id];

            Transform global_H_DOF0 = Transform(global_0.getCenter(),global_0.getOrientation());
            Transform global_H_DOF1 = Transform(global_1.getCenter(),global_1.getOrientation());

            Transform DOF0_H_local0, DOF1_H_local1;
            getDOFtoLocalTransform(edgeInList, DOF0_H_local0,  DOF1_H_local1);

            Transform global_H_local_0 = global_H_DOF0*DOF0_H_local0;
            Transform global_H_local_1 = global_H_DOF1*DOF1_H_local1;


            Transform global_H_local_middle;
            Real baryX = 0.5;
            Real L = this->getLength(edgeInList);

            this->InterpolateTransformUsingSpline(global_H_local_middle, baryX,  global_H_local_0, global_H_local_1, L);



            /*Vec3 result = global_H_local_0.getOrigin() * 0.5 + global_H_local_1.getOrigin() * 0.5;
            Quat slerp;
            slerp.slerp( global_H_local_0.getOrientation(), global_H_local_1.getOrientation(), 0.5, true );
            slerp.normalize();

            global_H_local_middle.setOrigin(result);
            global_H_local_middle.setOrientation(slerp);
                    */

            local_H_local0_rest = global_H_local_middle.inversed() * global_H_local_0;
            local_H_local1_rest = global_H_local_middle.inversed() * global_H_local_1;
        }
        else
            serr<<"This component needs a context mechanical state if straight is set to false."<<sendl;
    }
}


template<class DataTypes>
int BeamInterpolation<DataTypes>::getNodeIndices(unsigned int edgeInList,
                                                 unsigned int &node0Idx,
                                                 unsigned int &node1Idx )
{
    if ( this->_topologyEdges==NULL)
    {
        serr<<"In BeamInterpolation::getNodeIndices() no edges topology defined"<<sendl;
        return -1;
    }


    // 1. Get the indices of element and nodes
    ElementID e = this->m_edgeList.getValue()[edgeInList] ;
    core::topology::BaseMeshTopology::Edge edge=  (*this->_topologyEdges)[e];
    node0Idx = edge[0];
    node1Idx = edge[1];

    return 1;
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getInterpolationParam(unsigned int edgeInList, Real &_L, Real &_A, Real &_Iy ,
                                                                 Real &_Iz, Real &_Asy, Real &_Asz, Real &_J)
{
    // get the length of the beam:
    _L = this->m_lengthList.getValue()[edgeInList];

    BeamSection bS = this->getBeamSection(edgeInList);
    _A=bS._A;
    _Iy=bS._Iy;
    _Iz=bS._Iz;
    _Asy=bS._Asy;
    _Asz=bS._Asz;
    _J=bS._J;
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getSplinePoints(unsigned int edgeInList, const VecCoord &x, Vec3& P0, Vec3& P1, Vec3& P2, Vec3 &P3)
{
    Transform global_H_local0, global_H_local1;
    if (computeTransform2(edgeInList,  global_H_local0,  global_H_local1, x) == -1)
    {
        serr << "[getSplinePoints] error with computeTransform2. Aborting...." << sendl;
        return;
    }

   // << " getSplinePoints  : global_H_local0 ="<<global_H_local0<<"    global_H_local1 ="<<global_H_local1<<std::endl;
    const Real& _L = this->m_lengthList.getValue()[edgeInList];

    P0=global_H_local0.getOrigin();
    P3=global_H_local1.getOrigin();

    P1= P0 + global_H_local0.getOrientation().rotate(Vec3(1.0,0,0))*(_L/3.0);
    P2= P3 + global_H_local1.getOrientation().rotate(Vec3(-1,0,0))*(_L/3.0);
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::computeActualLength(Real &length, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3 &P3)
{
    // the computation of integral Int[0,1] ||dP(x)||  dx = length
    // is done using Gauss Points

    // definition of the Gauss points
    Real x1, x2, x3, x4;
    Real A = 2*sqrt(6.0/5.0);
    x1 = -sqrt((3.0 - A)/7.0 )/2.0+ 0.5;
    x2 = sqrt((3.0 - A) /7.0 )/2.0+ 0.5;
    x3 = -sqrt((3.0 + A)/7.0 )/2.0+ 0.5;
    x4 = sqrt((3.0 + A) /7.0 )/2.0+ 0.5;

    Vec3 dP1, dP2, dP3, dP4;



    dP1 = P0*(-3*(1-x1)*(1-x1)) + P1*(3-12*x1+9*x1*x1) + P2*(6*x1-9*x1*x1) + P3*(3*x1*x1);
    dP2 = P0*(-3*(1-x2)*(1-x2)) + P1*(3-12*x2+9*x2*x2) + P2*(6*x2-9*x2*x2) + P3*(3*x2*x2);
    dP3 = P0*(-3*(1-x3)*(1-x3)) + P1*(3-12*x3+9*x3*x3) + P2*(6*x3-9*x3*x3) + P3*(3*x3*x3);
    dP4 = P0*(-3*(1-x4)*(1-x4)) + P1*(3-12*x4+9*x4*x4) + P2*(6*x4-9*x4*x4) + P3*(3*x4*x4);

    // formula with 4 Gauss Points
    Real B= sqrt(30.0);
    length = ((18.0 + B) /72.0 )*dP1.norm() + ((18.0 + B) /72.0 )*dP2.norm() + ((18.0 - B) /72.0 )*dP3.norm() + ((18.0 - B) /72.0 )*dP4.norm();



    // verification
    /*

    Real length_verif=0.0;
    Vec3 seg, pos;
    pos=P0;
    for (Real bx=0.02; bx<1.00001; bx+=0.02)
    {
        // compute length
        seg  = -pos;
        pos = P0*(1-bx)*(1-bx)*(1-bx) + P1*3*bx*(1-bx)*(1-bx) + P2*3*bx*bx*(1-bx) + P3*bx*bx*bx;
        seg += pos;
        length_verif += seg.norm();
    }


    std::cout<<"computeActualLength verif=  length="<<length<<"  length_verif"<<length_verif<<std::endl;
    */




}


template<class DataTypes>
void BeamInterpolation<DataTypes>::computeStrechAndTwist(unsigned int edgeInList, const VecCoord &x, Vec3 &ResultNodeO, Vec3 &ResultNode1)
{

    // ResultNode = [half length of the beam (m), geometrical Twist angle (rad), additional mechanical Twist angle (rad)]

    //spline:
    Vec3 P0, P1, P2, P3;
    this->getSplinePoints(edgeInList, x, P0, P1, P2, P3);
    ///////// TODO :
    unsigned int node0Idx, node1Idx;
    sout << "in computeStrechAndTwist" << sendl;
    this->getNodeIndices(edgeInList,node0Idx,node1Idx);

    Real length0, length1;
    length0=0.0;
    length1=0.0;

    Vec3 seg;
    Vec3 pos = P0;
    Vec3 n_x, n_y, n_z, x_b, y_b, z_b;
    Quat R0 =  x[node0Idx].getOrientation(); ///// carreful !!! not necessary true !!
    R0.normalize();
    n_x = R0.rotate(Vec3(1.0,0.0,0.0));
    n_y = R0.rotate(Vec3(0.0,1.0,0.0));
    n_z = R0.rotate(Vec3(0.0,0.0,1.0));


    for (Real bx=0.02; bx<1.00001; bx+=0.02)
    {
        // compute length
        seg  = -pos;
        pos = P0*(1-bx)*(1-bx)*(1-bx) + P1*3*bx*(1-bx)*(1-bx) + P2*3*bx*bx*(1-bx) + P3*bx*bx*bx;
        seg += pos;
        if(bx<0.50001)
            length0 += seg.norm();
        else
            length1 += seg.norm();

        // compute frame => Kenneth method
        n_x =  P0*(-3*(1-bx)*(1-bx)) + P1*(3-12*bx+9*bx*bx) + P2*(6*bx-9*bx*bx) + P3*(3*bx*bx);
        n_x.normalize();
        n_z = n_x.cross(n_y);
        n_z.normalize();
        n_y = n_z.cross(n_x);
        n_y.normalize();

        if(bx>0.49999 && bx < 0.50001)
        {
            //bx == 0.5 => frame of the beam (without mechanical twist)
            x_b = n_x;
            y_b = n_y;
            z_b = n_z;
        }
    }

    // computation of twist angle:
    Quat globalRgeom1;
    globalRgeom1 = globalRgeom1.createQuaterFromFrame(n_x,n_y,n_z);
    Vec3 y_geom1 = globalRgeom1.rotate(Vec3(0.0,1.0,0.0));
    Vec3 z_geom1 = globalRgeom1.rotate(Vec3(0.0,0.0,1.0));

    Vec3 x_1, y_1;
    Quat R1 =  x[node1Idx].getOrientation();
    R1.normalize();
    x_1 = R1.rotate(Vec3(1.0,0.0,0.0));
    y_1 = R1.rotate(Vec3(0.0,1.0,0.0));


    //<<" Test : x_1 = "<<x_1<< "  should be equal to ="<< globalRgeom1.rotate(Vec3(1.0,0.0,0.0))<<std::endl;
    Real cosTheta= y_1 * y_geom1;
    Real theta;
    if(cosTheta > 1.0 )
        theta=0.0;
    else
    {
        if (y_1*z_geom1 < 0.0)
            theta = -acos(cosTheta);
        else
            theta= acos(cosTheta);
    }

    ResultNodeO[0]=length0;    ResultNodeO[1]=0.0;  ResultNodeO[2]= theta/2.0;
    ResultNode1[0]=length1;    ResultNode1[1]=0.0;  ResultNode1[2]= theta/2.0;

    //<<" length 0 ="<<length0<<"  length1 ="<<length1 <<"   total length ="<<this->getLength(edgeInList)<<std::endl;

}


template<class DataTypes>
void BeamInterpolation<DataTypes>::interpolatePointUsingSpline(unsigned int edgeInList,
                                                               const Real& baryCoord,
                                                               const Vec3& localPos,
                                                               const VecCoord &x,
                                                               Vec3& posResult,
                                                               bool recompute,
                                                               const sofa::core::ConstVecCoordId &vecXId)
{
    if(recompute)
    {

        // <<" interpolatePointUsingSpline : "<< edgeInList<<"  xbary="<<baryCoord<<"  localPos"<<localPos<<std::endl;
        const Real& _L = this->m_lengthList.getValue()[edgeInList];

        if(localPos.norm() >_L*1e-10)
        {
            Vec3 x_local(0,0,0);
            Transform global_H_localInterpol;

            this->InterpolateTransformUsingSpline(edgeInList, baryCoord, x_local, x, global_H_localInterpol);

            posResult = global_H_localInterpol.getOrigin() + global_H_localInterpol.getOrientation().rotate(localPos);

        }
        else
        {
            /// \todo remove call to computeTransform2 => make something faster !
            Transform global_H_local0, global_H_local1;
            computeTransform2(edgeInList,  global_H_local0,  global_H_local1, x);
            Vec3 DP=global_H_local0.getOrigin() - global_H_local1.getOrigin();

            if( DP.norm()< _L*0.01 )
            {
                posResult = global_H_local0.getOrigin();
                return;
            }


            Vec3 P0,P1,P2,P3;

            P0=global_H_local0.getOrigin();
            P3=global_H_local1.getOrigin();

            P1= P0 + global_H_local0.getOrientation().rotate(Vec3(1.0,0,0))*(_L/3.0);
            P2= P3 + global_H_local1.getOrientation().rotate(Vec3(-1,0,0))*(_L/3.0);

            Real bx=baryCoord;

            posResult = P0*(1-bx)*(1-bx)*(1-bx) + P1*3*bx*(1-bx)*(1-bx) + P2*3*bx*bx*(1-bx) + P3*bx*bx*bx;


        }

    }
    else
    {
        // no need to recompute the positions of the Spline  => we will read their value in the vector which id is vecXId

        const Real& _L = this->m_lengthList.getValue()[edgeInList];

        if(localPos.norm() >_L*1e-10)
        {
            Vec3 x_local(0,0,0);
            Transform global_H_localInterpol;

            this->InterpolateTransformUsingSpline(edgeInList, baryCoord, x_local, x, global_H_localInterpol);

            posResult = global_H_localInterpol.getOrigin() + global_H_localInterpol.getOrientation().rotate(localPos);

        }
        else
        {
            const VecVec3d& v = mStateNodes->read(vecXId)->getValue();
            Vec3 P0,P1,P2,P3;

            P0=v[edgeInList*4];
            P1=v[edgeInList*4+1];
            P2=v[edgeInList*4+2];
            P3=v[edgeInList*4+3];

            Vec3 DP= P0 - P3;

            if( DP.norm()< _L*0.01 )
            {
                posResult = P0;
                return;
            }

            Real bx=baryCoord;

            posResult = P0*(1-bx)*(1-bx)*(1-bx) + P1*3*bx*(1-bx)*(1-bx) + P2*3*bx*bx*(1-bx) + P3*bx*bx*bx;
        }
    }

}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getTangentUsingSplinePoints(unsigned int edgeInList, const Real& baryCoord, const sofa::core::ConstVecCoordId &vecXId, Vec3& t )
{

    const VecVec3d& splinePos = mStateNodes->read(vecXId)->getValue();
    Vec3 P0,P1,P2,P3;

    P0=splinePos[edgeInList*4];
    P1=splinePos[edgeInList*4+1];
    P2=splinePos[edgeInList*4+2];
    P3=splinePos[edgeInList*4+3];

    Real bx=baryCoord;

    t  = P0*-3*(1-bx)*(1-bx) + P1*(3-12*bx+ 9*bx*bx) + P2*(6*bx-9*bx*bx) + P3*3*bx*bx;


}


template<class DataTypes>
void  BeamInterpolation<DataTypes>::updateBezierPoints( const VecCoord &x, sofa::core::VecCoordId &vId_Out){
    //Mechanical Object to stock Bezier points.

//	std::cout<<" in updateBezierPoints vId_Out ="<<vId_Out<<std::endl;

    mStateNodes->resize(this->m_edgeList.getValue().size()*4);
    Data<VecVec3d>* datax = mStateNodes->write(vId_Out);
    VecVec3d& bezierPosVec = *datax->beginEdit();
    bezierPosVec.resize(this->m_edgeList.getValue().size()*4);


    for(unsigned int i=0; i< this->m_edgeList.getValue().size(); i++){
            updateBezierPoints(x,i,bezierPosVec);

    }
    datax->endEdit();
}

template<class DataTypes>
void BeamInterpolation<DataTypes>::updateBezierPoints( const VecCoord &x,unsigned int index, VecVec3d& bezierPosVec){

    // <<" interpolatePointUsingSpline : "<< edgeInList<<"  xbary="<<baryCoord<<"  localPos"<<localPos<<std::endl;
    const Real& _L = this->m_lengthList.getValue()[index];


    /// \todo remove call to computeTransform2 => make something faster !
    Transform global_H_local0, global_H_local1;
    computeTransform2(index,  global_H_local0,  global_H_local1, x);

    //Mechanical Object to stock Bezier points


    bezierPosVec[index*4] =global_H_local0.getOrigin(); //P0
    bezierPosVec[index*4+3]=global_H_local1.getOrigin(); //P3
    bezierPosVec[index*4+1]= bezierPosVec[index*4] + global_H_local0.getOrientation().rotate(Vec3(1.0,0,0))*(_L/3.0); //P1
    bezierPosVec[index*4+2]= bezierPosVec[index*4+3] + global_H_local1.getOrientation().rotate(Vec3(-1,0,0))*(_L/3.0); //P2


}


template<class DataTypes>
void BeamInterpolation<DataTypes>::updateInterpolation(){

    // this function allow to use BeamInterpolation using data (TODO: see if it could fall into engine concept)
    // inputs :
    // -   1.Data<helper::OptionsGroup > d_vecID;
    //         VecID => (1) "current" Pos, Vel    (2) "free" PosFree, VelFree   (3) "rest" PosRest, V=0
    //
    // -   2. Data< vector< Vec2 > > d_InterpolationInputs;
    //         Vector of 2-tuples (indice of the beam   ,   barycentric between 0 and 1)

    // ouptus:
    // -   1. d_InterpolatedPos  => result interpolate the position (6 D0Fs type rigid: position/quaterion)
    // -   2. d_InterpolatedVel  => result interpolate the velocities (linear velocity / angular velocity)


    bool debug=true;

    if (debug)
        std::cout<<"entering updateInterpolation"<<std::endl;

    const vector< Vec2 > &interpolationInputs = d_InterpolationInputs.getValue();
    unsigned int numInterpolatedPositions = interpolationInputs.size();


    VecCoord &interpolatedPos= (*d_InterpolatedPos.beginEdit());
    VecDeriv &interpolatedVel = (*d_InterpolatedVel.beginEdit());

    interpolatedPos.resize(numInterpolatedPositions);
    interpolatedVel.resize(numInterpolatedPositions);



    ///////// select the input position / velocity DATA /////
    if (debug)
        std::cout<<"select the input position  "<<"  selected Item : "<<d_vecID.getValue().getSelectedItem() << std::endl;

    const Data< VecCoord > *x = NULL;
    //const Data< VecDeriv > *v;
    bool computeVel = true;
    if(d_vecID.getValue().getSelectedItem() == "current")
    {
        std::cout<<" position " << std::endl;
        std::cout<< "      ="<< _mstate->read( core::ConstVecCoordId::position() )->getValue( )<<std::endl;
        x=_mstate->read( core::ConstVecCoordId::position() );
        //v=_mstate->read( core::VecDerivId::velocity() ) ;
    }
    else if(d_vecID.getValue().getSelectedItem() == "free")
    {
        x=_mstate->read( core::ConstVecCoordId::freePosition() ) ;
        //v=_mstate->read( core::VecDerivId::freeVelocity() ) ;
    }
    else // rest position
    {
        x=_mstate->read( core::ConstVecCoordId::restPosition() ) ;
        computeVel = false;

    }

    if (debug)
        std::cout<<"selected Item : "<<d_vecID.getValue().getSelectedItem()<<" compute Vel"<< computeVel <<std::endl;


    for (unsigned int i=0; i<numInterpolatedPositions; i++)
    {
        unsigned int numBeam =  (unsigned int) interpolationInputs[i][0];
        Real baryCoord= interpolationInputs[i][1];




        // convert position (RigidTypes) in Transform
        Transform global_H_local0, global_H_local1;
        this->computeTransform2(numBeam, global_H_local0, global_H_local1, x->getValue() );

        // compute the length of the spline
        Vec3 P0,P1,P2,P3;
        this->getSplinePoints(numBeam, x->getValue() , P0,  P1, P2, P3);
        Real length;
        this->computeActualLength(length,P0,P1,P2,P3);

        // get the result of the transform:
        Transform global_H_localResult;


        if(computeVel)
        {
            //this->InterpolateTransformAndVelUsingSpline(global_H_localResult, baryCoord, );
            this->InterpolateTransformUsingSpline(global_H_localResult,baryCoord,global_H_local0,global_H_local1,length);
        }
        else
        {
            this->InterpolateTransformUsingSpline(global_H_localResult,baryCoord,global_H_local0,global_H_local1,length);
        }



        // assign the result in the output data
        interpolatedPos[i].getCenter() = global_H_localResult.getOrigin();
        interpolatedPos[i].getOrientation() = global_H_localResult.getOrientation();



    }
    d_InterpolatedPos.endEdit();
    d_InterpolatedVel.endEdit();

}

/// InterpolateTransformUsingSpline
/// This function provide an interpolation of a frame that is placed between node0 and node1
/// the location of the frame is given by baryCoord
template<class DataTypes>
void BeamInterpolation<DataTypes>::InterpolateTransformUsingSpline(Transform& global_H_localResult,
                                                                           const Real &baryCoord,
                                                                           const Transform &global_H_local0,
                                                                           const Transform &global_H_local1,
                                                                           const Real &L)
{

    Vec3 P0,P1,P2,P3,dP01, dP12, dP03;

    // find the spline points
    P0=global_H_local0.getOrigin();
    P3=global_H_local1.getOrigin();
    P1= P0 + global_H_local0.getOrientation().rotate(Vec3(1.0,0,0))*(L/3.0);
    P2= P3 + global_H_local1.getOrientation().rotate(Vec3(-1,0,0))*(L/3.0);



    Real bx=baryCoord;
    Vec3 posResult;
    Quat quatResult;

    dP01 = P1-P0;
    dP12 = P2-P1;
    dP03 = P3-P0;




    if(dP01*dP12<0.0 && dP03.norm()<0.4*L)
    {

        // The beam is very compressed => it leads to a non correct interpolation using spline
        // (the spline formulation is based on the rest length of the beam)
        // For a correct result, we use linear interpolation instead...
        // (for the quaternion, we use a "simple" slerp


        quatResult = global_H_local0.getOrientation();

        //quatResult.slerp(global_H_local0.getOrientation(),global_H_local1.getOrientation(),bx,true);
        posResult = P0 *(1-bx) + P3*bx;
    }
    else
    {

        // The position of the frame is computed using the interpolation of the spline

         posResult = P0*(1-bx)*(1-bx)*(1-bx) + P1*3*bx*(1-bx)*(1-bx) + P2*3*bx*bx*(1-bx) + P3*bx*bx*bx;

        // the tangent is computed by derivating the spline

        Vec3 n_x =  P0*(-3*(1-bx)*(1-bx)) + P1*(3-12*bx+9*bx*bx) + P2*(6*bx-9*bx*bx) + P3*(3*bx*bx);

        // try to interpolate the "orientation" (especially the torsion) the best possible way... (but other ways should exit...)
        Quat R0, R1, Rslerp;

        //      1. The frame at each node of the beam are rotated to align x along n_x
        this->RotateFrameForAlignX(global_H_local0.getOrientation(), n_x, R0);
        this->RotateFrameForAlignX(global_H_local1.getOrientation(), n_x, R1);

        //     2. We use the "slerp" interpolation to find a solution "in between" these 2 solution R0, R1
        Rslerp.slerp(R0,R1, (float)bx,true);
        Rslerp.normalize();

        //     3. The result is not necessarily alligned with n_x, so we re-aligned Rslerp to obtain a quatResult.
        this->RotateFrameForAlignX(Rslerp, n_x,quatResult);

    }


    global_H_localResult.set(posResult, quatResult);


}

///  getTangent : Computation of a Tangent for the beam (it is given by the derivative of the spline formulation)
template<class DataTypes>
void BeamInterpolation<DataTypes>::getTangent(Vec3& t, const Real& baryCoord, const Transform &global_H_local0, const Transform &global_H_local1,const Real &L)
{

        Vec3 P0,P1,P2,P3, t0, t1;

        P0=global_H_local0.getOrigin();
        P3=global_H_local1.getOrigin();
        P1= P0 + global_H_local0.getOrientation().rotate(Vec3(1.0,0,0))*(L/3.0);
        P2= P3 + global_H_local1.getOrientation().rotate(Vec3(-1,0,0))*(L/3.0);

        t =  P0*(-3*(1-baryCoord)*(1-baryCoord)) + P1*(3-12*baryCoord+9*baryCoord*baryCoord) + P2*(6*baryCoord-9*baryCoord*baryCoord) + P3*(3*baryCoord*baryCoord);

        t.normalize();

}


/// ComputeTotalBendingRotationAngle
/// This function compute a global bending angle value for the current position of the beam
/// Principle: computation of several tangent between node0 to the node1 (given by numComputationPoints)
/// and computation of an angle (acos) between these successive tangents

template<class DataTypes>
typename BeamInterpolation<DataTypes>::Real BeamInterpolation<DataTypes>::ComputeTotalBendingRotationAngle(const Real& dx_computation,
    const Transform &global_H_local0, const Transform &global_H_local1,const Real &L,
    const Real& baryCoordMin, const Real& baryCoordMax)
{

    Vec3 P0,P1,P2,P3, t0, t1;

    P0=global_H_local0.getOrigin();
    P3=global_H_local1.getOrigin();
    P1= P0 + global_H_local0.getOrientation().rotate(Vec3(1.0,0,0))*(L/3.0);
    P2= P3 + global_H_local1.getOrientation().rotate(Vec3(-1,0,0))*(L/3.0);


    t0= P0*(-3*(1-baryCoordMin)*(1-baryCoordMin)) + P1*(3-12*baryCoordMin+9*baryCoordMin*baryCoordMin) + P2*(6*baryCoordMin-9*baryCoordMin*baryCoordMin) + P3*(3*baryCoordMin*baryCoordMin);
    t0.normalize();

    Real BendingAngle=0.0;

    unsigned int numComputationPoints = (unsigned int) ceil((baryCoordMax-baryCoordMin)*(L/dx_computation));
    numComputationPoints = std::max(numComputationPoints, 3u);

    for (unsigned int i=0; i<numComputationPoints; i++)
    {
        Real bx= ((Real)((i+1)/numComputationPoints))*(baryCoordMax-baryCoordMin) + baryCoordMin;
        t1 =  P0*(-3*(1-bx)*(1-bx)) + P1*(3-12*bx+9*bx*bx) + P2*(6*bx-9*bx*bx) + P3*(3*bx*bx);
        t1.normalize();

        if(dot(t0,t1)<1.0)
            BendingAngle += acos(dot(t0,t1));

        t0=t1;
    }

    return BendingAngle;
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::InterpolateTransformUsingSpline(unsigned int edgeInList,
                                                                   const Real& baryCoord,
                                                                   const Vec3& localPos,
                                                                   const VecCoord &x,
                                                                   Transform &global_H_localInterpol)
{
    Transform global_H_local0, global_H_local1;
    if (computeTransform2(edgeInList,  global_H_local0,  global_H_local1, x) == -1)
    {
        serr << "[InterpolateTransformUsingSpline] error with computeTransform2. Aborting...." << sendl;
        return;
    }

    const Real& _L = this->m_lengthList.getValue()[edgeInList];

    if(localPos.norm() >1e-10*_L)
        serr<<"WARNING in InterpolateTransformUsingSpline interpolate frame only on the central curve of the beam"<<sendl;

    InterpolateTransformUsingSpline(global_H_localInterpol, baryCoord, global_H_local0, global_H_local1, _L);
}



template<class DataTypes>
void BeamInterpolation<DataTypes>::InterpolateTransformAndVelUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord &x, const VecDeriv &v,
                                                   Transform &global_H_localInterpol, Deriv &v_interpol)
{

    // 1. Get the indices of element and nodes
    unsigned int node0Idx, node1Idx;
    getNodeIndices( edgeInList,  node0Idx, node1Idx );

    //2. Computes the optional rigid transformation of DOF0_Transform_node0 and DOF1_Transform_node1
    Transform DOF0_H_local0, DOF1_H_local1;
    getDOFtoLocalTransform(edgeInList, DOF0_H_local0,  DOF1_H_local1);


    // 3. Computes the transformation global To local for both nodes
    Transform global_H_DOF0(x[node0Idx].getCenter(),x[node0Idx].getOrientation());
    Transform global_H_DOF1(x[node1Idx].getCenter(),x[node1Idx].getOrientation());

    // - add a optional transformation
    Transform global_H_local0, global_H_local1;
    global_H_local0 = global_H_DOF0*DOF0_H_local0;
    global_H_local1 = global_H_DOF1*DOF1_H_local1;



    // 4. Computes the transformation
    const Real& _L = this->m_lengthList.getValue()[edgeInList];

    if(localPos.norm() >1e-10*_L)
        serr<<"WARNING in InterpolateTransformUsingSpline interpolate frame only on the central curve of the beam"<<sendl;

    InterpolateTransformUsingSpline(global_H_localInterpol, baryCoord, global_H_local0, global_H_local1, _L);


    //--------------- compute velocity interpolation --------------------//


    //5. Computes the velocities of the dof (in their own reference frame...)
    SpatialVector VelDOF0inDOF0, VelDOF1inDOF1;
    VelDOF0inDOF0.setLinearVelocity(   global_H_DOF0.backProjectVector( v[node0Idx].getVCenter() )  );
    VelDOF0inDOF0.setAngularVelocity(  global_H_DOF0.backProjectVector( v[node0Idx].getVOrientation() )  );
    VelDOF1inDOF1.setLinearVelocity(   global_H_DOF1.backProjectVector( v[node1Idx].getVCenter() )  );
    VelDOF1inDOF1.setAngularVelocity(  global_H_DOF1.backProjectVector( v[node1Idx].getVOrientation() )  );


    //6. Computes the velocities of the nodes (in their own reference frame...)
    SpatialVector VelNode0inNode0, VelNode1inNode1;
    VelNode0inNode0= DOF0_H_local0.inversed()*VelDOF0inDOF0;
    VelNode1inNode1= DOF1_H_local1.inversed()*VelDOF1inDOF1;

    //7. Interpolate the result and put in "Deriv" vector...

    v_interpol.getVCenter()= global_H_local0.projectVector(VelNode0inNode0.getLinearVelocity() ) * (1-baryCoord) +
            global_H_local1.projectVector(VelNode1inNode1.getLinearVelocity() ) * baryCoord;

    v_interpol.getVOrientation()= global_H_local0.projectVector(VelNode0inNode0.getAngularVelocity() ) * (1-baryCoord) +
            global_H_local1.projectVector(VelNode1inNode1.getAngularVelocity() ) * baryCoord;


}



////////////////////
// BeamInterpolation<DataTypes>::MapForceOnNodeUsingSpline
// 3DoF (reverse) mapping of a force (finput) that is positionned / to the centerline of a beam
// + point on the centerline given by edgeInList (num of the beam) and baryCoord
// + localPos provides the position shift / to the point on the centerline
// + x provides the positions of the nodes
// Result: Forces (and Torques) on the nodes of this beam
///////////////////

template<class DataTypes>
void BeamInterpolation<DataTypes>::MapForceOnNodeUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord& x, const Vec3& finput,
                                                             SpatialVector& FNode0output, SpatialVector& FNode1output )
{

    //1. get the curvilinear abs and spline parameters
    Real bx = baryCoord;
    Real a0=(1-bx)*(1-bx)*(1-bx);
    Real a1=3*bx*(1-bx)*(1-bx);
    Real a2=3*bx*bx*(1-bx);
    Real a3=bx*bx*bx;

    //2. computes a force on the 4 points of the spline:
    Vec3 F0, F1, F2, F3;
    F0 = finput*a0;
    F1 = finput*a1;
    F2 = finput*a2;
    F3 = finput*a3;

    //3. influence of these forces on the nodes of the beam    => TODO : simplify the computations !!!
    Transform DOF0Global_H_local0, DOF1Global_H_local1;
    this->getDOFtoLocalTransformInGlobalFrame(edgeInList, DOF0Global_H_local0, DOF1Global_H_local1, x);

    //rotate back to local frame
    SpatialVector f0, f1,f2,f3;
    f0.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F0) );
    f1.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F1) );
    f2.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F2) );
    f3.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F3) );

    // computes the torque created on DOF0 and DOF1 by f1 and f2
    Real L = this->getLength(edgeInList);

    if(localPos.norm() > L*1e-10)
    {
            f0.setTorque(localPos.cross(f0.getForce()+f1.getForce()));
            f3.setTorque(localPos.cross(f2.getForce()+f3.getForce()));
    }
    else
    {
            f0.setTorque(Vec3(0,0,0));
            f3.setTorque(Vec3(0,0,0));

    }

    Vec3 lever(L/3,0,0);
    f1.setTorque(lever.cross(f1.getForce()));
    f2.setTorque(-lever.cross(f2.getForce()));


    // back to the DOF0 and DOF1 frame:
    FNode0output = DOF0Global_H_local0 * (f0+f1);
    FNode1output = DOF1Global_H_local1 * (f2+f3);
}


////////////////////
// BeamInterpolation<DataTypes>::MapForceOnNodeUsingSpline
// 6DoF (reverse) mapping of a force (f6DofInput) that is positionned / to the centerline of a beam
// + point on the centerline given by edgeInList (num of the beam) and baryCoord
// + localPos provides the position shift / to the point on the centerline
// + x provides the positions of the nodes
// Result: Forces (and Torques) on the nodes of this beam
///////////////////

template<class DataTypes>
void BeamInterpolation<DataTypes>::MapForceOnNodeUsingSpline(unsigned int edgeInList, const Real& baryCoord, const Vec3& localPos, const VecCoord& x, const SpatialVector& f6DofInput,
                                   SpatialVector& FNode0output, SpatialVector& FNode1output )
{
    //1. get the curvilinear abs and spline parameters
    Real bx = baryCoord;
    Real a0=(1-bx)*(1-bx)*(1-bx);
    Real a1=3*bx*(1-bx)*(1-bx);
    Real a2=3*bx*bx*(1-bx);
    Real a3=bx*bx*bx;

    //2. computes a force on the 4 points of the spline:
    Vec3 F0, F1, F2, F3, C0, C3;
    F0 = f6DofInput.getForce()*a0;
    F1 = f6DofInput.getForce()*a1;
    F2 = f6DofInput.getForce()*a2;
    F3 = f6DofInput.getForce()*a3;
    C0 = f6DofInput.getTorque()*(a0+a1);
    C3 = f6DofInput.getTorque()*(a2+a3);

    //3. influence of these forces on the nodes of the beam    => TODO : simplify the computations !!!
    Transform DOF0Global_H_local0, DOF1Global_H_local1;
    this->getDOFtoLocalTransformInGlobalFrame(edgeInList, DOF0Global_H_local0, DOF1Global_H_local1, x);

    //rotate back to local frame
    SpatialVector f0, f1,f2,f3;
    f0.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F0) );
    f1.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F1) );
    f2.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F2) );
    f3.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F3) );

    // computes the torque created on DOF0 and DOF1 by f1 and f2
    Real L = this->getLength(edgeInList);

    if(localPos.norm() > L*1e-10)
    {
            f0.setTorque(localPos.cross(f0.getForce()+f1.getForce()) + DOF0Global_H_local0.getOrientation().inverseRotate(C0) );
            f3.setTorque(localPos.cross(f2.getForce()+f3.getForce()) + DOF1Global_H_local1.getOrientation().inverseRotate(C3) );
    }
    else
    {
            f0.setTorque( DOF0Global_H_local0.getOrientation().inverseRotate(C0) );
            f3.setTorque( DOF1Global_H_local1.getOrientation().inverseRotate(C3) );

    }

    Vec3 lever(L/3,0,0);
    f1.setTorque(lever.cross(f1.getForce()));
    f2.setTorque(-lever.cross(f2.getForce()));


    // back to the DOF0 and DOF1 frame:
    FNode0output = DOF0Global_H_local0 * (f0+f1);
    FNode1output = DOF1Global_H_local1 * (f2+f3);

}

template<class DataTypes>
bool BeamInterpolation<DataTypes>::breaksInTwo(const Real& /*x_min_out*/,  Real& /*x_break*/, int& /*numBeamsNotUnderControlled*/ )
{
    return false;
}

} /// namespace _beaminterpolation_

} /// namespace fem

} /// namespace component

} /// namespace sofa

#endif  /* SOFA_COMPONENT_FEM_BEAMINTERPOLATION_INL */
