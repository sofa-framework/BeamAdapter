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
#pragma once

#include <BeamAdapter/component/BeamInterpolation.h>
#include <BeamAdapter/component/BaseBeamInterpolation.inl>


namespace beamadapter
{

#define BEAMADAPTER_WITH_VERIFICATION false

using sofa::core::topology::BaseMeshTopology ;
using sofa::core::objectmodel::ComponentState ;
using sofa::core::behavior::MechanicalState;
using sofa::core::objectmodel::BaseContext ;
using sofa::helper::ReadAccessor ;


//////////////////////////////////// BREAMINTERPOLATION ////////////////////////////////////////////
template <class DataTypes>
BeamInterpolation<DataTypes>::BeamInterpolation() :
    crossSectionShape(initData(&crossSectionShape,
                               {"circular","elliptic (not available)","rectangular"},
                               "crossSectionShape",
                               "shape of the cross-section. Can be: circular, elliptic, square, rectangular. Default is circular" ))
  , m_defaultRadius(Real(1.0))
  , d_radius(initData(&d_radius, type::vector<Real>(1, m_defaultRadius), "radius", "radius of the beam (if circular cross-section is considered)"))

  , m_defaultInnerRadius(Real(0.0))
  , d_innerRadius(initData(&d_innerRadius, type::vector<Real>(1, m_defaultInnerRadius), "innerRadius", "inner radius of the beam if it applies"))

  , d_sideLength(initData(&d_sideLength, Real(1.0), "sideLength", "side length of the beam (if square cross-section is considered)"))

  , d_smallRadius(initData(&d_smallRadius, Real(1.0), "smallRadius", "small radius of the beam (if elliptic cross-section is considered)"))
  , d_largeRadius(initData(&d_largeRadius, Real(1.0), "largeRadius", "large radius of the beam (if elliptic cross-section is considered)"))

  , m_defaultLengthY(Real(1.0))
  , d_lengthY(initData(&d_lengthY, type::vector<Real>(1, m_defaultLengthY), "lengthY", "length of the beam section along Y (if rectangular cross-section is considered)"))

  , m_defaultLengthZ(Real(1.0))
  , d_lengthZ(initData(&d_lengthZ, type::vector<Real>(1, m_defaultLengthZ), "lengthZ", "length of the beam section along Z (if rectangular cross-section is considered)"))

  , m_defaultYoungModulus(Real(1e5))
  , d_defaultYoungModulus(initData(&d_defaultYoungModulus, type::vector<Real>(1, m_defaultYoungModulus), "defaultYoungModulus",
                                   "value of the young modulus if not defined in an other component"))

  , m_defaultPoissonRatio(Real(0.4))
  , d_poissonRatio(initData(&d_poissonRatio, type::vector<Real>(1, m_defaultPoissonRatio), "defaultPoissonRatio",
                                   "value of the poisson ratio if not defined in an other component"))
  , d_straight(initData(&d_straight,true,"straight","If true, will consider straight beams for the rest position"))  
  , d_vecID(initData(&d_vecID, {"current","free","rest"}, "vecID",
                     "input pos and vel (current, free pos/vel, rest pos)" ))
  , d_InterpolationInputs(initData(&d_InterpolationInputs, "InterpolationInputs", "vector containing (beamID, baryCoord)"))
  , d_InterpolatedPos(initData(&d_InterpolatedPos, "InterpolatedPos", "output Interpolated Position"))
  , d_InterpolatedVel(initData(&d_InterpolatedVel, "InterpolatedVel", "output Interpolated Velocity")) 
{
    
}

template <class DataTypes>
void BeamInterpolation<DataTypes>::computeRectangularCrossSectionInertiaMatrix(const Real& Ly, const Real& Lz, BeamSection& section)
{
    section._Iy=Ly*Lz*Lz*Lz/12.0;
    section._Iz=Lz*Ly*Ly*Ly/12.0;
    section._J=section._Iy + section._Iz;
    section._A = Ly*Lz;

    section._Asy = 0.0;
    section._Asz = 0.0;
}

template <class DataTypes>
void BeamInterpolation<DataTypes>::computeCircularCrossSectionInertiaMatrix(const Real &r, const Real &rInner, BeamSection &section)
{
    section._r = r;
    section._rInner = rInner;

    section._Iz = M_PI*(r*r*r*r - rInner*rInner*rInner*rInner)/4.0;

    ///_Iz = M_PI*(r*r*r*r)/4.0;
    section._Iy = section._Iz ;
    section._J = section._Iz + section._Iy;
    section._A = M_PI*(r*r - rInner*rInner);

    section._Asy = 0.0;
    section._Asz = 0.0;
}

template <class DataTypes>
void BeamInterpolation<DataTypes>::computeCrossSectionInertiaMatrix()
{
    msg_info() << "Cross section shape:" << crossSectionShape.getValue().getSelectedItem() ;
    if (m_constantBeam)
    {
        if ( crossSectionShape.getValue().getSelectedItem() == "elliptic")
        {
            /* code */
        }
        else if ( crossSectionShape.getValue().getSelectedItem() == "rectangular" )
        {
            computeRectangularCrossSectionInertiaMatrix(m_defaultLengthY, m_defaultLengthZ, m_constantSection);
        }
        else
        {
            computeCircularCrossSectionInertiaMatrix(m_defaultRadius, m_defaultInnerRadius, m_constantSection);
        }
    }
    else
    {
        Size nbEdges = this->m_topology->getNbEdges();
        m_section.resize(nbEdges);
        if ( crossSectionShape.getValue().getSelectedItem() == "elliptic")
        {
            /* code */
        }
        else if ( crossSectionShape.getValue().getSelectedItem() == "rectangular" )
        {
            const auto& lengthY = helper::getReadAccessor(d_lengthY);
            const auto& lengthZ = helper::getReadAccessor(d_lengthZ);
            for (sofa::Size beamId=0; beamId<nbEdges; beamId++)
                computeRectangularCrossSectionInertiaMatrix(lengthY[beamId], lengthZ[beamId], m_section[beamId]);
        }
        else
        {   const auto& radius = helper::getReadAccessor(d_radius);
            const auto& innerRadius = helper::getReadAccessor(d_innerRadius);
            for (sofa::Size beamId=0; beamId<nbEdges; beamId++)
                computeCircularCrossSectionInertiaMatrix(radius[beamId], innerRadius[beamId], m_section[beamId]);
        }
    }
}


////////////////////////////////// ADAPTIVE INTERPOLATION //////////////////////////////////////////
template <class DataTypes>
void BeamInterpolation<DataTypes>::init()
{
}

template <class DataTypes>
void BeamInterpolation<DataTypes>::checkDataSize(Real& defaultValue, Data<type::vector<Real>>& dataList, const size_t& nbEdges)
{
    auto values = sofa::helper::getWriteOnlyAccessor(dataList);
    if (values.size() != nbEdges)
    {
        Real value = defaultValue;
        if (values.size() == 0)
        {
            msg_warning() << "Empty data field for " << dataList.getName() <<". Set default " << value;
        } else
        {
            value = values[0];
        }
        values.resize(nbEdges);
        for (auto& _value: values)
            _value = value;
    }
    else
    {
        m_constantBeam = false;
    }
    defaultValue = values[0]; // if the sizes mismatch again at runtime, will use this default value
}

template <class DataTypes>
void BeamInterpolation<DataTypes>::bwdInit()
{
    this->d_componentState.setValue(ComponentState::Loading);
    BaseBeamInterpolation<DataTypes>::bwdInit();

    if (this->d_componentState.getValue() == ComponentState::Invalid)
        return;

    Size nbEdges = this->m_topology->getNbEdges();
    checkDataSize(m_defaultRadius, d_radius, nbEdges);
    checkDataSize(m_defaultInnerRadius, d_innerRadius, nbEdges);
    checkDataSize(m_defaultLengthY, d_lengthY, nbEdges);
    checkDataSize(m_defaultLengthZ, d_lengthZ, nbEdges);
    checkDataSize(m_defaultYoungModulus, d_defaultYoungModulus, nbEdges);
    checkDataSize(m_defaultPoissonRatio, d_poissonRatio, nbEdges);
    computeCrossSectionInertiaMatrix();

    if (!interpolationIsAlreadyInitialized())
    {
        auto edgeList = sofa::helper::getWriteOnlyAccessor(this->d_edgeList);
        edgeList.clear();

        for (unsigned int i=0; i<this->m_topology->getNbEdges(); i++)
        {
            edgeList.push_back(i);
        }

        auto DOF0TransformNode0 = sofa::helper::getWriteOnlyAccessor(this->d_DOF0TransformNode0);
        auto DOF1TransformNode1 = sofa::helper::getWriteOnlyAccessor(this->d_DOF1TransformNode1);

        if (!d_dofsAndBeamsAligned.getValue())
        {
            DOF0TransformNode0.resize(edgeList.size());
            DOF1TransformNode1.resize(edgeList.size());
        }

        ReadAccessor<Data<VecCoord> > statePos = this->m_mstate->read(sofa::core::vec_id::read_access::position) ;

        auto lengthList = sofa::helper::getWriteOnlyAccessor(this->d_lengthList);
        lengthList.clear();

        const unsigned int edgeListSize = this->d_edgeList.getValue().size();
        unsigned int nd0Id=0, nd1Id=0;


        for (unsigned int i = 0; i < edgeListSize; i++)
        {

            //// Copute the distance between beam nodes to obtain a first value for the beam length
            // Indeed the function "computeActualLength" needs a first approximation of this length

            this->getNodeIndices(i, nd0Id, nd1Id);

            if(this->d_dofsAndBeamsAligned.getValue())
            {
                //  when the dof and the beams are aligned it means that we can take the distance between nodes:
                Vec3 beam_segment = statePos[nd1Id].getCenter() - statePos[nd0Id].getCenter();
                lengthList.push_back(beam_segment.norm());
            }
            else
            {
                // if not aligned, we need to account for the transformation between node and dof:
                // this transforamtion is given by global_H_local0 for node 0 (and dof0)
                // and global_H_local1 for node 1 (and dof1)
                Transform global_H_local0, global_H_local1;
                this->computeTransform(i, nd0Id, nd1Id, global_H_local0, global_H_local1, statePos.ref());
                Vec3 beam_segment = global_H_local1.getOrigin() -  global_H_local0.getOrigin();
                lengthList.push_back(beam_segment.norm());

            }

            
            // finding the real length can not be done in one step
            // we do it in 3 iterations
            for (unsigned it=0; it<3; it++){
                 // now that we have an approximation of the length, we can estimate the position of the spline control point;
                 Vec3 P0,P1,P2,P3;
                this->getSplinePoints(i, statePos.ref(), P0, P1, P2, P3);

                // and we can compute more precisely the length
                Real ActualLength;
                this->computeActualLength(ActualLength,P0,P1,P2,P3);
                lengthList[i]=ActualLength;
            }


        }
    }

    if(!verifyTopology())
    {
        d_componentState.setValue(ComponentState::Invalid);
        return ;
    }
}

template<class DataTypes>
void BeamInterpolation<DataTypes>::reinit()
{
    init(); 
    bwdInit();
}

template<class DataTypes>
void BeamInterpolation<DataTypes>::storeResetState()
{
    if(d_componentState.getValue()==ComponentState::Invalid)
        return ;

    updateInterpolation();
}

template<class DataTypes>
void BeamInterpolation<DataTypes>::reset()
{
    if(d_componentState.getValue()==ComponentState::Invalid)
        return ;

    bwdInit();
}

template<class DataTypes>
bool BeamInterpolation<DataTypes>::interpolationIsAlreadyInitialized()
{
    if (this->d_edgeList.getValue().size() == 0)
        return false;

    const unsigned int nbEdges = this->d_edgeList.getValue().size();

    if (this->d_DOF0TransformNode0.getValue().size() != nbEdges)
        return false;

    if (this->d_DOF1TransformNode1.getValue().size() != nbEdges)
        return false;

    if (this->d_lengthList.getValue().size() != nbEdges)
        return false;

    return true;
}


/// verify that we got a non-null pointer on the topology
/// verify that the m_edgeList do not contain non existing edges
template<class DataTypes>
bool BeamInterpolation<DataTypes>::verifyTopology()
{
    //TODO(dmarchal) This contains "code" specific slang that cannot be understood by user.
    dmsg_info() << "The vector _topologyEdges is now set with " << this->m_topologyEdges->size() << " edges" ;


    const VecElementID &edgeList = this->d_edgeList.getValue();
    for (unsigned int j = 0; j < edgeList.size(); j++)
    {
        if(edgeList[j] > this->m_topologyEdges->size())
        {
            msg_error() << "The provided edge index '" << edgeList[j] << "'is larger than '"
                        << this->m_topologyEdges->size() << "' the amount of edges in the topology. " ;
            return false;
        }
    }

    return true;
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::clear()
{
    auto edgeList = sofa::helper::getWriteOnlyAccessor(this->d_edgeList);
    auto lengthList = sofa::helper::getWriteOnlyAccessor(this->d_lengthList);
    auto DOF0TransformNode0 = sofa::helper::getWriteOnlyAccessor(this->d_DOF0TransformNode0);
    auto DOF1TransformNode1 = sofa::helper::getWriteOnlyAccessor(this->d_DOF1TransformNode1);
    auto curvAbsList = sofa::helper::getWriteOnlyAccessor(this->d_curvAbsList);

    edgeList.clear();
    lengthList.clear();
    DOF0TransformNode0.clear();
    DOF1TransformNode1.clear();
    curvAbsList.clear();
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1, const Real &angle)
{
    auto edgeList = sofa::helper::getWriteOnlyAccessor(this->d_edgeList);
    auto lengthList = sofa::helper::getWriteOnlyAccessor(this->d_lengthList);
    auto DOF0TransformNode0 = sofa::helper::getWriteOnlyAccessor(this->d_DOF0TransformNode0);
    auto DOF1TransformNode1 = sofa::helper::getWriteOnlyAccessor(this->d_DOF1TransformNode1);
    auto curvAbsList = sofa::helper::getWriteOnlyAccessor(this->d_curvAbsList);

    curvAbsList.push_back(Vec2(x0, x1));

    edgeList.push_back(eID);
    lengthList.push_back(length);

    Quat QuatX ;
    QuatX.axisToQuat(Vec3(1,0,0), angle);
    QuatX.normalize();

    // as an angle is set between DOFs and Beam, they are no more aligned
    d_dofsAndBeamsAligned.setValue(false);
    DOF0TransformNode0.push_back(Transform(Vec3(0, 0, 0), QuatX ));
    DOF1TransformNode1.push_back(Transform(Vec3(0, 0, 0), QuatX ));
}


template <class DataTypes>
void BeamInterpolation<DataTypes>::getSamplingParameters(type::vector<Real>& /*xP_noticeable*/, type::vector< int>& /*nbP_density*/)
{
    msg_error()<<"getSamplingParameters is not implemented when _restShape== nullptr : TODO !! ";
}

template <class DataTypes>
typename BeamInterpolation<DataTypes>::Real BeamInterpolation<DataTypes>::getRestTotalLength()
{
    Real le(0.0);
    const type::vector< double > &lengthList = this->d_lengthList.getValue();

    for (unsigned int i = 0; i < lengthList.size(); i++)
        le += lengthList[i];

    return le;
}

template <class DataTypes>
void BeamInterpolation<DataTypes>::getCollisionSampling(Real &dx, const Real& /*x_localcurv_abs*/)
{
    unsigned int numLines = 30;
    dx = getRestTotalLength()/numLines;
}

template <class DataTypes>
void BeamInterpolation<DataTypes>::getNumberOfCollisionSegment(Real &dx, unsigned int &numLines)
{
    numLines = 30;
    dx = getRestTotalLength()/numLines;
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getInterpolationParameters(sofa::Index beamId, Real& _L, Real& _A, Real& _Iy,
    Real& _Iz, Real& _Asy, Real& _Asz, Real& _J)
{
    /// get the length of the beam:
    _L = this->d_lengthList.getValue()[beamId];

    const auto& section = getBeamSection(beamId);
    _A = section._A;
    _Iy = section._Iy;
    _Iz = section._Iz;
    _Asy = section._Asy;
    _Asz = section._Asz;
    _J = section._J;
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getMechanicalParameters(sofa::Index beamId, Real& youngModulus, Real& cPoisson, Real& massDensity)
{
    SOFA_UNUSED(massDensity);

    const auto& defaultYoungModuli = d_defaultYoungModulus.getValue();
    youngModulus = (beamId < static_cast<sofa::Index>(defaultYoungModuli.size()))? defaultYoungModuli[beamId]: m_defaultYoungModulus;

    const auto& poissonRatios = d_poissonRatio.getValue();
    cPoisson = (beamId < static_cast<sofa::Index>(defaultYoungModuli.size()))? poissonRatios[beamId]: m_defaultPoissonRatio;

    //TODO: massDensity??
}



template <class DataTypes>
void BeamInterpolation<DataTypes>::setTransformBetweenDofAndNode(int beam, const Transform &DOF_H_Node, unsigned int zeroORone )
{
    if (beam > int(this->d_DOF0TransformNode0.getValue().size()-1) || beam > int(this->d_DOF1TransformNode1.getValue().size()-1))
    {
        msg_error()<<"WARNING setTransformBetweenDofAndNode on non existing beam";
        return;
    }

    if (!zeroORone)
    {
        auto DOF0TransformNode0 = sofa::helper::getWriteOnlyAccessor(this->d_DOF0TransformNode0);
        DOF0TransformNode0[beam] = DOF_H_Node;
    }
    else
    {
        auto DOF1TransformNode1 = sofa::helper::getWriteOnlyAccessor(this->d_DOF1TransformNode1);
        DOF1TransformNode1[beam] = DOF_H_Node;
    }
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest)
{
    if(d_straight.getValue())
    {
        /// the beam is straight: local is in the middle of local0 and local1
        /// the transformation between local0 and local1 is provided by the length of the beam
        local_H_local0_rest.set(-Vec3(this->d_lengthList.getValue()[edgeInList]/2,0,0), Quat());
        local_H_local1_rest.set(Vec3(this->d_lengthList.getValue()[edgeInList]/2,0,0), Quat());
    }
    else
    {
        MechanicalState<DataTypes>* state = dynamic_cast<MechanicalState<DataTypes>*>(this->getContext()->getMechanicalState());

        if(state)
        {
            unsigned int node0Id, node1Id;
            this->getNodeIndices(edgeInList,node0Id,node1Id);

            Coord global_0 = state->read(sofa::core::vec_id::read_access::restPosition)->getValue()[node0Id];
            Coord global_1 = state->read(sofa::core::vec_id::read_access::restPosition)->getValue()[node1Id];

            Transform global_H_DOF0 = Transform(global_0.getCenter(),global_0.getOrientation());
            Transform global_H_DOF1 = Transform(global_1.getCenter(),global_1.getOrientation());

            Transform DOF0_H_local0, DOF1_H_local1;
            this->getDOFtoLocalTransform(edgeInList, DOF0_H_local0,  DOF1_H_local1);

            Transform global_H_local_0 = global_H_DOF0*DOF0_H_local0;
            Transform global_H_local_1 = global_H_DOF1*DOF1_H_local1;


            Transform global_H_local_middle;
            Real baryX = 0.5;
            Real L = this->getLength(edgeInList);

            this->InterpolateTransformUsingSpline(global_H_local_middle, baryX,  global_H_local_0, global_H_local_1, L);

            local_H_local0_rest = global_H_local_middle.inversed() * global_H_local_0;
            local_H_local1_rest = global_H_local_middle.inversed() * global_H_local_1;
        }
        else
            msg_error() <<"This component needs a context mechanical state if the 'straight' parameter is set to false." ;
    }
}


template<class DataTypes>
void BeamInterpolation<DataTypes>::getTangentUsingSplinePoints(unsigned int edgeInList, const Real& baryCoord,
                                                               const sofa::core::ConstVecCoordId &vecXId, Vec3& t )
{

    const VectorVec3& splinePos = this->m_StateNodes->read(vecXId)->getValue();
    Vec3 P0,P1,P2,P3;

    P0=splinePos[edgeInList*4];
    P1=splinePos[edgeInList*4+1];
    P2=splinePos[edgeInList*4+2];
    P3=splinePos[edgeInList*4+3];

    Real bx=baryCoord;

    t  = P0*-3*(1-bx)*(1-bx) + P1*(3-12*bx+ 9*bx*bx) + P2*(6*bx-9*bx*bx) + P3*3*bx*bx;
}




template<class DataTypes>
void BeamInterpolation<DataTypes>::updateInterpolation(){
    /// this function allow to use BeamInterpolation using data (TODO: see if it could fall into engine concept)
    /// inputs :
    /// -   1.Data<helper::OptionsGroup > d_vecID;
    ///         VecID => (1) "current" Pos, Vel    (2) "free" PosFree, VelFree   (3) "rest" PosRest, V=0
    ///
    /// -   2. Data< type::vector< Vec2 > > d_InterpolationInputs;
    ///         Vector of 2-tuples (indice of the beam   ,   barycentric between 0 and 1)
    /// ouptus:
    /// -   1. d_InterpolatedPos  => result interpolate the position (6 D0Fs type rigid: position/quaterion)
    /// -   2. d_InterpolatedVel  => result interpolate the velocities (linear velocity / angular velocity)

    if (BEAMADAPTER_WITH_VERIFICATION)
        dmsg_info() <<"entering updateInterpolation" ;

    const type::vector< Vec2 > &interpolationInputs = d_InterpolationInputs.getValue();
    unsigned int numInterpolatedPositions = interpolationInputs.size();

    auto interpolatedPos = sofa::helper::getWriteOnlyAccessor(d_InterpolatedPos);
    auto interpolatedVel = sofa::helper::getWriteOnlyAccessor(d_InterpolatedVel);

    interpolatedPos.resize(numInterpolatedPositions);
    interpolatedVel.resize(numInterpolatedPositions);

    ///////// select the input position / velocity DATA /////
    if (BEAMADAPTER_WITH_VERIFICATION)
        dmsg_info() << "select the input position  "<<"  selected Item : "<<d_vecID.getValue().getSelectedItem() ;

    const Data< VecCoord > *x = nullptr;

    ///const Data< VecDeriv > *v;
    bool computeVel = true;
    if(d_vecID.getValue().getSelectedItem() == "current")
    {
        dmsg_info() <<" position " << msgendl
                   << "      ="<< this->m_mstate->read( sofa::core::vec_id::read_access::position )->getValue( ) ;
        x=this->m_mstate->read( sofa::core::vec_id::read_access::position );
    }
    else if(d_vecID.getValue().getSelectedItem() == "free")
    {
        x=this->m_mstate->read( sofa::core::vec_id::read_access::freePosition ) ;
    }
    else /// rest position
    {
        x=this->m_mstate->read( sofa::core::vec_id::read_access::restPosition ) ;
        computeVel = false;
    }

    if (BEAMADAPTER_WITH_VERIFICATION)
        dmsg_info() << "selected Item : "<<d_vecID.getValue().getSelectedItem()<<" compute Vel"<< computeVel ;

    for (unsigned int i=0; i<numInterpolatedPositions; i++)
    {
        unsigned int numBeam =  (unsigned int) interpolationInputs[i][0];
        Real baryCoord= interpolationInputs[i][1];

        /// convert position (RigidTypes) in Transform
        Transform global_H_local0, global_H_local1;
        this->computeTransform(numBeam, global_H_local0, global_H_local1, x->getValue() );

        /// compute the length of the spline
        Vec3 P0,P1,P2,P3;
        this->getSplinePoints(numBeam, x->getValue() , P0,  P1, P2, P3);
        Real length;
        this->computeActualLength(length,P0,P1,P2,P3);

        /// get the result of the transform:
        Transform global_H_localResult;

        if(computeVel)
        {
            this->InterpolateTransformUsingSpline(global_H_localResult,baryCoord,global_H_local0,global_H_local1,length);
        }
//        else
//        {
//            InterpolateTransformUsingSpline(global_H_localResult,baryCoord,global_H_local0,global_H_local1,length);
//        }

        /// assign the result in the output data
        interpolatedPos[i].getCenter() = global_H_localResult.getOrigin();
        interpolatedPos[i].getOrientation() = global_H_localResult.getOrientation();
    }
}


} // namespace beamadapter
