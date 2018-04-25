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
// C++ Implementation : AdaptiveBeamMapping
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_MAPPING_BEAMLENGTHMAPPING_INL
#define SOFA_COMPONENT_MAPPING_BEAMLENGTHMAPPING_INL

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include "BeamLengthMapping.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <string>
#include <sofa/core/Mapping.inl>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/helper/AdvancedTimer.h>
#include <iomanip>


namespace sofa
{

namespace component
{

namespace mapping
{

namespace _beamlengthmapping_
{

using namespace sofa::defaulttype;
using sofa::core::State;
using helper::ReadAccessor;
using helper::WriteAccessor;
using sofa::core::ConstVecCoordId;
using sofa::helper::AdvancedTimer;
using sofa::core::MultiVecCoordId;
using sofa::core::VecCoordId;
using sofa::core::ConstMultiVecCoordId;
using core::MechanicalParams;

template <class TIn, class TOut>
BeamLengthMapping<TIn,TOut>::BeamLengthMapping(State< In >* from, State< Out >* to,
                                                   BeamInterpolation< TIn >* interpolation)
    : Inherit(from, to)
    , l_adaptativebeamInterpolation(initLink("interpolation", "Path to the Interpolation component on scene"), interpolation)
{


}


template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::init()
{

    BaseContext *context= dynamic_cast<BaseContext *>(this->getContext());

   // BInterpolation* i_test = context->get<BInterpolation>();

    if (!l_adaptativebeamInterpolation)
        l_adaptativebeamInterpolation.set(context->get<BInterpolation>());

    if (!l_adaptativebeamInterpolation)
        msg_error() <<"No Beam Interpolation found, the component can not work.";

}

template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::bwdInit()
{
    // todo ?


}

template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::reinit()
{
    init();

}


template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::reset()
{
    reinit();
}



template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::apply(const MechanicalParams* mparams, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    AdvancedTimer::stepBegin("AdaptiveBeamMappingApply");
    VecCoord& out = *dOut.beginEdit();
    const InVecCoord& in = dIn.getValue();



    unsigned int s = l_adaptativebeamInterpolation->getNumBeams();

    out.resize(s);

    for (unsigned int i=0; i<s; i++)
    {
        Vec<3, InReal> P0,P1,P2,P3;
        l_adaptativebeamInterpolation->getSplinePoints(i, in , P0,  P1, P2, P3);

        InReal length;
        l_adaptativebeamInterpolation->computeActualLength(length,P0,P1,P2,P3);

        out[i][0]=length;

    }


    dOut.endEdit();
    AdvancedTimer::stepEnd("AdaptiveBeamMappingApply");
}


template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::applyJ(const core::MechanicalParams* mparams, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    SOFA_UNUSED(mparams);

    AdvancedTimer::stepBegin("AdaptiveBeamMappingApplyJ");
    VecDeriv& out = *dOut.beginEdit();
    const InVecDeriv& in= dIn.getValue();

    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(VecCoordId::position());
    const InVecCoord& x_in = dataInX.getValue();

    unsigned int s = l_adaptativebeamInterpolation->getNumBeams();
    out.resize(s);
    for (unsigned int i=0; i<s; i++)
    {
        //1. get the indices of the Dofs of the beam a
        unsigned int IdxNode0, IdxNode1;
        l_adaptativebeamInterpolation->getNodeIndices(i,IdxNode0,IdxNode1);

        //2. get the velocity (or variation of position) of the Dofs
        SpatialVector vDOF0, vDOF1;
        vDOF0.setLinearVelocity (In::getDPos(in[IdxNode0]));
        vDOF0.setAngularVelocity(In::getDRot(in[IdxNode0]));
        vDOF1.setLinearVelocity (In::getDPos(in[IdxNode1]));
        vDOF1.setAngularVelocity(In::getDRot(in[IdxNode1]));


        // 3. get the transform to DOF in global frame from local frame
        Transform DOF0Global_H_local0, DOF1Global_H_local1;
        l_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(i, DOF0Global_H_local0, DOF1Global_H_local1, x_in);

        // 4. project the velocities in local frame:
        SpatialVector v_local0, v_local1;
        v_local0 = DOF0Global_H_local0.inversed()*vDOF0;
        v_local1 = DOF1Global_H_local1.inversed()*vDOF1;


        // 5. Computes the local velocities of the 4 points of the spline (equivalent to applyJ in rigid Mapping)
        Real L = l_adaptativebeamInterpolation->getLength(i);
        Vec3 lever(L/3.0,0,0);
        Vec3 v0, v1, v2, v3;
        v0 = v_local0.getLinearVelocity();
        v1 = v0 - lever.cross(v_local0.getAngularVelocity());
        v3 = v_local1.getLinearVelocity();
        v2 = v3 + lever.cross(v_local1.getAngularVelocity());

        // 6. rotate back velocities to the global frame
        Vec3 V0, V1, V2, V3;
        V0= DOF0Global_H_local0.getOrientation().rotate(v0) ;
        V1= DOF0Global_H_local0.getOrientation().rotate(v1) ;
        V2= DOF1Global_H_local1.getOrientation().rotate(v2) ;
        V3= DOF1Global_H_local1.getOrientation().rotate(v3) ;

        // 7. Compute the "velocity" of the length
        Vec<3, InReal> P0,P1,P2,P3;
        l_adaptativebeamInterpolation->getSplinePoints(i, x_in , P0,  P1, P2, P3);


        Real dlength;
        computeJSpline(dlength, P0,P1,P2,P3,V0,V1,V2,V3);

        out[i][0]=dlength;


    }


    dOut.endEdit();
    AdvancedTimer::stepEnd("AdaptiveBeamMappingApplyJ");
}


template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::applyJT(const core::MechanicalParams* mparams, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{
    SOFA_UNUSED(mparams);

    AdvancedTimer::stepBegin("AdaptiveBeamMappingMechanicalApplyJT");
    InVecDeriv& out = *dOut.beginEdit();
    const VecDeriv& in= dIn.getValue();

    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(VecCoordId::position());
    const InVecCoord& x_in = dataInX.getValue();

    unsigned int s = l_adaptativebeamInterpolation->getNumBeams();
    for (unsigned int i=0; i<s; i++)
    {
        //1. get the indices of the Dofs of the beam a
        unsigned int IdxNode0, IdxNode1;
        l_adaptativebeamInterpolation->getNodeIndices(i,IdxNode0,IdxNode1);

        //2. get the force on the mapped dof
        Real mappedF = in[i][0];

        //3. get the spline points
        Vec<3, InReal> P0,P1,P2,P3;
        l_adaptativebeamInterpolation->getSplinePoints(i, x_in , P0,  P1, P2, P3);

        //4. compute the equivalent forces on the spline control point (apply Jt on spline map)
        Vec3 F0, F1, F2, F3;

        computeJtSpline(mappedF, P0,P1,P2,P3, F0, F1, F2, F3);

        // 5. apply the forces to the nodes ot the beams (equivalent to applyJt on spline control point)

        Transform DOF0Global_H_local0, DOF1Global_H_local1;
        l_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(i, DOF0Global_H_local0, DOF1Global_H_local1, x_in);

        // rotate back the force to the local frame
        SpatialVector f0, f1,f2,f3;
        f0.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F0) );
        f1.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F1) );
        f2.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F2) );
        f3.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F3) );

        // computes the torque created on DOF0 and DOF1 by f1 and f2
        Real L = l_adaptativebeamInterpolation->getLength(i);
        Vec3 lever(L/3.0,0,0);
        f1.setTorque(lever.cross(f1.getForce()));
        f2.setTorque(-lever.cross(f2.getForce()));

        // back to the DOF0 and DOF1 frame:
        SpatialVector FNode0output, FNode1output;
        FNode0output = DOF0Global_H_local0 * (f0+f1);
        FNode1output = DOF1Global_H_local1 * (f2+f3);



        //2. put the result in out vector computes the equivalent forces on nodes + rotate to Global Frame from DOF frame
        In::setDPos(out[IdxNode0], In::getDPos(out[IdxNode0]) + FNode0output.getForce() ); // add  trans forces to out
        In::setDPos(out[IdxNode1], In::getDPos(out[IdxNode1]) + FNode1output.getForce());

        In::setDRot(out[IdxNode0], In::getDRot(out[IdxNode0]) + FNode0output.getTorque()); // add torques to out
        In::setDRot(out[IdxNode1], In::getDRot(out[IdxNode1]) + FNode1output.getTorque());

    }


    dOut.endEdit();
    AdvancedTimer::stepEnd("AdaptiveBeamMappingMechanicalApplyJT");
}


/// AdaptiveBeamMapping::applyJT(InMatrixDeriv& out, const OutMatrixDeriv& in)
/// this function propagate the constraint through the Adaptive Beam mapping :
/// if one constraint along (vector n) with a value (v) is applied on the childModel (like collision model)
/// then this constraint is transformed by (Jt.n) with value (v) for the rigid model
/// note : the value v is not propagated through the mapping
template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::applyJT(const core::ConstraintParams* cparams, Data<InMatrixDeriv>& dOut, const Data<OutMatrixDeriv>& dIn)
{
    SOFA_UNUSED(cparams);

    AdvancedTimer::stepBegin("AdaptiveBeamMappingConstrainApplyJT");

    InMatrixDeriv& out = *dOut.beginEdit();
    const OutMatrixDeriv& in = dIn.getValue();
    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(ConstVecCoordId::position());
    const InVecCoord& x_in = dataInX.getValue();


    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();
    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();
        for (typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin(); colIt != colItEnd; ++colIt)
        {
            typename In::MatrixDeriv::RowIterator o = out.writeLine(rowIt.index());
            unsigned int indexIn = colIt.index();
            const Deriv data = colIt.val();


            if (indexIn >= l_adaptativebeamInterpolation->getNumBeams())
            {
                msg_warning() <<"Wrong index in VecConst in  the beam does not exist";
                break;
            }


            //1. get the indices of the Dofs of the beam a
            unsigned int IdxNode0, IdxNode1;
            l_adaptativebeamInterpolation->getNodeIndices(indexIn,IdxNode0,IdxNode1);

            //2. get the force on the mapped dof
            Real mappedF = data[0];

            //3. get the spline points
            Vec<3, InReal> P0,P1,P2,P3;
            l_adaptativebeamInterpolation->getSplinePoints(indexIn, x_in , P0,  P1, P2, P3);

            //4. compute the equivalent forces on the spline control point (apply Jt on spline map)
            Vec3 F0, F1, F2, F3;

            computeJtSpline(mappedF, P0,P1,P2,P3, F0, F1, F2, F3);

            // 5. apply the forces to the nodes ot the beams (equivalent to applyJt on spline control point)

            Transform DOF0Global_H_local0, DOF1Global_H_local1;
            l_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(indexIn, DOF0Global_H_local0, DOF1Global_H_local1, x_in);

            // rotate back the force to the local frame
            SpatialVector f0, f1,f2,f3;
            f0.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F0) );
            f1.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F1) );
            f2.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F2) );
            f3.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F3) );

            // computes the torque created on DOF0 and DOF1 by f1 and f2
            Real L = l_adaptativebeamInterpolation->getLength(indexIn);
            Vec3 lever(L/3,0,0);
            f1.setTorque(lever.cross(f1.getForce()));
            f2.setTorque(-lever.cross(f2.getForce()));

            // back to the DOF0 and DOF1 frame:
            SpatialVector FNode0output, FNode1output;
            FNode0output = DOF0Global_H_local0 * (f0+f1);
            FNode1output = DOF1Global_H_local1 * (f2+f3);

            // Compute the mapped Constraint on the beam nodes
            InDeriv direction0;
            In::setDPos(direction0,FNode0output.getForce());
            In::setDRot(direction0,FNode0output.getTorque());
            InDeriv direction1;
            In::setDPos(direction1,FNode1output.getForce());
            In::setDRot(direction1,FNode1output.getTorque());

            o.addCol(IdxNode0, direction0);
            o.addCol(IdxNode1, direction1);


        }
    }


    dOut.endEdit();
    AdvancedTimer::stepEnd("AdaptiveBeamMappingConstraintApplyJT");
}




template <class TIn, class TOut>
void BeamLengthMapping<TIn, TOut>::computeJSpline(Real &dlength, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3& P3,
                                                const Vec3& dP0, const Vec3& dP1, const Vec3& dP2, const Vec3& dP3)
{
    /// the computation of the derivation of the
    /// integral Int[0,1]_||dP(x)||_ dx = length
    ///  done using Gauss Points

    /// definition of the Gauss points

    static Real A = 2*sqrt(6.0/5.0);
    static Real x1 = -sqrt((3.0 - A)/7.0 )/2.0+ 0.5;
    static Real x2 = sqrt((3.0 - A) /7.0 )/2.0+ 0.5;
    static Real x3 = -sqrt((3.0 + A)/7.0 )/2.0+ 0.5;
    static Real x4 = sqrt((3.0 + A) /7.0 )/2.0+ 0.5;



    Vec3 IP1, IP2, IP3, IP4, dIP1, dIP2, dIP3, dIP4;

    IP1 = P0*(-3*(1-x1)*(1-x1)) + P1*(3-12*x1+9*x1*x1) + P2*(6*x1-9*x1*x1) + P3*(3*x1*x1);
    IP2 = P0*(-3*(1-x2)*(1-x2)) + P1*(3-12*x2+9*x2*x2) + P2*(6*x2-9*x2*x2) + P3*(3*x2*x2);
    IP3 = P0*(-3*(1-x3)*(1-x3)) + P1*(3-12*x3+9*x3*x3) + P2*(6*x3-9*x3*x3) + P3*(3*x3*x3);
    IP4 = P0*(-3*(1-x4)*(1-x4)) + P1*(3-12*x4+9*x4*x4) + P2*(6*x4-9*x4*x4) + P3*(3*x4*x4);

    dIP1 = dP0*(-3*(1-x1)*(1-x1)) + dP1*(3-12*x1+9*x1*x1) + dP2*(6*x1-9*x1*x1) + dP3*(3*x1*x1);
    dIP2 = dP0*(-3*(1-x2)*(1-x2)) + dP1*(3-12*x2+9*x2*x2) + dP2*(6*x2-9*x2*x2) + dP3*(3*x2*x2);
    dIP3 = dP0*(-3*(1-x3)*(1-x3)) + dP1*(3-12*x3+9*x3*x3) + dP2*(6*x3-9*x3*x3) + dP3*(3*x3*x3);
    dIP4 = dP0*(-3*(1-x4)*(1-x4)) + dP1*(3-12*x4+9*x4*x4) + dP2*(6*x4-9*x4*x4) + dP3*(3*x4*x4);

    /// formula with 4 Gauss Points
    static  Real B= sqrt(30.0);
    dlength = ((18.0 + B) /72.0 )*(dot(IP1,dIP1))/IP1.norm() + ((18.0 + B) /72.0 )*(dot(IP2,dIP2))/IP2.norm() + ((18.0 - B) /72.0 )*(dot(IP3,dIP3))/IP3.norm() + ((18.0 - B) /72.0 )*(dot(IP4,dIP4))/IP4.norm();

}



template <class TIn, class TOut>
void BeamLengthMapping<TIn, TOut>::computeJtSpline(const Real &f_input, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3& P3,
                                                 Vec3& F0,  Vec3& F1,  Vec3& F2, Vec3& F3)
{
    /// In computeDLength, corresponds to a part of the applicaiton of the jacobian
    /// Here, we compute the transposed of the jacobian
    /// formula with 4 Gauss Points
    /// definition of the Gauss points

    static Real A = 2*sqrt(6.0/5.0);
    static Real x1 = -sqrt((3.0 - A)/7.0 )/2.0+ 0.5;
    static Real x2 = sqrt((3.0 - A) /7.0 )/2.0+ 0.5;
    static Real x3 = -sqrt((3.0 + A)/7.0 )/2.0+ 0.5;
    static Real x4 = sqrt((3.0 + A) /7.0 )/2.0+ 0.5;


    static Real a10= (-3*(1-x1)*(1-x1));
    static Real a11=(3-12*x1+9*x1*x1);
    static Real a12=(6*x1-9*x1*x1);
    static Real a13= (3*x1*x1);
    static Real a20=(-3*(1-x2)*(1-x2));
    static Real a21=(3-12*x2+9*x2*x2);
    static Real a22=(6*x2-9*x2*x2);
    static Real a23=(3*x2*x2);
    static Real a30=(-3*(1-x3)*(1-x3));
    static Real a31= (3-12*x3+9*x3*x3);
    static Real a32=(6*x3-9*x3*x3);
    static Real a33=(3*x3*x3);
    static Real a40= (-3*(1-x4)*(1-x4));
    static Real a41=(3-12*x4+9*x4*x4);
    static Real a42=(6*x4-9*x4*x4);
    static Real a43=(3*x4*x4);


    Vec3 IP1, IP2, IP3, IP4;

    IP1 = P0*a10 + P1*a11 + P2*a12 + P3*a13;
    IP2 = P0*a20 + P1*a21 + P2*a22 + P3*a23;
    IP3 = P0*a30 + P1*a31 + P2*a32 + P3*a33;
    IP4 = P0*a40 + P1*a41 + P2*a42 + P3*a43;

    Vec3 F_IP1, F_IP2, F_IP3, F_IP4;
    Real B= sqrt(30.0);
    F_IP1 = IP1*f_input*(18.0 + B)/(72.0*IP1.norm());
    F_IP2 = IP2*f_input*(18.0 + B)/(72.0*IP2.norm());
    F_IP3 = IP3*f_input*(18.0 - B)/(72.0*IP3.norm());
    F_IP4 = IP4*f_input*(18.0 - B)/(72.0*IP4.norm());

    F0 = F_IP1 * a10  + F_IP2 * a20  + F_IP3 * a30  + F_IP4 * a40 ;
    F1 = F_IP1 * a11  + F_IP2 * a21  + F_IP3 * a31  + F_IP4 * a41 ;
    F2 = F_IP1 * a12  + F_IP2 * a22  + F_IP3 * a32  + F_IP4 * a42 ;
    F3 = F_IP1 * a13  + F_IP2 * a23  + F_IP3 * a33  + F_IP4 * a43 ;

}


template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::draw(const VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMappings())
        return;
}


#ifndef SOFA_FLOAT
template class SOFA_BEAMADAPTER_API BeamLengthMapping<Rigid3dTypes, Vec1dTypes   >;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_BEAMADAPTER_API BeamLengthMapping< Rigid3fTypes, Vec1fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_BEAMADAPTER_API BeamLengthMapping< Rigid3dTypes, Vec1fTypes >;
template class SOFA_BEAMADAPTER_API BeamLengthMapping< Rigid3fTypes, Vec1dTypes >;
#endif
#endif


} /// namespace _beamlengthmapping_

} /// namespace mapping

} /// namespace component

} /// namespace sofa

#endif  /* SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_INL */
