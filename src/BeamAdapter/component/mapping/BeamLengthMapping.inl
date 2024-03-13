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
#include <BeamAdapter/component/mapping/BeamLengthMapping.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <string>
#include <sofa/core/Mapping.inl>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
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
using sofa::core::MultiVecCoordId;
using sofa::core::VecCoordId;
using sofa::core::VecDerivId;
using sofa::core::ConstMultiVecCoordId;
using core::MechanicalParams;

template <class TIn, class TOut>
BeamLengthMapping<TIn,TOut>::BeamLengthMapping(State< In >* from, State< Out >* to,
                                                   BeamInterpolation< TIn >* interpolation)
    : Inherit(from, to)
    , l_adaptativebeamInterpolation(initLink("interpolation", "Path to the Interpolation component on scene"), interpolation)
    , d_geometricStiffness(initData(&d_geometricStiffness, 2u, "geometricStiffness", "0 -> no GS, 1 -> exact GS, 2 -> stabilized GS (default)"))

{


}


template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::init()
{
    BaseContext *context= dynamic_cast<BaseContext *>(this->getContext());

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
    SOFA_UNUSED(mparams);
    SCOPED_TIMER("AdaptiveBeamMappingApply");

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
}


template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::applyJ(const core::MechanicalParams* mparams, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    SOFA_UNUSED(mparams);

    SCOPED_TIMER("AdaptiveBeamMappingApplyJ");

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
}


template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::applyJT(const core::MechanicalParams* mparams, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{
    SOFA_UNUSED(mparams);

    SCOPED_TIMER("AdaptiveBeamMappingMechanicalApplyJT");
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

        /* debug
        static int compteur=0;
        if(compteur==2 || compteur==6)
        {
            std::cout<<" comntribution spline freeze"<<std::endl;
            F0 = F0_buf;
            F1 = F1_buf;
            F2 = F2_buf;
            F3 = F3_buf;
        }
        else
        {
            F0_buf = F0;
            F1_buf = F1;
            F2_buf = F2;
            F3_buf = F3;
        }
        compteur++;


        std::cout<<" _F0  = "<<F0<<std::endl;
        std::cout<<" _F1  = "<<F1<<std::endl;
        std::cout<<" _F2  = "<<F2<<std::endl;
        std::cout<<" _F3  = "<<F3<<std::endl;

        std::cout<<" _F0+ _F1 = "<<F0+F1<<std::endl;
        std::cout<<" _F0+ _F1 = "<<F2+F3<<std::endl;

        std::cout<<std::setprecision(20)<<"P0 = ["<<P0<<"];"<<std::endl;
        std::cout<<"P1 = ["<<P1<<"];"<<std::endl;
        std::cout<<"P2 = ["<<P2<<"];"<<std::endl;
        std::cout<<"P3 = ["<<P3<<"];"<<std::endl;

         */
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

        /* debug
        std::cout<<" FNode0output = "<<FNode0output<<"  f0+f1 = "<<(f0+f1)<<std::endl;
        std::cout<<" FNode1output = "<<FNode0output<<"  f2+f3 = "<<(f2+f3)<<std::endl;
        */





        //2. put the result in out vector computes the equivalent forces on nodes + rotate to Global Frame from DOF frame
        In::setDPos(out[IdxNode0], In::getDPos(out[IdxNode0]) + FNode0output.getForce() ); // add  trans forces to out
        In::setDPos(out[IdxNode1], In::getDPos(out[IdxNode1]) + FNode1output.getForce());

        In::setDRot(out[IdxNode0], In::getDRot(out[IdxNode0]) + FNode0output.getTorque()); // add torques to out
        In::setDRot(out[IdxNode1], In::getDRot(out[IdxNode1]) + FNode1output.getTorque());

    }


    dOut.endEdit();
}


/// BeamLengthMapping::applyJT(InMatrixDeriv& out, const OutMatrixDeriv& in)
/// this function propagate the constraint through BeamLengthMapping :
template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::applyJT(const core::ConstraintParams* cparams, Data<InMatrixDeriv>& dOut, const Data<OutMatrixDeriv>& dIn)
{
    SOFA_UNUSED(cparams);

    SCOPED_TIMER("AdaptiveBeamMappingConstrainApplyJT");

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
}

/// BeamLengthMapping::applyDJT(MultiVecDerivId parentDfId, const ConstMultiVecDerivId childDfId)
/// this function computes the additional stiffness force created on the parents by the force on the child (due to mapping non-linearity)
template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::applyDJT(const MechanicalParams* mparams, core::MultiVecDerivId parentDfId, core::ConstMultiVecDerivId childDfId)
{
    const unsigned& geometricStiffness = d_geometricStiffness.getValue();
    if( !geometricStiffness ) return;

    const SReal kfactor = mparams->kFactor();

    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(VecCoordId::position());
    const InVecCoord& x_in = dataInX.getValue();

    const Data<InVecDeriv>& dataIndX = *this->getFromModel()->read(VecDerivId::dx());
    const InVecDeriv& parentDisplacement = dataIndX.getValue();

    helper::WriteAccessor<Data<InVecDeriv> > parentForce (*parentDfId[this->fromModel.get()].write());

    /*    const VecDeriv& childForce = this->getToModel()->readForces().ref();*/
    helper::ReadAccessor<Data<VecDeriv> > childForce( *childDfId[this->toModel.get()].read() );

    msg_info() << "********applyDJT*********\n childForce=" << childForce.ref() <<"\n*************";


    unsigned int s = l_adaptativebeamInterpolation->getNumBeams();
    for (unsigned int i=0; i<s; i++)
    {
        // force in compression (>0) can lead to negative eigen values in geometric stiffness
        // this results in a undefinite implicit matrix that causes instabilies
        // if stabilized GS (geometricStiffness==2) -> keep only force in extension
        if( childForce[i][0] < 0 || geometricStiffness==1 )
        {

            //1. get the indices of the Dofs of the beam a
            unsigned int IdxNode0, IdxNode1;
            l_adaptativebeamInterpolation->getNodeIndices(i,IdxNode0,IdxNode1);

            //2. get the force on the mapped dof
            Real childF = childForce[i][0];

            //3. get the spline points
            Vec<3, InReal> P0,P1,P2,P3;
            l_adaptativebeamInterpolation->getSplinePoints(i, x_in , P0,  P1, P2, P3);

            ////////////////////////////
            //4. compute the equivalent stiffness on the spline control point (apply DJt on spline map)
            Mat<4,4,Mat3> K;
            computeDJtSpline(childF, P0,P1,P2,P3, K);

            //- K must be transfered from the control points to the DOFs
            //  - get the transformation of the DOFs
            Transform DOF0Global_H_local0, DOF1Global_H_local1;
            l_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(i, DOF0Global_H_local0, DOF1Global_H_local1, x_in);

            // - create a matrix of the lever in the global frame
            Real L = l_adaptativebeamInterpolation->getLength(i);
            // - rotate the levers in the global frame
            Vec3 lev(-L/3.0,0.0,0.0);
            Vec3 Lev00_global = -DOF0Global_H_local0.getOrigin();
            Vec3 Lev01_global = DOF0Global_H_local0.getOrientation().rotate(lev) - DOF0Global_H_local0.getOrigin();
            lev[0]=L/3;
            Vec3 Lev12_global = DOF1Global_H_local1.getOrientation().rotate(lev) - DOF1Global_H_local1.getOrigin();
            Vec3 Lev13_global = -DOF1Global_H_local1.getOrigin();
            // create matrices:
            Mat3 Lev00_mat, Lev01_mat, Lev12_mat,Lev13_mat ;
            createCrossMatrix(Lev01_global, Lev01_mat);
            createCrossMatrix(Lev12_global, Lev12_mat);
            createCrossMatrix(Lev00_global, Lev00_mat);
            createCrossMatrix(Lev13_global, Lev13_mat);



            // Compute the  parent foces due to input parent displacement

            InDeriv Dx_DOF0 = parentDisplacement[IdxNode0];
            InDeriv Dx_DOF1 = parentDisplacement[IdxNode1];

            // displacement of the control points
            Vec3 dP[4];
            dP[0]= Dx_DOF0.getVCenter() + Lev00_mat*Dx_DOF0.getVOrientation();
            dP[1]= Dx_DOF0.getVCenter() + Lev01_mat*Dx_DOF0.getVOrientation();
            dP[2]= Dx_DOF1.getVCenter() + Lev12_mat*Dx_DOF1.getVOrientation();
            dP[3]= Dx_DOF1.getVCenter() + Lev13_mat*Dx_DOF1.getVOrientation();;

            /* debug
            for (unsigned int k=0; k<4; k++)
            {
                std::cout<<"dP_"<<k<<"="<<dP[k]<<std::endl;
            }
            */

            Vec3 dF[4];
            for (unsigned int k=0; k<4; k++)
            {
                dF[k].clear();
                for (unsigned int l=0; l<4; l++)
                {
                    dF[k] += K[k][l]*dP[l]*kfactor;
                }

                /* debug
                std::cout<<"dF_"<<k<<"="<<dF[k]<<std::endl;
                */
            }



            InDeriv F_dof0,F_dof1;
            F_dof0.clear();
            F_dof1.clear();

            /* DEBUG : A COMMENTER quand on veut freezer la contribution de spline dans ApplyJt et comparer numeriquement*/
            F_dof0.getVCenter()= dF[0]+dF[1];
            F_dof0.getVOrientation()=Lev00_mat.transposed()*dF[0] + Lev01_mat.transposed()*dF[1];
            F_dof1.getVCenter()=dF[2]+dF[3];
            F_dof1.getVOrientation()=Lev12_mat.transposed()*dF[2] + Lev13_mat.transposed()*dF[3];






             ////////////////////////////
            // 5. compute the equivalent stiffness on the rigid rotation of control point (apply DJt on rigid at fixed forces on spline)
            // -compute the equivalent forces on the spline control point (apply Jt on spline map)

            Vec3 F0, F1, F2, F3;
            computeJtSpline(childF, P0,P1,P2,P3, F0, F1, F2, F3);



            Mat3 F0_mat, F1_mat, F2_mat, F3_mat;
            createCrossMatrix(F0, F0_mat);
            createCrossMatrix(F1, F1_mat);
            createCrossMatrix(F2, F2_mat);
            createCrossMatrix(F3, F3_mat);


            Vec3 F0k = F0_mat*Lev00_mat*Dx_DOF0.getVOrientation()+ F1_mat*Lev01_mat*Dx_DOF0.getVOrientation();
            Vec3 F1k = F2_mat*Lev12_mat*Dx_DOF1.getVOrientation()+ F3_mat*Lev13_mat*Dx_DOF1.getVOrientation();




            /* debug
            std::cout<< "F0k ="<<F0k<< " rotDOF0 = "<<Dx_DOF0.getVOrientation()<<std::endl;
            std::cout<< "F1k ="<<F1k<< " rotDOF1 = "<<Dx_DOF1.getVOrientation()<<std::endl;
            */


            F_dof0.getVOrientation()-=F0k*kfactor;
            F_dof1.getVOrientation()-=F1k*kfactor;

            parentForce[IdxNode0]+=F_dof0;
            parentForce[IdxNode1]+=F_dof1;

        }
    }


}

template <class TIn, class TOut>
void BeamLengthMapping<TIn, TOut>::updateK(const core::MechanicalParams* mparams, core::ConstMultiVecDerivId childForceId )
{
    SOFA_UNUSED(mparams);
    const unsigned& geometricStiffness = d_geometricStiffness.getValue();
    if( !geometricStiffness ) { K_geom.resize(0,0); return; }
    //helper::ReadAccessor<Data<VecDeriv> > childForce( *childForceId[(const core::State<TOut>*)this->getToModels()[0]].read() );

    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(VecCoordId::position());
    const InVecCoord& x_in = dataInX.getValue();

    //const VecDeriv& childForce = this->getToModel()->readForces().ref();
    helper::ReadAccessor<Data<VecDeriv> > childForce( *childForceId[this->toModel.get()].read() );

    unsigned int s = l_adaptativebeamInterpolation->getNumBeams();

    msg_info() << "********updateK*********\n childForce=" << childForce.ref() <<"\n*************";

    K_geom.resize(Nin*x_in.size(), Nin*x_in.size());
    for (unsigned int i=0; i<s; i++)
    {
        // force in compression (>0) can lead to negative eigen values in geometric stiffness
        // this results in a undefinite implicit matrix that causes instabilies
        // if stabilized GS (geometricStiffness==2) -> keep only force in extension
        if( childForce[i][0] < 0 || geometricStiffness==1 )
        {
            //1. get the indices of the Dofs of the beam a
            unsigned int IdxNode[2];;
            l_adaptativebeamInterpolation->getNodeIndices(i,IdxNode[0],IdxNode[1]);

            //2. get the force on the mapped dof
            Real childF = childForce[i][0];

            //3. get the spline points
            Vec<3, InReal> P0,P1,P2,P3;
            l_adaptativebeamInterpolation->getSplinePoints(i, x_in , P0,  P1, P2, P3);

            ////////////////////////////
            //4. compute the equivalent stiffness on the spline control point (apply DJt on spline map)
            Mat<4,4,Mat3> K;
            computeDJtSpline(childF, P0,P1,P2,P3, K);

            //- K must be transfered from the control points to the DOFs
            //  - get the transformation of the DOFs
            Transform DOF0Global_H_local0, DOF1Global_H_local1;
            l_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(i, DOF0Global_H_local0, DOF1Global_H_local1, x_in);

            // - create a matrix of the lever in the global frame
            Real L = l_adaptativebeamInterpolation->getLength(i);
            // - rotate the levers in the global frame
            Vec3 lev(-L/3.0,0.0,0.0);
            Vec3 Lev00_global = -DOF0Global_H_local0.getOrigin();
            Vec3 Lev01_global = DOF0Global_H_local0.getOrientation().rotate(lev) - DOF0Global_H_local0.getOrigin();
            lev[0]=L/3;
            Vec3 Lev12_global = DOF1Global_H_local1.getOrientation().rotate(lev) - DOF1Global_H_local1.getOrigin();
            Vec3 Lev13_global = -DOF1Global_H_local1.getOrigin();

            // create matrices:
            Mat3 Lev00_mat, Lev01_mat, Lev12_mat,Lev13_mat ;
            createCrossMatrix(Lev01_global, Lev01_mat);
            createCrossMatrix(Lev12_global, Lev12_mat);
            createCrossMatrix(Lev00_global, Lev00_mat);
            createCrossMatrix(Lev13_global, Lev13_mat);



            //////////// => NEED TO CHANGE FOR MATRICES  ///////////////////
            //  Matrix SplineP_J_DOFs:
            //_        _   _  _         _                 _   _    _
            //|  dP[0] |   |  |I   Lev00 |                 |  | DX0 |
            //|        |   |  |          |                 |  |     |
            //|  dP[1] |   |  |I   Lev01 |                 |  | DT0 |
            //|        | = |  -          -   _         _   |  |     |
            //|  dP[2] |   |                 |I   Lev12 |  |  | DX1 |
            //|        |   |                 |          |  |  |     |
            //|  dP[3] |   |                 |I   Lev13 |  |  | DT1 |
            // -      -     -                             -    -   -

            // dF[i] += K[i][j] dP[j] (* kfactor)

            // F_DOF  = (SplinePoint_J_DOFs)^T * [dF[0]    dF[1]  dF[2]  dF[2] ]^T

            Mat<4,4,Mat3> SplineP_J_DOFs;
            SplineP_J_DOFs[0][0].identity(); SplineP_J_DOFs[0][1] = Lev00_mat; SplineP_J_DOFs[0][2].clear(); SplineP_J_DOFs[0][3].clear();
            SplineP_J_DOFs[1][0].identity(); SplineP_J_DOFs[1][1] = Lev01_mat; SplineP_J_DOFs[1][2].clear(); SplineP_J_DOFs[1][3].clear();
            SplineP_J_DOFs[2][0].clear(); SplineP_J_DOFs[2][1].clear(); SplineP_J_DOFs[2][2].identity(); SplineP_J_DOFs[2][3] = Lev12_mat;
            SplineP_J_DOFs[3][0].clear(); SplineP_J_DOFs[3][1].clear(); SplineP_J_DOFs[3][2].identity(); SplineP_J_DOFs[3][3] = Lev13_mat;


            Mat<4,4,Mat3> Result;
            for (unsigned int l=0; l<4;l++)// block lines
            {
                for (unsigned int c=0; c<4;c++) // block columns
                {
                    Result[l][c].clear();
                    for (unsigned int j=0; j<4;j++)
                    {
                        for (unsigned int k=0; k<4;k++)
                        {
                            Result[l][c] += SplineP_J_DOFs[j][l].transposed() * K[j][k] * SplineP_J_DOFs[k][c];
                        }
                    }

                }
            }





            ////////////////////////////
           // 5. compute the equivalent stiffness on the rigid rotation of control point (apply DJt on rigid at fixed forces on spline)
           // -compute the equivalent forces on the spline control point (apply Jt on spline map)

           Vec3 F0, F1, F2, F3;
           computeJtSpline(childF, P0,P1,P2,P3, F0, F1, F2, F3);
           Mat3 F0_mat, F1_mat, F2_mat, F3_mat;
           createCrossMatrix(F0, F0_mat);
           createCrossMatrix(F1, F1_mat);
           createCrossMatrix(F2, F2_mat);
           createCrossMatrix(F3, F3_mat);




           // force the symmetry
            /*
          Mat3 test1, test3;
           test1= (F0_mat*Lev00_mat+F1_mat*Lev01_mat);
           std::cout<<"^^^^$$$$$$$$$$$$$$$$$^^^^\n test1 ="<<test1<< "  symetrized  = " <<(test1+test1.transposed())*0.5<<std::endl;
           Result[1][1] -= (test1+test1.transposed())*0.5;
           test3 = (F2_mat*Lev12_mat+F3_mat*Lev13_mat);
           Result[3][3] -= (test3+test3.transposed())*0.5;



       // real derivation
           Result[1][1] -= (F0_mat*Lev00_mat+F1_mat*Lev01_mat);
           Result[3][3] -= (F2_mat*Lev12_mat+F3_mat*Lev13_mat);

*/
           ////////////////////////////
          // 6. put in matrix K_geom


           for (unsigned j=0; j<2; j++)
           {

               for (unsigned int k=0; k<2;k++)
               {
                   for(unsigned int l=0; l<3;l++)
                   {
                       for(unsigned int c=0; c<3;c++)
                       {

                           // translation translation
                           K_geom.add(Nin*IdxNode[j]+l  , Nin*IdxNode[k]+c  , Result[2*j  ][2*k  ][l][c]);
                           // translation rotation
                           K_geom.add(Nin*IdxNode[j]+l  , Nin*IdxNode[k]+c+3, Result[2*j  ][2*k+1][l][c]);
                           // rotation translation
                           K_geom.add(Nin*IdxNode[j]+l+3, Nin*IdxNode[k]+c  , Result[2*j+1][2*k  ][l][c]);
                           // rotation rotation
                           K_geom.add(Nin*IdxNode[j]+l+3, Nin*IdxNode[k]+c+3, Result[2*j+1][2*k+1][l][c]);
                       }
                   }

               }

           }






       } // if WE Compute the stiffness on this beam
   }//iterate on beams

   K_geom.compress();
}


template <class TIn, class TOut>
const sofa::linearalgebra::BaseMatrix* BeamLengthMapping<TIn, TOut>::getK()
{
    return &K_geom;
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
void BeamLengthMapping<TIn, TOut>::computeDJtSpline(const Real &f_input, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3& P3, Mat<4,4,Mat3> &mat)
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

    static Real a[][4] ={ {(-3*(1-x1)*(1-x1)), (3-12*x1+9*x1*x1), (6*x1-9*x1*x1), (3*x1*x1)},
                             {(-3*(1-x2)*(1-x2)), (3-12*x2+9*x2*x2), (6*x2-9*x2*x2), (3*x2*x2)},
                             {(-3*(1-x3)*(1-x3)), (3-12*x3+9*x3*x3), (6*x3-9*x3*x3), (3*x3*x3)},
                             {(-3*(1-x4)*(1-x4)), (3-12*x4+9*x4*x4), (6*x4-9*x4*x4), (3*x4*x4)} };



    static Real B= sqrt(30.0);


    Vec3 IP[4];

    IP[0] = P0*a[0][0] + P1*a[0][1] + P2*a[0][2] + P3*a[0][3]; //IP1
    IP[1] = P0*a[1][0] + P1*a[1][1] + P2*a[1][2] + P3*a[1][3]; //IP2
    IP[2] = P0*a[2][0] + P1*a[2][1] + P2*a[2][2] + P3*a[2][3]; //IP3
    IP[3] = P0*a[3][0] + P1*a[3][1] + P2*a[3][2] + P3*a[3][3]; //IP4


    // 1 compute the derivatives of these relations:
    //F_IP[0] = IP[0]*f_input*(18.0 + B)/(72.0*IP[0].norm());

    Mat3 K_IP[4];

    for(unsigned int k=0; k<4; k++){
        K_IP[k].clear();
        Vec3 n = IP[k]/IP[k].norm();
        for (unsigned int i =0; i<3; i++){
            for (unsigned int j=0; j<3; j++){

                if( i==j )
                    K_IP[k][i][j] = 1.0 - n[i]*n[j];
                else
                    K_IP[k][i][j] =    - n[i]*n[j];
            }
        }

        if(k<2)
            K_IP[k] *= f_input*(18.0 + B)/(72.0*IP[k].norm());
        else
            K_IP[k] *= f_input*(18.0 - B)/(72.0*IP[k].norm());

    }

    for(unsigned int k=0; k<4; k++){

        // 2 apply the same linear relationship:
        //    F0 = F_IP0 * a00  + F_IP1 * a10  + F_IP2 * a20  + F_IP3 * a30 ;

        for(unsigned int l=0; l<4;l++){
            mat[k][l].clear();

            for(unsigned int m=0; m<4;m++){
                mat[k][l] +=  K_IP[m]* a[m][k] * a[m][l];  // A^T * K * A
            }
        }

    }
}




template <class TIn, class TOut>
void BeamLengthMapping< TIn, TOut>::draw(const VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowMappings())
        return;
}


} /// namespace _beamlengthmapping_

} /// namespace mapping

} /// namespace component

} /// namespace sofa

#endif  /* SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_INL */
