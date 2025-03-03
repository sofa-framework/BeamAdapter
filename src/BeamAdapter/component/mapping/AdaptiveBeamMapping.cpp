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
#define SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEBEAMMAPPING_CPP

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/behavior/MechanicalState.h>
#include <BeamAdapter/config.h>
#include <sofa/core/ObjectFactory.h>

#include <BeamAdapter/component/mapping/AdaptiveBeamMapping.inl>

using namespace sofa::defaulttype;

namespace beamadapter
{

using namespace core;
using namespace core::behavior;

template<>
SOFA_BEAMADAPTER_API void AdaptiveBeamMapping<Rigid3Types, Rigid3Types >::apply(const MechanicalParams*, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn )
{
    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    const InVecCoord& in= dIn.getValue();

    m_isXBufferUsed=false;

    // When using an adaptatif controller, one need to redistribute the points at each time step
    if (d_useCurvAbs.getValue() && !d_contactDuplicate.getValue())
        computeDistribution();

    out.resize(m_pointBeamDistribution.size());

    for (unsigned int i=0; i<m_pointBeamDistribution.size(); i++)
    {
        PosPointDefinition  ppd = m_pointBeamDistribution[i];
        Transform posTransform;
        Vec3 localPos(0.,ppd.baryPoint[1],ppd.baryPoint[2]);
        l_adaptativebeamInterpolation->InterpolateTransformUsingSpline(ppd.beamId,ppd.baryPoint[0],localPos, in, posTransform );
        out[i].getCenter() = posTransform.getOrigin();
        out[i].getOrientation() = posTransform.getOrientation();
    }
}


template<>
SOFA_BEAMADAPTER_API void AdaptiveBeamMapping<Rigid3Types, Rigid3Types >::applyJonPoint(unsigned int i, SpatialVector& VNode0input, SpatialVector& VNode1input, Deriv& vOutput, const  InVecCoord& x)
{
    //1. get the curvilinear abs;
    PosPointDefinition  ppd = m_pointBeamDistribution[i];

    //2. get the indices
    unsigned int IdxNode0, IdxNode1;
    l_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

    //3. get the transform to DOF in global frame from local frame
    Transform DOF0Global_H_local0, DOF1Global_H_local1;
    l_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(ppd.beamId, DOF0Global_H_local0, DOF1Global_H_local1, x);

    //4. project the velocities in local frame:
    SpatialVector v_local0, v_local1;
    v_local0 = DOF0Global_H_local0.inversed()*VNode0input;
    v_local1 = DOF1Global_H_local1.inversed()*VNode1input;

    //5. Computes the local velocities of the 4 points of the spline
    Real L = l_adaptativebeamInterpolation->getLength(ppd.beamId);
    Vec3 lever(L/3,0,0);
    Vec3 V0, V1, V2, V3;
    V0 = v_local0.getLinearVelocity();
    V1 = V0 - lever.cross(v_local0.getAngularVelocity());
    V3 = v_local1.getLinearVelocity();
    V2 = V3 + lever.cross(v_local1.getAngularVelocity());

    //6. Rotate back the vectors in the global frame
    V0 = DOF0Global_H_local0.getOrientation().rotate(V0);
    V1 = DOF0Global_H_local0.getOrientation().rotate(V1);
    V2 = DOF1Global_H_local1.getOrientation().rotate(V2);
    V3 = DOF1Global_H_local1.getOrientation().rotate(V3);

    // 7. Rotate back the angular velocities in the global frame
    Vec3 W0, W3;
    W0 = DOF0Global_H_local0.getOrientation().rotate(v_local0.getAngularVelocity());
    W3 = DOF1Global_H_local1.getOrientation().rotate(v_local1.getAngularVelocity());

    // uses spline to interpolate:
    Real bx = ppd.baryPoint[0];
    Real a0=(1-bx)*(1-bx)*(1-bx);
    Real a1=3*bx*(1-bx)*(1-bx);
    Real a2=3*bx*bx*(1-bx);
    Real a3=bx*bx*bx;
    Rigid3Types::setDPos(vOutput,V0*a0 + V1*a1 + V2*a2 + V3*a3);
    Rigid3Types::setDRot(vOutput,W0*(a0+a1) + W3*(a2+a3));
}


template<>
SOFA_BEAMADAPTER_API void AdaptiveBeamMapping<Rigid3Types, Rigid3Types >::applyJTonPoint(unsigned int i, const Deriv& finput, SpatialVector& FNode0output, SpatialVector& FNode1output, const  InVecCoord& x)
{
    //1. get the curvilinear abs;
    PosPointDefinition  ppd = m_pointBeamDistribution[i];
    Real bx = ppd.baryPoint[0];

    SpatialVector f6DofInput;
    f6DofInput.setForce(Rigid3Types::getDPos(finput));
    f6DofInput.setTorque(Rigid3Types::getDRot(finput));

    l_adaptativebeamInterpolation->MapForceOnNodeUsingSpline(ppd.beamId, bx, Vec3(0,ppd.baryPoint[1],ppd.baryPoint[2]), x,
            f6DofInput, FNode0output, FNode1output);
}


template <>
SOFA_BEAMADAPTER_API void AdaptiveBeamMapping<Rigid3Types, Rigid3Types >::computeJacobianOnPoint(unsigned int i, const  InVecCoord& x)
{
    /////// TEST : calcul d'une jacobienne:
    Mat6x12 J;
    Mat12x6 Jt;

    for (unsigned int j=0; j<6; j++)
    {
        Deriv Id, Vresult;
        Id[j]=1.0;
        SpatialVector v_DOF0, v_DOF1;

        //  6 colonnes
        v_DOF0.clear();
        //v_DOF0.setLinearVelocity(Id.getVCenter());
        //v_DOF0.setAngularVelocity(Id.getVOrientation());
        v_DOF0.setLinearVelocity(Rigid3Types::getDPos(Id));
        v_DOF0.setAngularVelocity(Rigid3Types::getDRot(Id));
        v_DOF1.clear();
        applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
        J(0,j)=Vresult[0]; J(1,j)=Vresult[1]; J(2,j)=Vresult[2]; J(3,j)=Vresult[3]; J(4,j)=Vresult[4]; J(5,j)=Vresult[5];
        //3 colonnes
        //        v_DOF0.clear();
        //        v_DOF0.setLinearVelocity(Id.getVCenter());
        //        v_DOF0.setAngularVelocity(Id.getVOrientation());
        //        v_DOF1.clear();
        //        applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
        //        J(0,j+3)=Vresult[0]; J(1,j+3)=Vresult[1]; J(2,j+3)=Vresult[2]; J(3,j+3)=Vresult[3]; J(4,j+3)=Vresult[4]; J(5,j+3)=Vresult[5];
        //  6 colonnes
        v_DOF0.clear();
        v_DOF1.clear();
        v_DOF1.setLinearVelocity(Rigid3Types::getDPos(Id));
        v_DOF1.setAngularVelocity(Rigid3Types::getDRot(Id));
        applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
        J(0,j+6)=Vresult[0]; J(1,j+6)=Vresult[1]; J(2,j+6)=Vresult[2]; J(3,j+6)=Vresult[3]; J(4,j+6)=Vresult[4]; J(5,j+6)=Vresult[5];
        //    //3 colonnes
        //        v_DOF0.clear();
        //        v_DOF1.clear();
        //        v_DOF1.setLinearVelocity(Id.getVCenter());
        //        v_DOF1.setAngularVelocity(Id.getVOrientation());
        //        applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
        //        J(0,j+9)=Vresult[0]; J(1,j+9)=Vresult[1]; J(2,j+9)=Vresult[2]; J(3,j+9)=Vresult[3]; J(4,j+9)=Vresult[4]; J(5,j+9)=Vresult[5];


        SpatialVector F_DOF0, F_DOF1;
        applyJTonPoint(i, Id, F_DOF0, F_DOF1, x);
        Jt(0,j)=F_DOF0.getForce()[0]; Jt(1,j)=F_DOF0.getForce()[1];  Jt(2,j) =F_DOF0.getForce()[2];
        Jt(3,j)=F_DOF0.getTorque()[0];Jt(4,j)=F_DOF0.getTorque()[1]; Jt(5,j) =F_DOF0.getTorque()[2];
        Jt(6,j)=F_DOF1.getForce()[0]; Jt(7,j)=F_DOF1.getForce()[1];  Jt(8,j) =F_DOF1.getForce()[2];
        Jt(9,j)=F_DOF1.getTorque()[0];Jt(10,j)=F_DOF1.getTorque()[1];Jt(11,j)=F_DOF1.getTorque()[2];

    }
    Mat6x12 Test=J-Jt.transposed();

    dmsg_info()<<" ********** TEST J-Jt(transposed): ********** \n"<<Test;
}

template class SOFA_BEAMADAPTER_API AdaptiveBeamMapping<Rigid3Types, Vec3Types>;
template class SOFA_BEAMADAPTER_API AdaptiveBeamMapping<Rigid3Types, Rigid3Types>;

void registerAdaptiveBeamMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs.")
                             .add< AdaptiveBeamMapping<Rigid3Types, Vec3Types   > >(true) //default template
                             .add< AdaptiveBeamMapping<Rigid3Types, Rigid3Types > >());
}

} // namespace beamadapter
