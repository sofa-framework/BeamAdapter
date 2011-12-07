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

#include "AdaptiveBeamMapping.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/Mapping.inl>

namespace sofa
{

namespace component
{

namespace mapping
{

SOFA_DECL_CLASS(AdaptiveBeamMapping)

using namespace defaulttype;
using namespace core;
using namespace core::behavior;





// Register in the Factory
int AdaptiveBeamMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs")
#ifndef SOFA_FLOAT
.add< AdaptiveBeamMapping<Rigid3dTypes, Vec3dTypes   > >()
.add< AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes > >()
#endif
#ifndef SOFA_DOUBLE
.add< AdaptiveBeamMapping< Rigid3fTypes, Vec3fTypes > >()
#endif
#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
.add< AdaptiveBeamMapping< Rigid3dTypes, Vec3fTypes > >()
.add< AdaptiveBeamMapping< Rigid3fTypes, Vec3dTypes > >()
#endif
#endif
//.add< AdaptiveBeamMapping< Rigid3fTypes, Vec3dTypes > >()
//.add< AdaptiveBeamMapping< Rigid3fTypes, Vec3fTypes > >()
;





#ifndef SOFA_FLOAT
template<>
void AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes >::apply(const core::MechanicalParams* /* PARAMS FIRST */, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn )
{
        VecCoord& out = *dOut.beginEdit();
    const InVecCoord& in= dIn.getValue();

    x_buf_used=false;

// => dans le cas où on utilise un controller adaptatif il faut redistribuer les points à chaque pas de temps...
    if (useCurvAbs.getValue() && !contactDuplicate.getValue())
        computeDistribution();



    out.resize(pointBeamDistribution.size());

    for (unsigned int i=0; i<pointBeamDistribution.size(); i++)
    {
        PosPointDefinition  ppd = pointBeamDistribution[i];
        Transform posTransform;
    	Vec3 localPos(0.,ppd.baryPoint[1],ppd.baryPoint[2]);
        m_adaptativebeamInterpolation->InterpolateTransformUsingSpline(ppd.beamId,ppd.baryPoint[0],localPos, in, posTransform );
        out[i].getCenter() = posTransform.getOrigin();
        out[i].getOrientation() = posTransform.getOrientation();

    }

    dOut.endEdit();
}

template<>
void AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes >::applyJonPoint(unsigned int i, SpatialVector& VNode0input, SpatialVector& VNode1input, Deriv& vOutput, const  InVecCoord& x)
{
    //1. get the curvilinear abs;
    PosPointDefinition  ppd = pointBeamDistribution[i];

    //2. get the indices
    unsigned int IdxNode0, IdxNode1;
    m_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

    //3. get the transform to DOF in global frame from local frame
    Transform DOF0Global_H_local0, DOF1Global_H_local1;
    m_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(ppd.beamId, DOF0Global_H_local0, DOF1Global_H_local1, x);

    //4. project the velocities in local frame:
    SpatialVector v_local0, v_local1;
    v_local0 = DOF0Global_H_local0.inversed()*VNode0input;
    v_local1 = DOF1Global_H_local1.inversed()*VNode1input;

    //5. Computes the local velocities of the 4 points of the spline
    Real L = m_adaptativebeamInterpolation->getLength(ppd.beamId);
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
    //vOutput.getVCenter() = V0*a0 + V1*a1 + V2*a2 + V3*a3;
    //vOutput.getVOrientation() = W0*(a0+a1) + W3*(a2+a3);
    Rigid3dTypes::setDPos(vOutput,V0*a0 + V1*a1 + V2*a2 + V3*a3);
    Rigid3dTypes::setDRot(vOutput,W0*(a0+a1) + W3*(a2+a3));
}

template<>
void AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes >::applyJTonPoint(unsigned int i, const Deriv& finput, SpatialVector& FNode0output, SpatialVector& FNode1output, const  InVecCoord& x)
{
	//1. get the curvilinear abs;
	PosPointDefinition  ppd = pointBeamDistribution[i];
	Real bx = ppd.baryPoint[0];
	Real a0=(1-bx)*(1-bx)*(1-bx);
	Real a1=3*bx*(1-bx)*(1-bx);
	Real a2=3*bx*bx*(1-bx);
	Real a3=bx*bx*bx;

	//2. computes a force on the 4 points of the spline:
	Vec3 F0, F1, F2, F3, C0, C3;

//	F0 = finput.getVCenter()*a0;
//	F1 = finput.getVCenter()*a1;
//	F2 = finput.getVCenter()*a2;
//	F3 = finput.getVCenter()*a3;
//	C0 = finput.getVOrientation()*(a0+a1);
//	C3 = finput.getVOrientation()*(a2+a3);

	F0 = Rigid3dTypes::getDPos(finput)*a0;
	F1 = Rigid3dTypes::getDPos(finput)*a1;
	F2 = Rigid3dTypes::getDPos(finput)*a2;
	F3 = Rigid3dTypes::getDPos(finput)*a3;
	C0 = Rigid3dTypes::getDRot(finput)*(a0+a1);
	C3 = Rigid3dTypes::getDRot(finput)*(a2+a3);


	//  std::cout<<" ************ applyJTonPoint *************"<<std::endl;
	//  std::cout<<" finput="<<finput<<"  -  F0 = "<<F0<<"  -  F1 = "<<F1<<"  - F2 = "<<F2<<"  - F3 = "<<F3<<std::endl;


	//3. influence of these forces on the nodes of the beam    => TODO : simplify the computations !!!
	Transform DOF0Global_H_local0, DOF1Global_H_local1;
	m_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(ppd.beamId, DOF0Global_H_local0, DOF1Global_H_local1, x);

	//std::cout<<" applyJTonPoint : DOF0Global_H_local0="<<DOF0Global_H_local0<<"  -  DOF1Global_H_local1="<<DOF1Global_H_local1<<std::endl;


	// rotate back to local frame
	SpatialVector f0, f1,f2,f3;
	f0.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F0) );
	f1.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F1) );
	f2.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F2) );
	f3.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F3) );


	// computes the torque created on DOF0 and DOF1 by f1 and f2
	Real L = m_adaptativebeamInterpolation->getLength(ppd.beamId);
	Vec3 lever(L/3,0,0);
	f0.setTorque( DOF0Global_H_local0.getOrientation().inverseRotate(C0) );
	f1.setTorque(lever.cross(f1.getForce()));
	f2.setTorque(-lever.cross(f2.getForce()));
	f3.setTorque( DOF1Global_H_local1.getOrientation().inverseRotate(C3) );

	//   std::cout<<" f0="<<f0<<"  -  f1 ="<<f1<<" - f2="<<f2<<"  -  f3 ="<<f3<<std::endl;

	//    // back to the DOF0 and DOF1 frame:
	FNode0output = DOF0Global_H_local0 * (f0+f1);
	FNode1output = DOF1Global_H_local1 * (f2+f3);

	//   std::cout<<" FNode0output="<<FNode0output<<"  -  FNode1output ="<<FNode1output<<std::endl;

}

template <>
void AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes >::computeJacobianOnPoint(unsigned int i, const  InVecCoord& x)
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
        v_DOF0.setLinearVelocity(Rigid3dTypes::getDPos(Id));
        v_DOF0.setAngularVelocity(Rigid3dTypes::getDRot(Id));
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
        v_DOF1.setLinearVelocity(Rigid3dTypes::getDPos(Id));
        v_DOF1.setAngularVelocity(Rigid3dTypes::getDRot(Id));
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

    //std::cout<<" ********** TEST J: ********** \n"<<J<<std::endl;
    //std::cout<<" ********** TEST Jt(transposed): ********** \n"<<Jt.transposed()<<std::endl;
    std::cout<<" ********** TEST J-Jt(transposed): ********** \n"<<Test<<std::endl;


}



template <>
int AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes >::addPoint (const Coord& c, int )
{
    int i = points.getValue().size();
    Vec3 test = c.getCenter();

    points.beginEdit()->push_back(test);
    return i;
}


#endif




#ifndef SOFA_FLOAT
template class SOFA_BEAMADAPTER_API AdaptiveBeamMapping<Rigid3dTypes, Vec3dTypes   >;
template class SOFA_BEAMADAPTER_API AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes >;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_BEAMADAPTER_API AdaptiveBeamMapping< Rigid3fTypes, Vec3fTypes >;
#endif

#ifndef SOFA_FLOAT
#ifndef SOFA_DOUBLE
template class SOFA_BEAMADAPTER_API AdaptiveBeamMapping< Rigid3dTypes, Vec3fTypes >;
template class SOFA_BEAMADAPTER_API AdaptiveBeamMapping< Rigid3fTypes, Vec3dTypes >;
#endif
#endif


} // namespace mapping

} // namespace component

} // namespace sofa

