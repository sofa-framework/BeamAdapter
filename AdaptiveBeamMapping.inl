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
#ifndef SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_INL
#define SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_INL

#include "AdaptiveBeamMapping.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <string>
#include <sofa/core/Mapping.inl>

#include <sofa/simulation/common/AnimateBeginEvent.h>


namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;





template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::apply(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
	VecCoord& out = *dOut.beginEdit();
	const InVecCoord& in= dIn.getValue();

	x_buf_used=false;

	// => dans le cas où on utilise un controller adaptatif il faut redistribuer les points à chaque pas de temps...
	if (useCurvAbs.getValue() && !contactDuplicate.getValue())
		computeDistribution();


	if (!isSubMapping)
	{
		this->toModel->resize( pointBeamDistribution.size() );
		out.resize(pointBeamDistribution.size());
	}

	for (unsigned int i=0; i<pointBeamDistribution.size(); i++)
	{
		PosPointDefinition  ppd = pointBeamDistribution[i];
		sofa::defaulttype::Vec<3, InReal> pos;
		const Vec3 localPos(0.,ppd.baryPoint[1],ppd.baryPoint[2]);
		m_adaptativebeamInterpolation->interpolatePointUsingSpline(ppd.beamId, ppd.baryPoint[0], localPos, in, pos );
		if(isSubMapping){
			if(idPointSubMap.size()>0)
				out[idPointSubMap[i]] = pos;}
		else
			out[i] = pos;

	}


	dOut.endEdit();

}

template <class TIn, class TOut>
//void AdaptiveBeamMapping< TIn, TOut>::applyJ( typename Out::VecDeriv& out, const typename In::VecDeriv& in )
void AdaptiveBeamMapping< TIn, TOut>::applyJ(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{

	VecDeriv& out = *dOut.beginEdit();
	const InVecDeriv& in= dIn.getValue();

	//const typename In::VecCoord& x = *this->fromModel->getX();
	Data<InVecCoord>& dataInX = *this->getFromModel()->write(sofa::core::VecCoordId::position());
	InVecCoord& x = *dataInX.beginEdit();

	InVecCoord x_buf2;

	if(x_buf_used)
	{
		// TODO : solve this problem during constraint motion propagation !!
		x_buf2 = x;
		x = x_buf;
	}

	if (out.size() != pointBeamDistribution.size() && !isSubMapping)
		out.resize(pointBeamDistribution.size());

	for (unsigned int i=0; i<pointBeamDistribution.size(); i++)
	{
		PosPointDefinition  ppd = pointBeamDistribution[i];
		//1. get the indices
		unsigned int IdxNode0, IdxNode1;
		m_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);


		SpatialVector v_DOF0, v_DOF1;
		//        v_DOF0.setLinearVelocity ( in[IdxNode0].getVCenter() );
		//        v_DOF0.setAngularVelocity(in[IdxNode0].getVOrientation());
		//        v_DOF1.setLinearVelocity ( in[IdxNode1].getVCenter() );
		//        v_DOF1.setAngularVelocity(in[IdxNode1].getVOrientation());
		v_DOF0.setLinearVelocity (In::getDPos(in[IdxNode0]) );
		v_DOF0.setAngularVelocity(In::getDRot(in[IdxNode0]));
		v_DOF1.setLinearVelocity (In::getDPos(in[IdxNode1]) );
		v_DOF1.setAngularVelocity(In::getDRot(in[IdxNode1]));

		Deriv Vresult;

		applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);

		if(isSubMapping){
			if(idPointSubMap.size()>0)
				out[idPointSubMap[i]] = Vresult;}
		else
			out[i] = Vresult;
	}
	if(x_buf_used)
	{
		x = x_buf2;
		x_buf_used = false;
	}


	dOut.endEdit();
	dataInX.endEdit();
}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJT(const core::MechanicalParams* /*mparams*/ /* PARAMS FIRST */, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{

	InVecDeriv& out = *dOut.beginEdit();
	const VecDeriv& in= dIn.getValue();

	//const typename In::VecCoord& x = *this->fromModel->getX();
	const Data<InVecCoord>& dataInX = *this->getFromModel()->read(sofa::core::ConstVecCoordId::position());
	const InVecCoord& x = dataInX.getValue();


	for (unsigned int i=0; i<pointBeamDistribution.size(); i++)
	{
		PosPointDefinition  ppd = pointBeamDistribution[i];
		//1. get the indices
		unsigned int IdxNode0, IdxNode1;
		m_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);


		Deriv finput;
		if(isSubMapping){
			if(idPointSubMap.size()>0)
				finput = in[idPointSubMap[i]];}
		else
			finput = in[i];

		SpatialVector FNode0, FNode1;
		applyJTonPoint(i, finput, FNode0, FNode1, x);

		//2. put the result in out vector computes the equivalent forces on nodes + rotate to Global Frame from DOF frame

		In::setDPos(out[IdxNode0], In::getDPos(out[IdxNode0]) + FNode0.getForce());// out[IdxNode0].getVCenter()+=FNode0.getForce();
		In::setDPos(out[IdxNode1], In::getDPos(out[IdxNode1]) + FNode1.getForce());// out[IdxNode1].getVCenter()+=FNode1.getForce();
		In::setDRot(out[IdxNode0], In::getDRot(out[IdxNode0]) + FNode0.getTorque());// out[IdxNode0].getVOrientation()+=FNode0.getTorque();
		In::setDRot(out[IdxNode1], In::getDRot(out[IdxNode1]) + FNode1.getTorque());// out[IdxNode1].getVOrientation()+=FNode1.getTorque();

	}

	dOut.endEdit();

}


// AdaptiveBeamMapping::applyJT(InMatrixDeriv& out, const OutMatrixDeriv& in) //
// this function propagate the constraint through the Adaptive Beam mapping :
// if one constraint along (vector n) with a value (v) is applied on the childModel (like collision model)
// then this constraint is transformed by (Jt.n) with value (v) for the rigid model
// note : the value v is not propagated through the mapping
template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJT(const core::ConstraintParams* /*cparams*/ /* PARAMS FIRST */, Data<InMatrixDeriv>& dOut, const Data<OutMatrixDeriv>& dIn)
{

	InMatrixDeriv& out = *dOut.beginEdit();
	const OutMatrixDeriv& in = dIn.getValue();


	//const typename In::VecCoord& x = *this->fromModel->getX();
	const Data<InVecCoord>& dataInX = *this->getFromModel()->read(sofa::core::ConstVecCoordId::position());
	const InVecCoord& x = dataInX.getValue();

	x_buf_used=false;
	x_buf = x ;


	//////////// What's for ?? it seems not useful//////////
	bool proximity_lever = false;
	if (proximity.getValue() > 0.0)
		proximity_lever = true;

	if (proximity_lever && this->f_printLog.getValue() )
	{
		serr<<" the constraints are contact at the surface of the beam (not at their center)"<<sendl;
	}
	/////////////////////////////////////


	typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();

	for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
	{

		typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();

		for (typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin(); colIt != colItEnd; ++colIt)
		{
			typename In::MatrixDeriv::RowIterator o = out.writeLine(rowIt.index());

			unsigned int indexIn = colIt.index();
			const Deriv data = colIt.val();

			if (isSubMapping){
				// look if we get the indexIn in the idPointSubMap:
				unsigned int i=0;
				while( i<idPointSubMap.size() && idPointSubMap[i]!=indexIn)
				{
					i++;
				}

				if (i<idPointSubMap.size())
				{
					indexIn = i;
				}
				else
				{
					continue;
				}
			}


			if(indexIn<pointBeamDistribution.size()){
				PosPointDefinition  ppd = pointBeamDistribution[indexIn];
				unsigned int IdxNode0, IdxNode1;
				m_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

				SpatialVector FNode0, FNode1;

				applyJTonPoint(indexIn, data, FNode0, FNode1, x);

				// Compute the mapped Constraint on the beam nodes ///
				InDeriv direction0;
				In::setDPos(direction0,FNode0.getForce());// direction0.getVCenter() = FNode0.getForce();
				In::setDRot(direction0,FNode0.getTorque());// direction0.getVOrientation() = FNode0.getTorque();
				InDeriv direction1;
				In::setDPos(direction1,FNode1.getForce());// direction1.getVCenter() = FNode1.getForce();
				In::setDRot(direction1,FNode1.getTorque());// direction1.getVOrientation() = FNode1.getTorque();

				o.addCol(IdxNode0, direction0);
				o.addCol(IdxNode1, direction1);

			}
			else
			{
				if ( this->f_printLog.getValue() )
					serr<<" WARNING = wrong index in VecConst in"<<sendl;
				break;
			}


		}
	}


	dOut.endEdit();
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::init()
{
	if (!m_adaptativebeamInterpolation) 
	   	m_adaptativebeamInterpolation.set(dynamic_cast<core::objectmodel::BaseContext *>(this->getContext())->get<BInterpolation>());

	if (!m_adaptativebeamInterpolation)
		serr<<"No Beam Interpolation found, the component can not work!"<<sendl;
}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::bwdInit()
{

	const sofa::helper::vector<Vec3>& pts = points.getValue();

	if (pts.size() == 0)
	{
		const typename Out::VecCoord& xTo = *this->toModel->getX();

		if( xTo.size()==0)
		{
			serr<<" Warning no point defined in the AdaptiveBeamMapping "<<sendl;
			// _problem = true;
		}
		else
		{
			sofa::helper::vector<Vec3>& pts2 = *(points.beginEdit());

			sout<<"no point defined in the AdaptiveBeamMapping - uses positions defined by Mechanical State"<<sendl;
			for(unsigned int i=0; i<xTo.size();i++)
			{
				Vec3 p(xTo[i][0], xTo[i][1], xTo[i][2]);
				pts2.push_back(p);
			}
			points.endEdit();
		}
	}


	bool curvAbs = this->useCurvAbs.getValue();
	if (curvAbs)
	{
		for (unsigned int i=0; i<pts.size()-1; i++)
		{

			// TODO : really necessary ?
			if( pts[i][0]>pts[i+1][0])
				sout<<" when using useCurvAbs==true, points must be sorted according to their curvAbs "<<sendl;
		}
	}

	m_adaptativebeamInterpolation->bwdInit();
	computeDistribution();
	if (!isSubMapping)
	{
		core::Mapping< TIn, TOut>::init();
	}
	if(contactDuplicate.getValue()==true)
	{
		const std::string path = nameOfInputMap.getValue();
		this->fromModel->getContext()->get(m_inputMapping, sofa::core::objectmodel::BaseContext::SearchRoot);
		if(m_inputMapping==NULL)
			serr<<"WARNING : can not found the input  Mapping"<<sendl;
		else
			sout<<"input Mapping named "<<m_inputMapping->getName()<<" is found"<<sendl;
	}

}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::beginAddContactPoint()
{
	//std::cerr<<"beginAddContactPoint called on mapping "<<std::endl;
	points.beginEdit()->clear();
	points.endEdit();
	pointBeamDistribution.clear();
}


template <class TIn, class TOut>
int AdaptiveBeamMapping< TIn, TOut>::addContactPoint(const Vec3& bary)
{
	unsigned int index = points.getValue().size();
	points.beginEdit()->push_back(bary);
	points.endEdit();
	if(this->toModel->getX()->size() <= index)
		this->toModel->resize(index+1);
	return index;
}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::computeIdxAndBaryCoordsForAbs(unsigned int &b, Real &x_bary, const Real &x_abs )
{

	InReal x_abs_input = (InReal) x_abs;
	InReal x_bary_output = (InReal) x_bary;
	// new :
	m_adaptativebeamInterpolation->getBeamAtCurvAbs(x_abs_input,b,x_bary_output);
	x_bary = (Real) x_bary_output;

}



template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::computeDistribution()
{
	//The normal procedure is givng "points", computeDistribution comppute theses points and put
	// on the baryPoints which are "pointBeamDistribution"
	//If the mapping is a barycentric one, that mean a list of baryPoints are already setted
	//Then no need to recompute this
	if(!isBarycentricMapping)
	{
		bool curvAbs = useCurvAbs.getValue();
		const sofa::helper::vector<Vec3>& pts = points.getValue();
		pointBeamDistribution.clear();

		unsigned int numBeams = m_adaptativebeamInterpolation->getNumBeams();
		if (numBeams==0)
		{
			if (this->f_printLog.getValue())
				serr<<" no beams found in adaptBeamInterpolation in BeamInterpolation  named"<<m_adaptativebeamInterpolation->getName()<<sendl;
			return;
		}

		if (curvAbs)
		{
			int ptsPerBeam = nbPointsPerBeam.getValue();
			if(ptsPerBeam)
			{	// Recreating the distribution based on the current sampling of the beams
				unsigned int nbBeams = m_adaptativebeamInterpolation->getNumBeams();
				for(unsigned int b=0; b<nbBeams; ++b)
				{
					for(int i=0; i<ptsPerBeam; ++i)
					{
						PosPointDefinition beamDistrib;
						beamDistrib.beamId = b;
						beamDistrib.baryPoint[0] = (Real)i / ptsPerBeam;
						beamDistrib.baryPoint[1] = 0.0;
						beamDistrib.baryPoint[2] = 0.0;

						pointBeamDistribution.push_back(beamDistrib);
					}
				}

				if(nbBeams)
				{	// Last point
					PosPointDefinition beamDistrib;
					beamDistrib.beamId = nbBeams-1;
					beamDistrib.baryPoint[0] = 1.0;
					beamDistrib.baryPoint[1] = 0.0;
					beamDistrib.baryPoint[2] = 0.0;

					pointBeamDistribution.push_back(beamDistrib);
				}

				BaseMeshTopology* topo = getContext()->getMeshTopology();
				if(topo)
				{
					topo->clear();
					int nbEdges = pointBeamDistribution.size() - 1;
					for(int i=0; i<nbEdges; ++i)
						topo->addEdge(i, i+1);
				}
			}
			else
			{	// We use the points Data
				for (unsigned int i=0; i<pts.size(); i++)
				{
					unsigned int b=0;
					Real x_abs = pts[i][0];
					Real x_bary = 0.0;

					computeIdxAndBaryCoordsForAbs(b, x_bary,  x_abs );

					PosPointDefinition beamDistrib;
					beamDistrib.beamId = b;
					beamDistrib.baryPoint[0] = x_bary;
					beamDistrib.baryPoint[1] = pts[i][1];
					beamDistrib.baryPoint[2] = pts[i][2];

					pointBeamDistribution.push_back(beamDistrib);
				}
			}
		}
		else
		{
			for (unsigned int i=0; i<pts.size(); i++)
			{
				PosPointDefinition beamDistrib;
				beamDistrib.beamId = (int) floor(pts[i][0]);
				if ( (beamDistrib.beamId>numBeams-1 || beamDistrib.beamId<0.0 ) && this->f_printLog.getValue()  )
					serr<<"WARNING : points["<<i<<"][0] = "<<beamDistrib.baryPoint[0]<<" is defined outside of the beam length"<<sendl;
				beamDistrib.baryPoint[0] = pts[i][0] - floor(pts[i][0]);
				beamDistrib.baryPoint[1] = pts[i][1];
				beamDistrib.baryPoint[2] = pts[i][2];
				pointBeamDistribution.push_back(beamDistrib);
			}
		}
	}
}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::computeJacobianOnPoint(unsigned int i, const typename In::VecCoord& x)
{
	//	/////// TEST : calcul d'une jacobienne:
	Mat3x12 J;
	Mat12x3 Jt;

	for (unsigned int j=0; j<3; j++)
	{
		Deriv Id(0,0,0), Vresult;
		Id[j]=1.0;
		SpatialVector v_DOF0, v_DOF1;

		//3 colonnes
		v_DOF0.clear();
		v_DOF0.setLinearVelocity(Id);
		v_DOF1.clear();
		applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
		J(0,j)=Vresult[0]; J(1,j)=Vresult[1]; J(2,j)=Vresult[2];
		//3 colonnes
		v_DOF0.clear();
		v_DOF0.setAngularVelocity(Id);
		v_DOF1.clear();
		applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
		J(0,j+3)=Vresult[0]; J(1,j+3)=Vresult[1]; J(2,j+3)=Vresult[2];
		//3 colonnes
		v_DOF0.clear();
		v_DOF1.clear();
		v_DOF1.setLinearVelocity(Id);

		applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
		J(0,j+6)=Vresult[0]; J(1,j+6)=Vresult[1]; J(2,j+6)=Vresult[2];
		//3 colonnes
		v_DOF0.clear();
		v_DOF1.clear();
		v_DOF1.setAngularVelocity(Id);
		applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
		J(0,j+9)=Vresult[0]; J(1,j+9)=Vresult[1]; J(2,j+9)=Vresult[2];


		SpatialVector F_DOF0, F_DOF1;
		applyJTonPoint(i, Id, F_DOF0, F_DOF1, x);
		Jt(0,j)=F_DOF0.getForce()[0]; Jt(1,j)=F_DOF0.getForce()[1];  Jt(2,j) =F_DOF0.getForce()[2];
		Jt(3,j)=F_DOF0.getTorque()[0];Jt(4,j)=F_DOF0.getTorque()[1]; Jt(5,j) =F_DOF0.getTorque()[2];
		Jt(6,j)=F_DOF1.getForce()[0]; Jt(7,j)=F_DOF1.getForce()[1];  Jt(8,j) =F_DOF1.getForce()[2];
		Jt(9,j)=F_DOF1.getTorque()[0];Jt(10,j)=F_DOF1.getTorque()[1];Jt(11,j)=F_DOF1.getTorque()[2];

	}
}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJonPoint(unsigned int i, SpatialVector& VNode0input, SpatialVector& VNode1input, Deriv& vOutput, const InVecCoord& x)
{

	// 1. get the curvilinear abs;
	PosPointDefinition  ppd = pointBeamDistribution[i];

	// 2. get the indices
	unsigned int IdxNode0, IdxNode1;
	m_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

	// 3. get the transform to DOF in global frame from local frame
	Transform DOF0Global_H_local0, DOF1Global_H_local1;
	m_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(ppd.beamId, DOF0Global_H_local0, DOF1Global_H_local1, x);

	// 4. project the velocities in local frame:
	SpatialVector v_local0, v_local1;
	v_local0 = DOF0Global_H_local0.inversed()*VNode0input;
	v_local1 = DOF1Global_H_local1.inversed()*VNode1input;

	// 5. Computes the local velocities of the 4 points of the spline
	Real L = m_adaptativebeamInterpolation->getLength(ppd.beamId);
	Vec3 lever(L/3,0,0);
	Vec3 V0, V1, V2, V3;
	V0 = v_local0.getLinearVelocity();
	V1 = V0 - lever.cross(v_local0.getAngularVelocity());
	V3 = v_local1.getLinearVelocity();
	V2 = V3 + lever.cross(v_local1.getAngularVelocity());

	// 6. Rotate back the vectors in the global frame
	V0 = DOF0Global_H_local0.getOrientation().rotate(V0);
	V1 = DOF0Global_H_local0.getOrientation().rotate(V1);
	V2 = DOF1Global_H_local1.getOrientation().rotate(V2);
	V3 = DOF1Global_H_local1.getOrientation().rotate(V3);

	const Vec3 localPos(0.,ppd.baryPoint[1],ppd.baryPoint[2]);
	if(localPos.norm() > L*1e-10)
	{
		lever = localPos;

		// 7. compute the effect of the angular velocities on the points that are not aligned along the center line
		Vec3 DV0, DV3;
		DV0 = DOF0Global_H_local0.getOrientation().rotate( - lever.cross( v_local0.getAngularVelocity() ) );
		DV3 = DOF1Global_H_local1.getOrientation().rotate( - lever.cross( v_local1.getAngularVelocity() ) );

		// uses spline to interpolate:
		Real bx = ppd.baryPoint[0];
		Real a0=(1-bx)*(1-bx)*(1-bx);
		Real a1=3*bx*(1-bx)*(1-bx);
		Real a2=3*bx*bx*(1-bx);
		Real a3=bx*bx*bx;
		vOutput = V0*a0 + V1*a1 + V2*a2 + V3*a3;
		vOutput += DV0*(a0+a1) + DV3*(a2+a3);

	}
	else
	{
		// uses spline to interpolate:
		Real bx = ppd.baryPoint[0];
		vOutput = V0*(1-bx)*(1-bx)*(1-bx) + V1*3*bx*(1-bx)*(1-bx) + V2*3*bx*bx*(1-bx) + V3*bx*bx*bx;
	}
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJTonPoint(unsigned int i, const Deriv& finput, SpatialVector& FNode0output, SpatialVector& FNode1output, const InVecCoord& x)
{

    PosPointDefinition  ppd = pointBeamDistribution[i];
    const Vec3 localPos(0.,ppd.baryPoint[1],ppd.baryPoint[2]);
    const Vec3 Fin(finput[0], finput[1], finput[2]);

    m_adaptativebeamInterpolation->MapForceOnNodeUsingSpline(ppd.beamId, ppd.baryPoint[0], localPos, x, Fin, FNode0output, FNode1output );
/*
	//1. get the curvilinear abs;
	PosPointDefinition  ppd = pointBeamDistribution[i];
	Real bx = ppd.baryPoint[0];
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
	m_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(ppd.beamId, DOF0Global_H_local0, DOF1Global_H_local1, x);

	//rotate back to local frame
	SpatialVector f0, f1,f2,f3;
	f0.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F0) );
	f1.setForce( DOF0Global_H_local0.getOrientation().inverseRotate(F1) );
	f2.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F2) );
	f3.setForce( DOF1Global_H_local1.getOrientation().inverseRotate(F3) );

	// computes the torque created on DOF0 and DOF1 by f1 and f2
	Real L = m_adaptativebeamInterpolation->getLength(ppd.beamId);

	const Vec3 localPos(0.,ppd.baryPoint[1],ppd.baryPoint[2]);
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
      */

}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::draw(const core::visual::VisualParams* vparams)
{
	if (!vparams->displayFlags().getShowMappings()) return;
}



template<>
void AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes >::apply(const core::MechanicalParams* /* PARAMS FIRST */, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn );

} // namespace mapping

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_INL */
