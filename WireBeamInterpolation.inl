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

#include "WireBeamInterpolation.h"
#include "BeamInterpolation.inl"

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
//#include <sofa/component/topology/GridTopology.h>
//#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/decompose.h>
//#include <sofa/helper/gl/template.h>
//#include <sofa/helper/gl/Axis.h>
//#include <sofa/helper/rmath.h>
//#include <assert.h>
//#include <iostream>
//#include <set>
//#include <sofa/helper/system/gl.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

//#include <sofa/defaulttype/SolidTypes.inl>
#include <sofa/helper/OptionsGroup.h>

#include <sofa/helper/gl/Cylinder.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/gl/Axis.h>
//#include <sofa/simulation/common/Node.h>

#include <algorithm>

#define _max_(a,b) (((a)>(b))?(a):(b))
#define _min_(a,b) (((a)<(b))?(a):(b))

namespace sofa
{

namespace component
{

namespace fem
{


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

	helper::vector<Real> xP_noticeable;
	helper::vector< int> nbP_density;

	m_restShape.get()->getSamplingParameters(xP_noticeable, nbP_density);
}


template <class DataTypes>
 void WireBeamInterpolation<DataTypes>::bwdInit()
{
    Inherited::bwdInit();

    if (!this->isControlled())
        serr << "not straightRestShape or this->Edge_List is assigned" << sendl;
    else
        sout << "external controller for this ForceField is detected" << sendl;
}


template<class DataTypes>
void WireBeamInterpolation<DataTypes>::addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1,
                                               const Transform &DOF0_H_Node0, const Transform &DOF1_H_Node1)
{	
	VecElementID &edgeList = *this->m_edgeList.beginEdit();
	vector< double > &lengthList = *this->m_lengthList.beginEdit();
	vector< Transform > &DOF0TransformNode0 = *this->m_DOF0TransformNode0.beginEdit();
	vector< Transform > &DOF1TransformNode1 = *this->m_DOF1TransformNode1.beginEdit();
    vector< Vec2 > &curvAbsList = *this->m_curvAbsList.beginEdit();

    edgeList.push_back(eID);
    lengthList.push_back(length);

    curvAbsList.push_back(Vec2(x0, x1));

    // as an angle is set between DOFs and Beam, they are no more aligned
    this->dofsAndBeamsAligned.setValue(false);
    DOF0TransformNode0.push_back(DOF0_H_Node0);
    DOF1TransformNode1.push_back(DOF1_H_Node1);

	this->m_edgeList.endEdit();
	this->m_lengthList.endEdit();
	this->m_DOF0TransformNode0.endEdit();
	this->m_DOF1TransformNode1.endEdit();
	this->m_curvAbsList.endEdit();
}



template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getRestTransform(unsigned int edgeInList, Transform &local0_H_local1_rest)
{
	this->f_printLog.setValue(true);
	serr<<"WARNING : getRestTransform not implemented for not straightRestShape"<<sendl;

	// the beam is straight: the transformation between local0 and local1 is provided by the length of the beam
	local0_H_local1_rest.set(Vec3(this->m_lengthList.getValue()[edgeInList], 0, 0), Quat());
}


template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest)
{
	if (this->isControlled() && this->m_restShape!=NULL)
	{
        const Vec2 &curvAbs = this->m_curvAbsList.getValue()[edgeInList];

		Real x_middle = (curvAbs.x() + curvAbs.y()) / 2;
		Transform global_H_local_middle, global_H_local_0, global_H_local_1;

		this->m_restShape.get()->getRestTransformOnX(global_H_local_middle, x_middle);
		this->m_restShape.get()->getRestTransformOnX(global_H_local_0, curvAbs.x());
		this->m_restShape.get()->getRestTransformOnX(global_H_local_1, curvAbs.y());

		local_H_local0_rest = global_H_local_middle.inversed() * global_H_local_0;
		local_H_local1_rest = global_H_local_middle.inversed() * global_H_local_1;

		return;
	}

	this->f_printLog.setValue(true);
	serr << "WARNING : getRestTransform not implemented for not straightRestShape" << sendl;


	// the beam is straight: local is in the middle of local0 and local1
	// the transformation between local0 and local1 is provided by the length of the beam
	double edgeMidLength = this->m_lengthList.getValue()[edgeInList] / 2.0;

	local_H_local0_rest.set(-Vec3(edgeMidLength,0,0), Quat());
	local_H_local1_rest.set(Vec3(edgeMidLength,0,0), Quat());
}



template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getBeamAtCurvAbs(const Real& x_input, unsigned int &edgeInList_output, Real& baryCoord_output, unsigned int start)
{
	////lTotalRest = total length of the
	//Real lTotalRest = this->getRestTotalLength();
	////LTotal =  length sum of the beams that are "out"
	//Real LTotal=0.0;
	//unsigned int start=0;


	if(this->brokenInTwo )
	{
		Real x_abs_broken = this->m_restShape.get()->getReleaseCurvAbs();

		////////// case 1.a : broken part !!
		if (x_input > x_abs_broken)
		{

			// x_i = curv_abs from the begining of the broken part
			Real x_i = x_input-x_abs_broken;
			Real x=0.0;

			for (unsigned int e=0; e<this->_numBeamsNotUnderControl; e++)
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

			edgeInList_output = this->_numBeamsNotUnderControl-1;
			baryCoord_output = 1.0;
			return;

		}
		////////// case 1.b : controlled part !!
		else
		{
			start = this->_numBeamsNotUnderControl;
		}
	}

	Inherited::getBeamAtCurvAbs(x_input, edgeInList_output, baryCoord_output, start);

	//// we find the length of the beam that is "out"
	//for (unsigned int e=start; e<this->Edge_List.size(); e++)
	//{
	//	LTotal += this->getLength(e);
	//}


	//// x_i = abs_curv from the begining of the instrument
	//Real  x_i = x_input + LTotal - lTotalRest;

	//if( x_i < 0.0)
	//{
	//	edgeInList_output = start;
	//	baryCoord_output = 0;
	//	return;
	//}

	//// we compute the x value of each node (the topology is supposed to be a regular seq of segment
	//Real x = 0;

	//for (unsigned int e=start; e<this->Edge_List.size(); e++)
	//{
	//	x += this->getLength(e);
	//	if(x > x_i)
	//	{
	//		edgeInList_output = e;
	//		Real x0 = x - this->getLength(e);
	//		baryCoord_output =(x_i-x0) / this->getLength(e);
	//		return;
	//	}
	//}

	//edgeInList_output = this->Edge_List.size()-1;
	//baryCoord_output = 1.0;
}

template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getCurvAbsAtBeam(const unsigned int &edgeInList_input, const Real& baryCoord_input, Real& x_output)
{	// TODO : version plus complete prenant en compte les coupures et autres particularites de ce modele ?
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
	ProjectionSearch<DataTypes> ps(this, x_input, vecX, xcurv_output, tolerance);
	return ps.doSearch(xcurv_output);
	
/*	unsigned int edge;
	Real bx;
	
	if(xcurv_output<0)
	{
		edge=0;
		bx=0;
	}
	else if(xcurv_output > this->getRestTotalLength())
	{
		edge=this->getNumBeams()-1;
		bx=1;
	}
	else
		getBeamAtCurvAbs(xcurv_output, edge, bx);

	Vec3 P0,P1,P2,P3;
	unsigned int it=0;
	Transform global_H_local0, global_H_local1;
	bool lastTry = false;

    while(it < this->getNumBeams()+1 ) // no more iteration than the number of beam....
    {
        this->getSplinePoints(edge, vecX, P0,P1,P2,P3);

        Vec3 P = P0*(1-bx)*(1-bx)*(1-bx) + P1*3*bx*(1-bx)*(1-bx) + P2*3*bx*bx*(1-bx) + P3*bx*bx*bx;
        Vec3 dP= P0*(-3*(1-bx)*(1-bx)) + P1*(3-12*bx+9*bx*bx) + P2*(6*bx-9*bx*bx) + P3*(3*bx*bx);
        Real f_x = dot( (x_input - P) , dP ) ;

        if (fabs(f_x)/dP.norm() < tolerance) // reach convergence
        {
            getCurvAbsAtBeam(edge,bx,xcurv_output);
            return true;
        }
        Vec3 d2P = P0*6*(1-bx) + P1*(-12 + 18*bx) + P2*(6-18*bx) + P3*6*bx;

        Real df_x = dot(-dP,dP) + dot((x_input-P), d2P);

        // debug
        getCurvAbsAtBeam(edge,bx,xcurv_output);
//        std::cout<<" test at xcurv ="<<xcurv_output<<"  f_x = "<< f_x<<" df_x ="<<df_x<<"  x_input-P ="<<x_input-P<<std::endl;

        if (fabs(df_x) < 1e-5*tolerance)
        {
//            serr<<"Problem in getCurvAbsOfProjection : local minimum without solution and f_x= "<<f_x<<" not null"<<sendl;
            continue;
        }

        Real d_bx = -f_x/df_x;

        if(bx+d_bx > 1)
        {
            if (edge == this->getNumBeams()-1)
            {
				if(lastTry) // Did not find a solution
					return false;
//                serr<<" Problem: no solution found on the thread..."<<sendl;
                // try a last iteration at the end of the thread...
				lastTry = true;
                it = this->getNumBeams()-1;
                edge=this->getNumBeams()-1;
                bx=1;
                xcurv_output=this->getRestTotalLength();
            }
            else
            {
                edge++;
                bx = _min_(bx+d_bx-1.0,1.0);
            }
        }
		else if(bx+d_bx< 0)
        {
            if (edge == 0)
            {
				if(lastTry) // Did not find a solution
					return false;
//                serr<<" Problem: no solution found on the thread..."<<sendl;
                // try a last iteration at the begining of the thread...
				lastTry = true;
                it = this->getNumBeams()-1;
                edge=0;
                bx=0.0;
                xcurv_output=0.0;
            }
            else
            {
                edge--;
                bx = _max_(bx+d_bx+1.0,0.0);
            }
        }
		else
			bx+=d_bx;

        it++;
    }

	return true;	// */
}

template<class DataTypes>
bool WireBeamInterpolation<DataTypes>::breaksInTwo(const Real &x_min_out,  Real &x_break, int &numBeamsNotUnderControlled )
{
    const Real eps = 0.0000000001;

    if (this->brokenInTwo)
    {
        serr << " already broken" << sendl;
        return false;
    }

    if (!this->isControlled() || this->m_restShape == NULL || x_min_out <= eps)
    {
        serr<<" problem with function breaksInTwo "<<sendl;
        return false;
    }

    // if the release point is not "out" (x_min_out> x_break), then the break is not possible
    x_break = m_restShape.get()->getReleaseCurvAbs();
    if (x_min_out > x_break)
        return false;


    // put the info of the "released" part of the beam in the beginning of the beams;
    this->_numBeamsNotUnderControl=0;
    unsigned int duplicatePoint=0;

    // browse the curvilinear abscissa to find the point that needs to be duplicate
    // put the info of the second part of the wire at the begining
    unsigned int i=0;

	VecElementID &edgeList = *this->m_edgeList.beginEdit();
	vector< double > &lengthList = *this->m_lengthList.beginEdit();
	vector< Transform > &DOF0TransformNode0 = *this->m_DOF0TransformNode0.beginEdit();
	vector< Transform > &DOF1TransformNode1 = *this->m_DOF1TransformNode1.beginEdit();
    vector< Vec2 > &curvAbsList = *this->m_curvAbsList.beginEdit();

	const unsigned int curvAbsListSize = curvAbsList.size();

    for (unsigned int e = 1; e < curvAbsListSize; e++)
    {
        if (fabs(curvAbsList[e].x() - x_break) < eps)
        {
            duplicatePoint = e;
            this->_numBeamsNotUnderControl = curvAbsListSize - e;
        }

        if (curvAbsList[e].x() > (x_break - eps))
        {
            edgeList[i] = edgeList[e];
            lengthList[i] = lengthList[e];
            curvAbsList[i] = curvAbsList[e];

            // When the instrument are rotated we apply a transformation between DOF and beam node
            // (should always be the case and dofsAndBeamsAligned should be false)
            if (!this->dofsAndBeamsAligned.getValue())
            {
                DOF0TransformNode0[i] = DOF0TransformNode0[e];
                DOF1TransformNode1[i] = DOF1TransformNode1[e];
            }

            i++;
        }
    }

	this->m_edgeList.endEdit();
	this->m_lengthList.endEdit();
	this->m_DOF0TransformNode0.endEdit();
	this->m_DOF1TransformNode1.endEdit();
	this->m_curvAbsList.endEdit();

    if (duplicatePoint == 0)
        serr << " Problem no point were found at the x_break position ! getReleaseCurvAbs() should provide a <<notable>> point" << sendl;

    m_restShape.get()->releaseWirePart();
    numBeamsNotUnderControlled = this->_numBeamsNotUnderControl;

    this->brokenInTwo = true;

    return true;
}

template<class DataTypes>
bool ProjectionSearch<DataTypes>::doSearch(Real& result)
{
	initSearch(e);

	// Do one pass of the newton method
	newtonMethod();

	Real dist = computeDistAtCurvAbs(e);

	// If the new estimate is good, go with it
	if(found)
	{
		result = e;
		return testForProjection(result);
	}

	// Look if the new estimate is closer to the target than the previous estimate
	if(dist < de)
	{
		// If it is, continue with the newton method, it should converge fast
		while(totalIterations < interpolation->getNumBeams()+1)
		{
			de = dist;
			newtonMethod();
			totalIterations++;

			if(found)
			{
				// Go back to global curv abs
				result = beamStart + (beamEnd - beamStart) * le;
                return testForProjection(result);
			}

			dist = computeDistAtCurvAbs(e);
			if(dist > de)	// Diverging
				break;

			if(le < 0 || le > 1)	// We have to change beam, use the dichotomic method instead
				break;
		}
	}

	// If the estimate is outside the beam, or is further than the previous one, change beam
	//  and use a dichotomic search until we find a solution
	if(!found)
	{
		totalIterations = 0;

		while(totalIterations < interpolation->getNumBeams() + 10 && dichotomicIterations < 10)
		{
			totalIterations++;
			dichotomicIterations++;

			// We will compute 10 samples
			range = segEnd - segStart;
			rangeSampling = range / sampling;
			for(unsigned int i=0; i<=sampling; i++)
			{
				Real curvAbs = segStart + rangeSampling * (Real)i;
				distTab[i] = computeDistAtCurvAbs(curvAbs);
			}
			unsigned int minIndex = std::min_element(distTab, distTab + sampling + 1) - distTab;	// Where is the minimum
			e = segStart + rangeSampling * minIndex;
			if(testForProjection(e))
			{
				result = e;
				return true;
			}

			// If the minium is at one extremity, change beam
			if(minIndex == 0)
				changeCurrentBeam(beamIndex - 1);
			else if(minIndex == sampling)
				changeCurrentBeam(beamIndex + 1);
			else
			{	// We continue the search with a smaller interval (keeping only 3 points)
				segStart = e - rangeSampling;
				segEnd = e + rangeSampling;
			}
		}
	}

	result = e;
	return testForProjection(result);
}

template<class DataTypes>
void ProjectionSearch<DataTypes>::initSearch(Real curvAbs)
{
	e = curvAbs;
	if(e < 0)
	{
		beamIndex = 0;
		le = 0;
	}
	else if(e > interpolation->getRestTotalLength())
	{
		beamIndex = interpolation->getNumBeams()-1;
		le = 1;
	}
	else
		interpolation->getBeamAtCurvAbs(e, beamIndex, le);

	searchDirection = 0;
	interpolation->getAbsCurvXFromBeam(beamIndex, beamStart, beamEnd);
	segStart = beamStart;
	segEnd = beamEnd;
	interpolation->getSplinePoints(beamIndex, x, P0, P1, P2, P3);
	de = computeDistAtCurvAbs(e);
}

template<class DataTypes>
void ProjectionSearch<DataTypes>::newtonMethod()
{
	Real bx = le, bx2 = bx*bx, bx3 = bx2*bx;
	Real obx = 1-bx, obx2 = obx*obx, obx3 = obx2*obx;
    Vec3 P = P0*obx3 + P1*3*bx*obx2 + P2*3*bx2*obx + P3*bx3;
    Vec3 dP = P0*(-3*obx2) + P1*(3-12*bx+9*bx2) + P2*(6*bx-9*bx2) + P3*(3*bx2);
    Real f_x = dot( (target - P) , dP ) ;

	if (f_x==0.0 || fabs(f_x)/dP.norm() < tolerance) // reach convergence
    {
        interpolation->getCurvAbsAtBeam(beamIndex, le, e);
		found = true;
        return;
    }
    Vec3 d2P = P0*6*(1-bx) + P1*(-12 + 18*bx) + P2*(6-18*bx) + P3*6*bx;

    Real df_x = dot(-dP,dP) + dot((target-P), d2P);

    if (fabs(df_x) < 1e-5*tolerance)
    {
//		std::cerr << "Problem in getCurvAbsOfProjection : local minimum without solution and f_x= "<<f_x<<" not null"<< std::endl;
        return;
    }

    Real d_bx = -f_x/df_x;
	le += d_bx;

	e = beamStart + (beamEnd - beamStart) * le;

	// NOTE : bx+d_bx-1.0 ne donne pas une estimation correcte de la position dans l'autre beam, puisque sa longueur peut �tre diff�rente !
}

template<class DataTypes>
bool ProjectionSearch<DataTypes>::changeCurrentBeam(int index)
{
	// If at the end of the thread
	if(index < 0)
	{
		segEnd = segStart + rangeSampling;
		return false;
	}
	else if(index > static_cast<int>(interpolation->getNumBeams())-1)
	{
		segStart = segEnd - rangeSampling;
		return false;
	}

	int dir = index - beamIndex;
	if(searchDirection * dir < 0)
	{	// Changing the direction of search means we are looking for a point near an extremity
		// We know we are looking for a point inside the interval [0.9;1.1]
		// but the ranges for each beam can be different
		Real nStart, nEnd;
		interpolation->getAbsCurvXFromBeam(index, nStart, nEnd);
		Real nRangeSampling = (nEnd - nStart) / sampling;
	
		if(dir < 0)
		{
			segStart = beamStart - nRangeSampling;
			segEnd = beamStart + rangeSampling;
		}
		else
		{
			segStart = beamEnd - rangeSampling;
			segEnd = beamEnd + nRangeSampling;
		}
		return false;
	}
	

	if(dir < 0) 
		searchDirection = -1;
	else if(dir > 0)
		searchDirection = 1;
	
	// Really changing beam
	beamIndex = index;
	interpolation->getAbsCurvXFromBeam(beamIndex, beamStart, beamEnd);
	segStart = beamStart;
	segEnd = beamEnd;
	interpolation->getSplinePoints(beamIndex, x, P0, P1, P2, P3);
	dichotomicIterations = 0;

	return true;
}

template<class DataTypes>
typename ProjectionSearch<DataTypes>::Real ProjectionSearch<DataTypes>::computeDistAtCurvAbs(Real curvAbs)
{
	if(curvAbs >= beamStart && curvAbs <= beamEnd)
	{	// We can use the control points we saved
		Real bx = (curvAbs - beamStart) / (beamEnd - beamStart), bx2 = bx*bx, bx3 = bx2*bx;
		Real obx = 1-bx, obx2 = obx*obx, obx3 = obx2*obx;
		Vec3 P = P0*obx3 + P1*3*bx*obx2 + P2*3*bx2*obx + P3*bx3;

		return (target-P).norm();
	}
	else
	{	// TODO : save all the control points so we don't have to compute them again
		Real bx;
		unsigned int index;
		Vec3 tP0, tP1, tP2, tP3;
		interpolation->getBeamAtCurvAbs(curvAbs, index, bx);
		interpolation->getSplinePoints(index, x, tP0, tP1, tP2, tP3);
		Real bx2 = bx*bx, bx3 = bx2*bx;
		Real obx = 1-bx, obx2 = obx*obx, obx3 = obx2*obx;
		Vec3 P = tP0*obx3 + tP1*3*bx*obx2 + tP2*3*bx2*obx + tP3*bx3;

		return (target-P).norm();
	}
}

template<class DataTypes>
bool ProjectionSearch<DataTypes>::testForProjection(Real curvAbs)
{
	Real bx;
	unsigned int index;
	Vec3 tP0, tP1, tP2, tP3;
	interpolation->getBeamAtCurvAbs(curvAbs, index, bx);
	interpolation->getSplinePoints(index, x, tP0, tP1, tP2, tP3);
	Real bx2 = bx*bx, bx3 = bx2*bx;
	Real obx = 1-bx, obx2 = obx*obx, obx3 = obx2*obx;
	Vec3 P = tP0*obx3 + tP1*3*bx*obx2 + tP2*3*bx2*obx + tP3*bx3;
    Vec3 dP = tP0*(-3*obx2) + tP1*(3-12*bx+9*bx2) + tP2*(6*bx-9*bx2) + tP3*(3*bx2);
    Real f_x = dot( (target - P) , dP ) ;

    if (fabs(f_x)/dP.norm() < tolerance)
		return true;

	return false;
}

} // namespace fem

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_FEM_WIREBEAMINTERPOLATION_INL */
