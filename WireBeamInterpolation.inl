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

#include <sofa/core/behavior/ForceField.inl>
#include <sofa/core/topology/BaseMeshTopology.h>
//#include <sofa/component/topology/GridTopology.h>
//#include <sofa/simulation/common/Simulation.h>
#include <sofa/helper/PolarDecompose.h>
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

    if(!this->isControlled())
        serr<<" not straightRestShape or this->Edge_List is assigned "<<sendl;
    else
        sout<<"external controller for this ForceField is detected"<<sendl;
}




template <class DataTypes>
void WireBeamInterpolation<DataTypes>::getAbsCurvXFromBeam(int& beam, Real& x_curv)
{
    //std::cout<<" getAbsCurvXFromBeam ("<<beam<<") - num Beam in Curv_abs_List: "<<Curv_abs_List.size()<<std::endl;
    x_curv = this->Curv_abs_List[beam].second;
    return;
}



template<class DataTypes>
void WireBeamInterpolation<DataTypes>::clear()
{
    if(this->brokenInTwo)
    {
        this->Edge_List.resize(this->_numBeamsNotUnderControl);
        this->Length_List.resize(this->_numBeamsNotUnderControl);
        this->Curv_abs_List.resize(this->_numBeamsNotUnderControl);
        this->DOF0_Transform_node0.resize(this->_numBeamsNotUnderControl);
        this->DOF1_Transform_node1.resize(this->_numBeamsNotUnderControl);
    }
    else
    {
        this->Edge_List.clear();
        this->Length_List.clear();
        this->Curv_abs_List.clear();
        this->DOF0_Transform_node0.clear();
        this->DOF1_Transform_node1.clear();
    }
}

template<class DataTypes>
void WireBeamInterpolation<DataTypes>::addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1, const Real &angle)
{

    this->Edge_List.push_back(eID);
    this->Length_List.push_back(length);
    std::pair<Real,Real> x_pair(x0,x1);
    this->Curv_abs_List.push_back(x_pair);

    Quat QuatX ;
    QuatX.axisToQuat(Vec3(1,0,0), angle);
    QuatX.normalize();

    // as a angle is set between DOFs and Beam, they are no more aligned
    this->dofsAndBeamsAligned.setValue(false);
    this->DOF0_Transform_node0.push_back(Transform(Vec3(0,0,0), QuatX ));
    this->DOF1_Transform_node1.push_back(Transform(Vec3(0,0,0), QuatX ));
}

template<class DataTypes>
void WireBeamInterpolation<DataTypes>::addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1,
                                               const Transform &DOF0_H_Node0, const Transform &DOF1_H_Node1)
{

    this->Edge_List.push_back(eID);
    this->Length_List.push_back(length);
    std::pair<Real,Real> x_pair(x0,x1);
    this->Curv_abs_List.push_back(x_pair);

    // as a angle is set between DOFs and Beam, they are no more aligned
    this->dofsAndBeamsAligned.setValue(false);
    this->DOF0_Transform_node0.push_back(DOF0_H_Node0);
    this->DOF1_Transform_node1.push_back(DOF1_H_Node1);
}



template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getRestTransform(unsigned int edgeInList, Transform &local0_H_local1_rest)
{
	this->f_printLog.setValue(true);
	serr<<"WARNING : getRestTransform not implemented for not straightRestShape"<<sendl;

	// the beam is straight: the transformation between local0 and local1 is provided by the length of the beam
	local0_H_local1_rest.set(Vec3(this->Length_List[edgeInList],0,0), Quat());
}


template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest)
{

	if(this->isControlled() && this->m_restShape!=NULL)
	{
		Real x_middle = (this->Curv_abs_List[edgeInList].first + this->Curv_abs_List[edgeInList].second)/2;
		Transform global_H_local_middle, global_H_local_0, global_H_local_1;
                this->m_restShape.get()->getRestTransformOnX(global_H_local_middle,x_middle );
                this->m_restShape.get()->getRestTransformOnX(global_H_local_0,this->Curv_abs_List[edgeInList].first );
                this->m_restShape.get()->getRestTransformOnX(global_H_local_1,this->Curv_abs_List[edgeInList].second );
		local_H_local0_rest = global_H_local_middle.inversed()*global_H_local_0;
		local_H_local1_rest = global_H_local_middle.inversed()*global_H_local_1;

		return;
	}

	this->f_printLog.setValue(true);
	serr<<"WARNING : getRestTransform not implemented for not straightRestShape"<<sendl;


	// the beam is straight: local is in the middle of local0 and local1
	// the transformation between local0 and local1 is provided by the length of the beam
	local_H_local0_rest.set(-Vec3(this->Length_List[edgeInList]/2,0,0), Quat());
	local_H_local1_rest.set(Vec3(this->Length_List[edgeInList]/2,0,0), Quat());

}



template<class DataTypes>
 void WireBeamInterpolation<DataTypes>::getBeamAtCurvAbs(const Real& x_input, unsigned int &edgeInList_output, Real& baryCoord_output)
 {
     //lTotalRest = total length of the
     Real lTotalRest = this->getRestTotalLength();
     //LTotal =  length sum of the beams that are "out"
     Real LTotal=0.0;
     unsigned int start=0;


     if(this->brokenInTwo )
     {

         Real x_abs_broken = this->m_restShape.get()->getReleaseCurvAbs();

         if (x_input > x_abs_broken)         ////////// case 1.a : broken part !!
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
         else                                   ////////// case 1.b : controlled part !!
         {

             start = this->_numBeamsNotUnderControl;

         }

     }

     // we find the length of the beam that is "out"
     for (unsigned int e=start; e<this->Edge_List.size(); e++)
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

     // we compute the x value of each node (the topology is supposed to be a regular seq of segment
     Real x = 0;

     for (unsigned int e=start; e<this->Edge_List.size(); e++)
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

     edgeInList_output = this->Edge_List.size()-1;
     baryCoord_output = 1.0;

}

template<class DataTypes>
void WireBeamInterpolation<DataTypes>::getCurvAbsAtBeam(unsigned int &edgeInList_input, Real& baryCoord_input, Real& x_output)
{	// TODO : version plus complete prenant en compte les coupures et autres particularites de ce modele ?
	x_output = 0;
	for(unsigned int i=0; i<edgeInList_input; i++)
		x_output += this->getLength(i);

	x_output += this->getLength(edgeInList_input) * baryCoord_input;
}

template<class DataTypes>
bool WireBeamInterpolation<DataTypes>::getApproximateCurvAbs(const Vec3& x_input, const VecCoord& x, Real& x_output)
{
	unsigned int beamIndex = 0;
	int dir = 0;
        Real beamBary = 0.0, closestDist = -1;
	bool projected = false;
	
	if(x.size() <= 1) 
	{ 
		x_output = 0.0; 
		return false; 
	}

	unsigned int nb = x.size() - 1;
	for(unsigned int i=0; i<nb; i++) // Let say we go from "left" to "right"
	{
		Transform globalHlocal0, globalHlocal1;
		computeTransform2(i, globalHlocal0, globalHlocal1, x);
		Vec3 A = globalHlocal0.getOrigin(), B = globalHlocal1.getOrigin();
		Real r = ((x_input-A) * (B-A)) / (B-A).norm2();
		
		if(r < 0)	// "left" of the segment
		{
			if(dir == 1) // We were on the right of the previous segment, test with the vertex itself
			{
				double dist = (x_input-x[i].getCenter()).norm2();
				if(closestDist<0 || dist<closestDist)
				{
					beamIndex = i;
					beamBary = 0.0;
					projected = true;	// This is a correct answer, even if it lie between 2 segments
					closestDist = dist;
				}
			}
			dir = -1;
		}
		else if(r > 1)	// "right" of the segment
		{
			dir = 1;
		}
		else if(r<=1) // on the segment
		{
			Vec3 proj = A + (B-A) * r;
			double dist = (x_input-proj).norm2();
			if(closestDist<0 || dist<closestDist)
			{
				beamIndex = i;
				beamBary = r;
				projected = true;
				closestDist = dist;
			}
		}
	}

	if(!projected)	// Not on the curve, test start and end vertices
	{
                Real d1 = (x_input-x[0].getCenter()).norm2();
                Real d2 = (x_input-x[nb].getCenter()).norm2();
		if(d1 < d2)
			x_output = 0.0;
		else
			x_output = this->getRestTotalLength();
		return false;
	}

	// We know the beam the point can be projected to, translate that to an abscissa
	getCurvAbsAtBeam(beamIndex, beamBary, x_output);
	return true;
}




template<class DataTypes>
bool WireBeamInterpolation<DataTypes>::getCurvAbsOfProjection(const Vec3& x_input, const VecCoord& x, Real& xcurv_output, const Real& tolerance)
{
    unsigned int edge;
    Real bx;

    if (xcurv_output>0 && xcurv_output<  this->getRestTotalLength())
        getBeamAtCurvAbs(xcurv_output, edge, bx);
    else
    {
        edge=0;
        bx=0;
    }


    Vec3 P0,P1,P2,P3;
    unsigned int it=0;
    Transform global_H_local0, global_H_local1;
	bool lastTry = false;

    while(it < this->getNumBeams()+1 ) // no more iteration than the number of beam....
    {
        this->getSplinePoints(edge, x, P0,P1,P2,P3);

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
            serr<<"Problem in getCurvAbsOfProjection : local minimum without solution and f_x= "<<f_x<<" not null"<<sendl;
            continue;
        }

        Real d_bx = -f_x/df_x;

        if(bx+d_bx > 1)
        {
            if (edge == this->getNumBeams()-1)
            {
				if(lastTry) // Did not find a solution
					return false;
                serr<<" Problem: no solution found on the thread..."<<sendl;
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
                serr<<" Problem: no solution found on the thread..."<<sendl;
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

	return true;
}


template<class DataTypes>
bool WireBeamInterpolation<DataTypes>::breaksInTwo(const Real &x_min_out,  Real &x_break, int &numBeamsNotUnderControlled )
{

    Real eps=0.0000000001;

    if(this->brokenInTwo)
    {
        serr<<" already broken"<<sendl;
        return false;
    }

    if( !this->isControlled() || this->m_restShape==NULL || x_min_out <= eps)
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


    for (unsigned int e=1; e<this->Curv_abs_List.size(); e++)
    {
        if( fabs(this->Curv_abs_List[e].first-x_break)<eps)
        {
            duplicatePoint = e;
            this->_numBeamsNotUnderControl = this->Curv_abs_List.size() - e;
        }

        if(this->Curv_abs_List[e].first > (x_break-eps) )
        {
            this->Edge_List[i] = this->Edge_List[e];
            this->Length_List[i] =  this->Length_List[e];
            this->Curv_abs_List[i]= this->Curv_abs_List[e];


            // When the instrument are rotated we apply a transformation between DOF and beam node
            // (should always be the case and dofsAndBeamsAligned should be false)
            if(!this->dofsAndBeamsAligned.getValue())
            {
                this->DOF0_Transform_node0[i] = this->DOF0_Transform_node0[e];
                this->DOF1_Transform_node1[i] = this->DOF1_Transform_node1[e];

            }

            i++;
        }


    }

    if(duplicatePoint==0)
        serr<<" Problem no point were found at the x_break position ! getReleaseCurvAbs() should provide a <<notable>> point"<<sendl;

    m_restShape.get()->releaseWirePart();
    numBeamsNotUnderControlled = this->_numBeamsNotUnderControl;

    this->brokenInTwo = true;

    return true;


}



}// namespace fem

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_FEM_WIREBEAMINTERPOLATION_INL */
