/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMLENGTHCONSTRAINT_INL
#define SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMLENGTHCONSTRAINT_INL

#include "AdaptiveBeamLengthConstraint.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/gl/template.h>
#include <sofa/core/visual/DrawTool.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

void AdaptiveBeamLengthConstraintResolution::init(int line, double** /*w*/, double* force) 
{ 
	if(_initF)
		force[line] = *_initF;
}
void AdaptiveBeamLengthConstraintResolution::resolution(int line, double** w, double* d, double* force)
{
    force[line] -= d[line] / w[line][line];
    if(force[line] < 0)
            force[line] = 0;
}

void AdaptiveBeamLengthConstraintResolution::store(int line, double* force, bool /*convergence*/)
{
	if(_initF)
		*_initF = force[line];
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::init()
{
    this->mstate= dynamic_cast< core::behavior::MechanicalState<DataTypes> *> (this->getContext()->getMechanicalState());
        assert(this->mstate);
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::reset()
{
	internalInit();
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::internalInit()
{	// We search for the closest segment, on which to project each point
	// Convention : object1 is the beam model, object2 is the list of point constraints

	if(!m_interpolation.get())
	{
		serr << "Could not find the beam interpolation" << sout;
		return;
	}

}


template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::detectElongation(const VecCoord& x, const VecCoord& xfree)
{
	Vec3 P0,P1,P2,P3;
	Real length, rest_length, length_interval, rest_length_interval;
	Real alarmLength = m_alarmLength.getValue();
	bool prev_stretch = false;

	fem::WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();

	// storage of the length (and the rest_length)  of the interval being stretched
	length_interval=0.0;
	rest_length_interval=0.0;

	// storage of the interval information
	IntervalDefinition intervalDef;
	Real angleInterval=0.0;

	for (unsigned int b=0; b<interpolation->getNumBeams(); b++)
	{
		// 1. compute the actual length and the rest length of the beams
		interpolation->getSplinePoints(b,x,P0,P1,P2,P3);
		//interpolation->computeActualLength(length, P0,P1,P2,P3);
		// TODO : optimization: finally it is not necessary to compute all the spline points
		length=(P0-P3).norm();
		rest_length = interpolation->getLength(b);

		// 2. compute the bending angle
		Transform Tnode0, Tnode1;
		Real angleBeam=0.0;
		interpolation->computeTransform2(b,Tnode0,Tnode1,x);
		interpolation->ComputeTotalBendingRotationAngle(angleBeam, rest_length/10.0, Tnode0, Tnode1,rest_length , 0.0, 1.0);


		// 3. treatment of the different case..
		unsigned n0, n1;
		interpolation->getNodeIndices(b, n0, n1);
		bool case1a = (n0==n1);

		if(prev_stretch)
		{
			////CASE 1: previous beam was stretched:
			// find the case that necessitates to stop the interval:
			// (a) rigidification
			// (b) current beam not stretched
			// (c) too large bending angle...
			// (d) last beam ! [ see after the loop "for" ]
			//// => store the information : index [begin end], Pos [begin end]

			bool case1b = (rest_length*alarmLength > length);
			bool case1c = (angleBeam >= 0.99*m_maxBendingAngle.getValue());// angleBeam+angleInterval > m_maxBendingAngle.getValue()
			case1c=false;
			case1b=false;

			if (case1a || case1b || case1c) // CASE 1 (a) + (b) + (c)
			{
				if(this->f_printLog.getValue())
				{
					std::cout<<" beam "<<b<<" case 1 detected (a):"<<n0<<" == ?"<<n1<<"  (b) :"<<
					rest_length*alarmLength<<" > ?"<<length<<"  (c) : "<<angleBeam+angleInterval<< " > ?" <<m_maxBendingAngle.getValue() <<std::endl;
				}

				Real angleBuf= angleBeam+angleInterval;

				// Stop the rigidification
				// the interval ends at the beginning of the current beam
				intervalDef.posEnd = P0; // store the position
				intervalDef.IdxEnd = n0; // store the index of the dof

				Transform DOF0_H_local0, DOF1_H_local1;
				interpolation->getDOFtoLocalTransform(b, DOF0_H_local0,  DOF1_H_local1);

				intervalDef.dof_H_end = DOF0_H_local0; // store the transform from dof to pos

				interpolation->getSplinePoints(b,xfree,P0,P1,P2,P3);

				Transform global_H_local0_free, global_H_local1_free;
				interpolation->computeTransform2(b,  global_H_local0_free,  global_H_local1_free, xfree);
				intervalDef.posFreeEnd  = global_H_local0_free.getOrigin(); // store the free position

				intervalDef.rest_length=rest_length_interval; // store the rest_length

				if(this->f_printLog.getValue())
				std::cout<<" rest_length_interval ="<<rest_length_interval<<std::endl;

				// ends the interval
				prev_stretch=false;
				angleInterval=0.0;
				rest_length_interval=0.0;

				// verify that the interval length is not null:
				if ((intervalDef.posBegin-intervalDef.posEnd).norm() < 0.0000001)
					serr<<"WARNING interval of size = 0 detected => 1 beam has a bendingAngle > m_maxBendingAngle"<<sendl;
				else
					_constraintIntervals.push_back(intervalDef);

				// isolate the case 1 (c)
				if (!case1a && !case1b && case1c) //
				{
					if(this->f_printLog.getValue())
						std::cout<<" isolate case 1(c)"<<std::endl;

					// create a new interval:
					prev_stretch=true;
					angleInterval = angleBeam;
					rest_length_interval = rest_length;

					// store the information of the beginning of the new interval
					// -> it corresponds to the end of the previous interval
					intervalDef.dof_H_begin = DOF0_H_local0;
					intervalDef.IdxBegin = n0;
					intervalDef.posBegin = intervalDef.posEnd;
					intervalDef.posFreeBegin = intervalDef.posFreeEnd;
				}
			}
			else
			{
				// Continue the rigidification
				angleInterval+=angleBeam;
				rest_length_interval += rest_length;
			}
		}
		else
		{
			//// CASE 2: previous beam was not stretched:
			// (a) the current beam is stretched: start the interval => store index_begin Pos_begin
			//     veriy the bending angle (in case the interval= the beam)
			// (b) the current beam is not stretched=> nothing to do !

			bool case2a = (rest_length*alarmLength < length);
			case2a=true;

			if (  case2a && !case1a ) // CASE 2 (a)
			{
				if( this->f_printLog.getValue())
					std::cout<<" beam "<<b<<" case 2 (a) detected "<<rest_length*alarmLength<<" < ?"<<length<<std::endl;
				// create a new interval:
				prev_stretch=true;
				angleInterval = angleBeam;
				rest_length_interval = rest_length;

				// store the information of the beginning of the new interval
				// -> it corresponds to the position of node 0 of the beam
				Transform DOF0_H_local0, DOF1_H_local1;
				interpolation->getDOFtoLocalTransform(b, DOF0_H_local0,  DOF1_H_local1);

				Transform global_H_local0, global_H_local1;
				interpolation->computeTransform2(b,  global_H_local0,  global_H_local1, x);


				Transform global_H_local0_free, global_H_local1_free;
				interpolation->computeTransform2(b,  global_H_local0_free,  global_H_local1_free, xfree);


				intervalDef.dof_H_begin = DOF0_H_local0;
				intervalDef.IdxBegin = n0;
				intervalDef.posBegin = global_H_local0.getOrigin();
				intervalDef.posFreeBegin = global_H_local0_free.getOrigin();
			}
			else
			{
				if( this->f_printLog.getValue())
					std::cout<<" beam "<<b<<" case 2 (b) detected "<<rest_length*alarmLength<<" > ?"<<length<<" or n0="<<n0<<" ==? "<<"n1="<<n1 <<std::endl;
			}
		}
	}

	unsigned int b = interpolation->getNumBeams()-1;
	unsigned n0, n1;
	interpolation->getNodeIndices(b, n0, n1);
	if(prev_stretch) // case 1(d)
	{
		if( this->f_printLog.getValue())
			std::cout<<" case 1 (d) detected on the last beam"<<std::endl;

		// store the information of the beginning of the new interval
		// -> it corresponds to the position of node 0 of the beam
		Transform DOF0_H_local0, DOF1_H_local1;
		interpolation->getDOFtoLocalTransform(b, DOF0_H_local0,  DOF1_H_local1);

		Transform global_H_local0, global_H_local1;

		interpolation->computeTransform2(b,  global_H_local0,  global_H_local1, x);


		Transform global_H_local0_free, global_H_local1_free;
		interpolation->computeTransform2(b,  global_H_local0_free,  global_H_local1_free, xfree);


		intervalDef.dof_H_end = DOF1_H_local1;
		intervalDef.IdxEnd = n1;
		intervalDef.posEnd = global_H_local1.getOrigin();
		intervalDef.posFreeEnd = global_H_local1_free.getOrigin();
		intervalDef.rest_length=rest_length_interval; // store the rest_length

		if(this->f_printLog.getValue())
		std::cout<<" rest_length_interval ="<<rest_length_interval<<std::endl;

		_constraintIntervals.push_back(intervalDef);
	}
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams * /*cParams*/ /* PARAMS FIRST */, DataMatrixDeriv &c_d, unsigned int &constraintId, const DataVecCoord &x_d)
{

  //  std::cout<<"begin buildConstraintMatrix "<<std::endl;
	violations.clear();
    Real constrainedLength = m_constrainedLength.getValue();

	nbConstraints = 0;
	cid = constraintId;
    const VecCoord& x= *this->mstate->getX();
    const VecCoord& xfree = *this->mstate->getXfree();
   // std::cout<< " x ="<<x<<" \n xfree ="<<xfree<<std::endl;

    MatrixDeriv& c = *c_d.beginEdit();

    _constraintIntervals.clear();
    detectElongation( x, xfree);

    if( this->f_printLog.getValue())
    {
        for (unsigned int i=0; i<_constraintIntervals.size();i++)
        {
            std::cout<<"constraint["<<i<<"] between pos: "<<_constraintIntervals[i].posBegin<<" and pos: "<<_constraintIntervals[i].posEnd<<std::endl;

        }
    }


    // for each constraint Interval, a constraint is created
    for (unsigned int i=0; i<_constraintIntervals.size();i++)
    {

        // get the indices of the dofs involved
        unsigned int n0,n1;
        n0=_constraintIntervals[i].IdxBegin;
        n1=_constraintIntervals[i].IdxEnd;


        // compute the violation and the local direction of the constraint
        Vec3 dir = _constraintIntervals[i].posEnd - _constraintIntervals[i].posBegin;
        Vec3 PbPeFree= _constraintIntervals[i].posFreeEnd - _constraintIntervals[i].posFreeBegin;
        Real length = dir.norm();
        dir.normalize();
        Real length_free= dot(dir, PbPeFree);


      //  std::cout<<" dot(dir, PbPeFree) ="<<dot(dir, PbPeFree)<<"  - length_free="<<length_free<<"  - length ="<<length<<std::endl;



        // put the violation in a buffer
        violations.push_back(_constraintIntervals[i].rest_length * constrainedLength - length_free);

        // project the  direction of the constraint (considered as a force) in the DOF frame

        Vec3 lever0 = x[n0].getOrientation().rotate( _constraintIntervals[i].dof_H_begin.getOrigin() );
        Vec3 lever1 = x[n1].getOrientation().rotate( _constraintIntervals[i].dof_H_end.getOrigin()  );


		if(this->f_printLog.getValue())
		{
			if(lever0.norm() > 0.0001)
				std::cout<<" lever0 ="<<lever0<<" dir ="<<dir<<std::endl;

			if(lever1.norm() > 0.0001)
				std::cout<<" lever1 ="<<lever1<<" -dir ="<<-dir<<std::endl;
		}



        MatrixDerivRowIterator c_it = c.writeLine(cid + nbConstraints);

        c_it.addCol(n0, Vec6(dir, cross(lever0, dir) ) );
        c_it.addCol(n1, Vec6(-dir, cross(lever1, -dir)  ) );

        nbConstraints++;

    }
    constraintId +=nbConstraints;
	if(this->f_printLog.getValue())
		std::cout<<"end buildConstraintMatrix "<<std::endl;
}


template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams* /* PARAMS FIRST */, defaulttype::BaseVector *v, const DataVecCoord &, const DataVecDeriv &)
{


	unsigned int nb = violations.size();

	for(unsigned int i=0; i<nb; i++)
    {

		v->set(cid+i, violations[i]);
    }




}


template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{


    unsigned int nb = violations.size();
    for(unsigned int i=0; i<nb; i++)
    {
        resTab[offset] = new AdaptiveBeamLengthConstraintResolution(); //&prevForces[activatedBeamsAbscissa[i]]
        offset++;
    }



}


template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{ 
    if (!vparams->displayFlags().getShowInteractionForceFields()) return;

    if(_constraintIntervals.size()==0)
        return;



    //trace a point at the beginning and at the end of each constraint interval
    glDisable(GL_LIGHTING);
    glPointSize(10);
    glBegin(GL_POINTS);
    glColor4f(0.0f,1.0f,0.0f,1.0f);
    for(unsigned int i=0; i<_constraintIntervals.size(); i++)
    {
        helper::gl::glVertexT(_constraintIntervals[i].posBegin);
        helper::gl::glVertexT(_constraintIntervals[i].posEnd);
    }
    glEnd();

    glPointSize(1);

    // trace a straight line between the beginning and the end of each interval
    // TODO: change color if constrained Length reached...
    glBegin(GL_LINES);
    for(unsigned int i=0; i<_constraintIntervals.size(); i++)
    {
        helper::gl::glVertexT(_constraintIntervals[i].posBegin);
        helper::gl::glVertexT(_constraintIntervals[i].posEnd);
    }

    glEnd();

    glEnable(GL_LIGHTING);




    /*
	if(!vparams->displayFlags().getShowInteractionForceFields()) return;

	glDisable(GL_LIGHTING);
	glPointSize(10);
	glBegin(GL_POINTS);
	unsigned int m = this->mstate2->getSize();
	const VecCoord& x = *this->mstate2->getX();
	for(unsigned int i=0; i<m; i++)
	{
        glColor4f(0.0f,1.0f,projected[i]?1:0.0f,1.0f);
		helper::gl::glVertexT(x[i]);
	}
	
	glEnd();
	glPointSize(1);
     */

}

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
