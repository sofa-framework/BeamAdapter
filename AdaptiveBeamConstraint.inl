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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMCONSTRAINT_INL
#define SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMCONSTRAINT_INL

#include "AdaptiveBeamConstraint.h"
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/constraintset/BilateralInteractionConstraint.h>

#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/gl/template.h>
namespace sofa
{

namespace component
{

namespace constraintset
{


void AdaptiveBeamConstraintResolution::resolution(int line, double** w, double* d, double* force)
{
	double f[2];
	f[0] = force[line]; f[1] = force[line+1];
	
	force[line] -= d[line] / w[line][line];
	d[line+1] += w[line+1][line] * (force[line]-f[0]);
	force[line+1] -= d[line+1] / w[line+1][line+1];
	d[line+2] += w[line+2][line] * (force[line]-f[0]) + w[line+2][line+1] * (force[line+1]-f[1]);
}

void AdaptiveBeamConstraintResolution::init(int line, double** w, double* /*force*/)
{
	slidingW = w[line+2][line+2];
}

void AdaptiveBeamConstraintResolution::store(int line, double* force, bool /*convergence*/)
{
	if(_slidingDisp)
		*_slidingDisp = force[line+2] * slidingW;
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::init()
{
	assert(this->mstate1);
	assert(this->mstate2);
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::reset()
{
	internalInit();
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::bwdInit()
{
	internalInit();
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::internalInit()
{	// We search for the closest segment, on which to project each point
	// Convention : object1 is the beam model, object2 is the list of point constraints

	if(!m_interpolation.get())
	{
		serr << "Could not find the beam interpolation" << sout;
		return;
	}

	unsigned int m2 = this->mstate2->getSize();	
	previousPositions.clear();
	previousPositions.resize(m2);
	projected.clear();
	projected.resize(m2);
	displacements.clear();
	displacements.resize(m2);
	
	fem::WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();
	const VecCoord& x1 = *this->mstate1->getX();
	const VecCoord& x2 = *this->mstate2->getX();

	for(unsigned int i=0; i<m2; i++)
	{
		Real r = -1;
		Pos pt = x2[i].getCenter();
		bool p = interpolation->getApproximateCurvAbs(pt, x1, r);

		previousPositions[i] = r;
		projected[i] = p;

		Real r2 = r;
		if(p)
		{
			interpolation->getCurvAbsOfProjection(pt, x1, r2);
			std::cout << "Point " << i << " projected at " << r << " -> " << r2 << std::endl; 
		}
		else
			std::cout << "Point " << i << " not projected near " << r << std::endl; 
	}
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams * /*cParams*/ /* PARAMS FIRST */, DataMatrixDeriv &c1_d, DataMatrixDeriv &c2_d, unsigned int &constraintId
		, const DataVecCoord &x1_d, const DataVecCoord &x2_d)
{
	violations.clear();
	nbConstraints = 0;

    Transform Tnode0, Tnode1, Tresult;
    Real baryCoord;
    unsigned int beam;

	unsigned int m2 = this->mstate2->getSize();
	fem::WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();
	const VecCoord& x1 = *this->mstate1->getXfree();
	const VecCoord& x2 = *this->mstate2->getXfree();
    for(unsigned int i=0; i<m2 ; i++)
    {
		if(!projected[i])
			continue;

		// Get new projection on the curve
		previousPositions[i] += displacements[i];
		interpolation->getCurvAbsOfProjection(x2[i].getCenter(), x1, previousPositions[i]);

		// Position and frame on the curve
        interpolation->getBeamAtCurvAbs(previousPositions[i], beam, baryCoord);
        interpolation->computeTransform2(beam, Tnode0, Tnode1, x1);
        interpolation->InterpolateTransformUsingSpline(Tresult, baryCoord, Tnode0,Tnode1,interpolation->getLength(beam));
		Pos p = Tresult.getOrigin();
		Pos dir, dir1, dir2;
		Rot rot = Tresult.getOrientation();
		dir = rot.rotate(Pos(1,0,0));
		dir1 = rot.rotate(Pos(0,1,0));
		dir2 = rot.rotate(Pos(0,0,1));

		// Compute violations
		Pos violation = x2[i].getCenter() - p;
		violations.push_back(violation * dir);
		violations.push_back(violation * dir1);
		violations.push_back(violation * dir2);

		// Define the constraint
	}
}


template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams* /* PARAMS FIRST */, defaulttype::BaseVector *v, const DataVecCoord &, const DataVecCoord &
		, const DataVecDeriv &, const DataVecDeriv &)
{
/*	unsigned int nb = violations.size();
	for(unsigned int i=0; i<nb; i++)
		v->set(cid+i, violations[i]);	*/
}


template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
/*	unsigned int nb = this->mstate2->getSize();
	for(unsigned int i=0; i<nb; i++)
	{
		if(!projected[i]) continue;
		resTab[offset] = new AdaptiveBeamConstraintResolution(&displacements[i]);
		offset += 3;
	}	*/
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{ 
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
}

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
