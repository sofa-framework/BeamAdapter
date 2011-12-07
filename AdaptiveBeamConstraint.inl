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

	int m2 = this->mstate2->getSize();	
	previousPositions.clear();
	previousPositions.resize(m2);
	projected.clear();
	projected.resize(m2);
	
	fem::WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();
	const VecCoord& x1 = *this->mstate1->getX();
	const VecCoord& x2 = *this->mstate2->getX();

	for(int i=0; i<m2; i++)
	{
		Real r = -1;
		bool p = interpolation->getApproximateCurvAbs(x2[i].getCenter(), x1, r);

		previousPositions[i] = r;
		projected[i] = p;
		std::cout << "Point " << i << (p?" projected at ":" not projected near ") << r << std::endl; 
	}
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams * /*cParams*/ /* PARAMS FIRST */, DataMatrixDeriv &c1_d, DataMatrixDeriv &c2_d, unsigned int &constraintId
		, const DataVecCoord &x1_d, const DataVecCoord &x2_d)
{
}


template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams* /* PARAMS FIRST */, defaulttype::BaseVector *v, const DataVecCoord &, const DataVecCoord &
		, const DataVecDeriv &, const DataVecDeriv &)
{
}


template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{ 
	glDisable(GL_LIGHTING);
	glPointSize(10);
	glBegin(GL_POINTS);
	int m = this->mstate2->getSize();
	const VecCoord& x = *this->mstate2->getX();
	for(int i=0; i<m; i++)
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
