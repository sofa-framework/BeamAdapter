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
void AdaptiveBeamLengthConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams * /*cParams*/ /* PARAMS FIRST */, DataMatrixDeriv &c_d, unsigned int &constraintId, const DataVecCoord &x_d)
{
	violations.clear();
//	activated_beams.clear();
	activatedBeamsAbscissa.clear();
	nbConstraints = 0;
	cid = constraintId;
    const VecCoord& x= x_d.getValue();
    const VecCoord& xfree = *this->mstate->getXfree();
    MatrixDeriv& c = *c_d.beginEdit();
    Vec3 P0,P1,P2,P3;
    Real length, length_free;
	Real alarmLength = m_alarmLength.getValue(), constrainedLength = m_constrainedLength.getValue();

    fem::WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();

    for (unsigned int b=0; b<interpolation->getNumBeams(); b++)
    {
        Real rest_length = interpolation->getLength(b);
        interpolation->getSplinePoints(b,x,P0,P1,P2,P3);
        interpolation->computeActualLength(length, P0,P1,P2,P3);

        if (length > rest_length * alarmLength)
        {
			unsigned int n0,n1;
			interpolation->getNodeIndices(b, n0, n1);
			if (n0==n1)
				continue; // no constraint => the beam is "rigidified"

//			activated_beams.push_back(b);
			Real abs = 0;
			interpolation->getAbsCurvXFromBeam(b, abs);
			activatedBeamsAbscissa.push_back(abs);

			interpolation->getSplinePoints(b,xfree,P0,P1,P2,P3);
			interpolation->computeActualLength(length_free, P0,P1,P2,P3);

			violations.push_back(rest_length * constrainedLength - length_free);

			Transform global_H_local0, global_H_local1;
			interpolation->computeTransform2(b, global_H_local0, global_H_local1, x);

			Vec3 dir_x0_global = global_H_local0.projectVector( Vec3(1.0,0.0,0.0) );
			Vec3 dir_x1_global = global_H_local1.projectVector( Vec3(1.0,0.0,0.0) );

			MatrixDerivRowIterator c_it = c.writeLine(cid + nbConstraints);

			c_it.addCol(n0, Vec6(dir_x0_global, Vec3(0.0,0.0,0.0) ) );
			c_it.addCol(n1, Vec6(-dir_x1_global, Vec3(0.0,0.0,0.0) ) );

			nbConstraints++;
        }
    }

    constraintId += nbConstraints;
    c_d.endEdit();
}


template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams* /* PARAMS FIRST */, defaulttype::BaseVector *v, const DataVecCoord &, const DataVecDeriv &)
{
	unsigned int nb = violations.size();
	for(unsigned int i=0; i<nb; i++)
		v->set(cid+i, violations[i]);
}


template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
	std::sort(activatedBeamsAbscissa.begin(), activatedBeamsAbscissa.end());
	unsigned int i=0, nb = activatedBeamsAbscissa.size();

	// Removing from the map old forces for beams that are not constrained anymore
	for(std::map<Real, double>::iterator iter=prevForces.begin(); iter!=prevForces.end(); )
	{
		Real key = iter->first;
		// Does this key exist in the beams list ?
		bool exist = false;
		for(unsigned int j=i; j<nb; j++)
			if(activatedBeamsAbscissa[j] == key)
			{
				exist = true;
				i = j;
				break;
			}

		if(!exist)
			prevForces.erase(iter++);
		else
			++iter;
	}

	// Make sure we have a force for each constrained beam
	for(i=0; i<nb; i++)
		prevForces[activatedBeamsAbscissa[i]];	// Will create a new pair if not existant

    nb = violations.size();
    for(unsigned int i=0; i<nb; i++)
    {
        resTab[offset] = new AdaptiveBeamLengthConstraintResolution(&prevForces[activatedBeamsAbscissa[i]]);
        offset++;
    }
}


template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::draw(const core::visual::VisualParams* /*vparams*/)
{ 
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
