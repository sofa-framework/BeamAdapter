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
//#include <Constraint/constraintset/BilateralInteractionConstraint.h>
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
//	force[line+2] = 0;
//	force[line+2] -= d[line+2] / w[line+2][line+2];
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
		Vec3 pt = x2[i].getCenter();
		bool p = interpolation->getApproximateCurvAbs(pt, x1, r);

		previousPositions[i] = r;
		projected[i] = p;
	}
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams * /*cParams*/ /* PARAMS FIRST */, DataMatrixDeriv &c1_d, DataMatrixDeriv &c2_d, unsigned int &constraintId
		, const DataVecCoord &x1_d, const DataVecCoord &x2_d)
{
	violations.clear();
	nbConstraints = 0;
	cid = constraintId;

    Transform Tnode0, Tnode1, Tresult;
    Real baryCoord;
    unsigned int beam;

	unsigned int m2 = this->mstate2->getSize();
	fem::WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();
	MatrixDeriv& c1 = *c1_d.beginEdit();
	MatrixDeriv& c2 = *c2_d.beginEdit();
	const VecCoord& x1= x1_d.getValue();
	const VecCoord& x2= x2_d.getValue();
	const VecCoord& x1free = *this->mstate1->getXfree();
	const VecCoord& x2free = *this->mstate2->getXfree();
    for(unsigned int i=0; i<m2 ; i++)
    {
		if(!projected[i])
			continue;

		// Get new projection on the curve
                previousPositions[i] += (Real) displacements[i];
		if(!interpolation->getCurvAbsOfProjection(x2[i].getCenter(), x1, previousPositions[i], 1e-5))
		{
			projected[i] = false;
			continue;
		}

		// Position and frame on the curve
        interpolation->getBeamAtCurvAbs(previousPositions[i], beam, baryCoord);
        interpolation->computeTransform2(beam, Tnode0, Tnode1, x1free);
        interpolation->InterpolateTransformUsingSpline(Tresult, baryCoord, Tnode0, Tnode1, interpolation->getLength(beam));
		Pos p = Tresult.getOrigin();
		Pos dir, dir1, dir2;
		Rot rot = Tresult.getOrientation();
		dir = rot.rotate(Pos(1,0,0));
		dir1 = rot.rotate(Pos(0,1,0));
		dir2 = rot.rotate(Pos(0,0,1));

		// Compute violations
		Pos violation = p - x2free[i].getCenter();
		violations.push_back(violation * dir1);
		violations.push_back(violation * dir2);
		violations.push_back(violation * dir);

		// Define the constraint
		unsigned int node0, node1;
		SpatialVector sv0, sv1;
		Vec3 nullRot(0,0,0);
		interpolation->getNodeIndices(beam, node0, node1);

		MatrixDerivRowIterator c1_it = c1.writeLine(cid + nbConstraints);
		MatrixDerivRowIterator c2_it = c2.writeLine(cid + nbConstraints);
		interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir1, sv0, sv1);
		c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
		c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
		c2_it.addCol(i, Deriv(-dir1, nullRot));
		nbConstraints++;

		c1_it = c1.writeLine(cid + nbConstraints);
		c2_it = c2.writeLine(cid + nbConstraints);
		interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir2, sv0, sv1);
		c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
		c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
		c2_it.addCol(i, Deriv(-dir2, nullRot));
		nbConstraints++;

		c1_it = c1.writeLine(cid + nbConstraints);
		c2_it = c2.writeLine(cid + nbConstraints);
		interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir, sv0, sv1);
		c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
		c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
		c2_it.addCol(i, Deriv(-dir, nullRot));
		nbConstraints++;
	}

	constraintId += nbConstraints;
	c1_d.endEdit();
	c2_d.endEdit();
}


template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams* /* PARAMS FIRST */, defaulttype::BaseVector *v, const DataVecCoord &, const DataVecCoord &
		, const DataVecDeriv &, const DataVecDeriv &)
{
	unsigned int nb = violations.size();
	for(unsigned int i=0; i<nb; i++)
		v->set(cid+i, violations[i]);
}


template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{
	unsigned int nb = this->mstate2->getSize();
	for(unsigned int i=0; i<nb; i++)
	{
		if(!projected[i]) continue;

		resTab[offset] = new AdaptiveBeamConstraintResolution(&displacements[i]);
		offset += 3;
	}
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{ 
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
}

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
