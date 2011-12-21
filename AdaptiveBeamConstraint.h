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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMCONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMCONSTRAINT_H

#include "WireBeamInterpolation.h"
#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <iostream>
#include <sofa/defaulttype/Vec.h>

namespace sofa
{

namespace component
{

namespace constraintset
{
using sofa::helper::vector;

template<class DataTypes>
class AdaptiveBeamConstraint : public core::behavior::PairInteractionConstraint<DataTypes>
{
public:
	SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamConstraint,DataTypes), SOFA_TEMPLATE(sofa::core::behavior::PairInteractionConstraint,DataTypes));

	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::MatrixDeriv MatrixDeriv;
	typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
	typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Deriv Deriv;
	typedef typename Coord::value_type Real;
	typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
	typedef typename core::behavior::PairInteractionConstraint<DataTypes> Inherit;
	typedef core::objectmodel::Data<VecCoord>		DataVecCoord;
	typedef core::objectmodel::Data<VecDeriv>		DataVecDeriv;
	typedef core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;
	typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
	typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;
	typedef typename DataTypes::Coord::Pos Pos;
	typedef typename DataTypes::Coord::Rot Rot;
	typedef sofa::defaulttype::Vec<3, Real> Vec3;
	typedef sofa::defaulttype::Vec<6, Real> Vec6;

protected:
	
	unsigned int cid;
		
	int nbConstraints; // number of constraints created
	std::vector<Real> violations;
        std::vector<Real> previousPositions;// the position on which each point was projected
        std::vector<double> displacements; 	// displacement=double for compatibility with constraint resolution
	std::vector<bool> projected;

	SingleLink<AdaptiveBeamConstraint<DataTypes>, fem::WireBeamInterpolation<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_interpolation;
	
	void internalInit();

	
	AdaptiveBeamConstraint(MechanicalState* object1, MechanicalState* object2)
	: Inherit(object1, object2)
	, m_interpolation(initLink("interpolation", "link to the interpolation component in the scene"))
	{
	}
	
	AdaptiveBeamConstraint(MechanicalState* object)
	: Inherit(object, object)
	, m_interpolation(initLink("interpolation", "link to the interpolation component in the scene"))
	{
	}

	AdaptiveBeamConstraint()
	: m_interpolation(initLink("interpolation", "link to the interpolation component in the scene"))
	{
	}
	
	virtual ~AdaptiveBeamConstraint()
	{
	}
public:
	virtual void reset();
	
	virtual void init();
	
	void buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST =core::ConstraintParams::defaultInstance()*/, DataMatrixDeriv &c1, DataMatrixDeriv &c2, unsigned int &cIndex, const DataVecCoord &x1, const DataVecCoord &x2);

	void getConstraintViolation(const core::ConstraintParams* cParams /* PARAMS FIRST =core::ConstraintParams::defaultInstance()*/, defaulttype::BaseVector *v, const DataVecCoord &x1, const DataVecCoord &x2
		, const DataVecDeriv &v1, const DataVecDeriv &v2);

	virtual void getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset);

	void draw(const core::visual::VisualParams* vparams);
};

class AdaptiveBeamConstraintResolution : public core::behavior::ConstraintResolution
{
public:
	AdaptiveBeamConstraintResolution(double* sliding = NULL)
	: _slidingDisp(sliding) { nbLines = 3; }
	virtual void init(int line, double** w, double* force);
	virtual void resolution(int line, double** w, double* d, double* force);
	virtual void store(int line, double* force, bool /*convergence*/);
	
protected:
	double slidingW;
	double* _slidingDisp;
};

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMCONSTRAINT_H
