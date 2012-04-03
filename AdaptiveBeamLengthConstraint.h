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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMLENGTHCONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMLENGTHCONSTRAINT_H

#include "WireBeamInterpolation.h"
#include <sofa/core/behavior/Constraint.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <iostream>
#include <sofa/defaulttype/Vec.h>
#include <map>

namespace sofa
{

namespace component
{

namespace constraintset
{
using sofa::helper::vector;

template<class DataTypes>
class AdaptiveBeamLengthConstraint : public core::behavior::Constraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamLengthConstraint,DataTypes), SOFA_TEMPLATE(sofa::core::behavior::Constraint,DataTypes));

	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::MatrixDeriv MatrixDeriv;
	typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
	typedef typename DataTypes::Coord Coord;
	typedef typename DataTypes::Deriv Deriv;
	typedef typename Coord::value_type Real;
	typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename core::behavior::Constraint<DataTypes> Inherit;
	typedef core::objectmodel::Data<VecCoord>		DataVecCoord;
	typedef core::objectmodel::Data<VecDeriv>		DataVecDeriv;
	typedef core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;
	typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
	typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;
	typedef typename DataTypes::Coord::Pos Pos;
	typedef typename DataTypes::Coord::Rot Rot;
	typedef sofa::defaulttype::Vec<3, Real> Vec3;
	typedef sofa::defaulttype::Vec<6, Real> Vec6;
    typedef typename  std::map<Real, double>::iterator MapIterator;

protected:
	
	unsigned int cid;
		
	int nbConstraints; // number of constraints created
	std::vector<Real> violations;
 //   std::vector<unsigned int> activated_beams;	// Not used ?
	std::vector<Real> activatedBeamsAbscissa;
	std::map<Real, double> prevForces;	// Map abscissa <-> previous constraint force

	Data<Real> m_alarmLength, m_constrainedLength;

    SingleLink<AdaptiveBeamLengthConstraint<DataTypes>, fem::WireBeamInterpolation<DataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_interpolation;
	
	void internalInit();
	
    AdaptiveBeamLengthConstraint(MechanicalState* object = NULL)
    : Inherit(object)
	, m_alarmLength(initData(&m_alarmLength, (Real)1.02, "alarmLength", "Elongation before creating a constraint"))
	, m_constrainedLength(initData(&m_constrainedLength, (Real)1.05, "constrainedLength", "Allowed elongation of a beam"))
    , m_interpolation(initLink("interpolation", "link to the interpolation component in the scene"))
	{
	}
	
    virtual ~AdaptiveBeamLengthConstraint()
	{
	}

public:
	virtual void reset();
	
	virtual void init();
	
	virtual  void buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST =core::ConstraintParams::defaultInstance()*/, DataMatrixDeriv &c, unsigned int &cIndex, const DataVecCoord &x);

	void getConstraintViolation(const core::ConstraintParams* cParams /* PARAMS FIRST =core::ConstraintParams::defaultInstance()*/, defaulttype::BaseVector *viol, const DataVecCoord &x,const DataVecDeriv &v);

	virtual void getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset);

	void draw(const core::visual::VisualParams* vparams);
};

class AdaptiveBeamLengthConstraintResolution : public core::behavior::ConstraintResolution
{
public:
	AdaptiveBeamLengthConstraintResolution(double* initF=NULL) : _initF(initF) { nbLines = 1; }
	virtual void init(int line, double** /*w*/, double* force);
	virtual void resolution(int line, double** w, double* d, double* force);
	virtual void store(int line, double* force, bool /*convergence*/);
	
protected:
	double* _initF;
};

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMLENGTHCONSTRAINT_H
