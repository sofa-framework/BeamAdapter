/******************************************************************************
*                              BeamAdapter plugin                             *
*                  (c) 2006 Inria, University of Lille, CNRS                  *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: see Authors.md                                                     *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#pragma once

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/helper/map.h>

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/behavior/Constraint.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/type/Vec.h>

#include <BeamAdapter/component/WireBeamInterpolation.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Forward declarations, see https://en.wikipedia.org/wiki/Forward_declaration
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa::component::constraintset
{

/////////////////////////////////// private namespace pattern //////////////////////////////////////
/// To avoid the lacking of names imported with with 'using' in the other's component namespace
/// you should use a private namespace and "export" only this one in the public namespace.
/// This is done at the end of this file, have a look if you are not used to this pattern.
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace _adaptivebeamlengthconstraint_
{
using sofa::core::behavior::ConstraintResolution ;
using sofa::core::behavior::Constraint ;
using sofa::core::behavior::MechanicalState ;
using sofa::core::ConstraintParams ;
using sofa::core::objectmodel::Data ;
using sofa::type::Vec ;
using sofa::type::vector;
using sofa::defaulttype::SolidTypes;
using sofa::linearalgebra::BaseVector ;


template<typename Real>
class IntervalDefinition
{
public:
    typedef typename SolidTypes<Real>::Transform Transform;
    typedef Vec<3, Real> Vec3;

    /// definition of an interval which length is "constrained"
    /// positions begin /  end of the interval
    Vec3 posBegin, posEnd;

    /// positions free : begin / end of the interval
    Vec3 posFreeBegin, posFreeEnd;

    /// index of the dofs: begin / end of the interval
    unsigned int IdxBegin, IdxEnd;

    /// transform from dof to begin/end of the interval
    Transform dof_H_begin, dof_H_end;

    /// rest length of the interval (if no stretching)
    Real rest_length;

    /// is it in elongation
    bool active;
} ;

/*!
 * \class AdaptiveBeamLengthConstraint
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class AdaptiveBeamLengthConstraint : public Constraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamLengthConstraint,DataTypes),
               SOFA_TEMPLATE(Constraint,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;
    typedef typename DataTypes::Coord::Pos Pos;
    typedef typename DataTypes::Coord::Rot Rot;
    typedef typename std::map<Real, double>::iterator MapIterator;

    typedef MechanicalState<DataTypes> TypedMechanicalState;
    typedef Constraint<DataTypes> Inherit;
    typedef Data<VecCoord>	 	  DataVecCoord;
    typedef Data<VecDeriv> 		  DataVecDeriv;
    typedef Data<MatrixDeriv>     DataMatrixDeriv;
    typedef Vec<3, Real> Vec3;
    typedef Vec<6, Real> Vec6;

public:
    virtual void init() override ;
    virtual void reset() override ;
    virtual void draw(const core::visual::VisualParams* vparams) override ;


    virtual  void buildConstraintMatrix(const ConstraintParams* cParams,
                                        DataMatrixDeriv &c, unsigned int &cIndex, const DataVecCoord &x) override ;

    virtual void getConstraintViolation(const ConstraintParams* cParams, BaseVector *viol,
                                        const DataVecCoord &x,const DataVecDeriv &v) override ;

    virtual void getConstraintResolution(const ConstraintParams* cParams, std::vector<ConstraintResolution*>& resTab, unsigned int& offset) override ;


protected:
    AdaptiveBeamLengthConstraint(TypedMechanicalState* object = nullptr) ;
    virtual ~AdaptiveBeamLengthConstraint() = default;

    void internalInit();

    unsigned int           m_cid {0} ;
    int                    m_nbConstraints {0};
    type::vector<Real>           m_violations;
    std::map<Real, double> m_prevForces;	/// Map abscissa <-> previous constraint force

    Data<Real>             m_alarmLength ;
    Data<Real>             m_constrainedLength ;
    Data<Real>             m_maxBendingAngle ;
    SingleLink<AdaptiveBeamLengthConstraint<DataTypes>,
               fem::WireBeamInterpolation<DataTypes>,
               BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_interpolation;

private:
    void detectElongation(const VecCoord &x, const VecCoord& xfree);
    type::vector<IntervalDefinition<Real>>                                   m_constraintIntervals;
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEBEAMLENGTHCONSTRAINT_CPP)
extern template class SOFA_BEAMADAPTER_API AdaptiveBeamLengthConstraint<defaulttype::Rigid3Types>;
#endif

} /// namespace _adaptivebeamlengthconstraint_

using _adaptivebeamlengthconstraint_::AdaptiveBeamLengthConstraint ;

} /// namespace sofa::component::constraintset
