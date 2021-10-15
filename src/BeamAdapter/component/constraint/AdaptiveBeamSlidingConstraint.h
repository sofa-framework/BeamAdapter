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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMSLIDINGCONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMSLIDINGCONSTRAINT_H

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <BeamAdapter/component/WireBeamInterpolation.h>


namespace sofa
{

namespace component
{

namespace constraintset
{

/////////////////////////////////// private namespace pattern //////////////////////////////////////
/// To avoid the lacking of names imported with with 'using' in the other's component namespace
/// you should use a private namespace and "export" only this one in the public namespace.
/// This is done at the end of this file, have a look if you are not used to this pattern.
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace _adaptiveBeamSlidingConstraint_
{

using sofa::core::behavior::PairInteractionConstraint ;
using sofa::core::behavior::ConstraintResolution ;
using sofa::core::visual::VisualParams ;
using sofa::core::ConstraintParams ;
using sofa::defaulttype::SolidTypes ;
using sofa::defaulttype::BaseVector ;
using sofa::core::objectmodel::Data ;
using sofa::type::Vec ;
using sofa::component::fem::WireBeamInterpolation ;
using std::vector;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
 * \class AdaptiveBeamSlidingConstraint
 * @brief AdaptiveBeamSlidingConstraint Class constrain a rigid to be attached to a beam (only in position, not the orientation)
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class SOFA_BEAMADAPTER_API AdaptiveBeamSlidingConstraint : public PairInteractionConstraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamSlidingConstraint, DataTypes),
               SOFA_TEMPLATE(PairInteractionConstraint, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename DataTypes::MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename Coord::value_type Real;
    typedef typename core::behavior::MechanicalState<DataTypes> TypedMechanicalState;
    typedef typename core::behavior::PairInteractionConstraint<DataTypes> Inherit;
    typedef Data<VecCoord>		DataVecCoord;
    typedef Data<VecDeriv>		DataVecDeriv;
    typedef Data<MatrixDeriv>    DataMatrixDeriv;
    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;
    typedef typename DataTypes::Coord::Pos Pos;
    typedef typename DataTypes::Coord::Rot Rot;
    typedef Vec<3, Real> Vec3;
    typedef Vec<6, Real> Vec6;

    /////////////////////////// Inherited from BaseObject //////////////////////////////////////////
    virtual void init() override ;
    virtual void reset() override ;
    virtual void draw(const VisualParams* vparams) override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////


    /////////////////////////// Inherited from BaseConstraint //////////////////////////////////////
    virtual void buildConstraintMatrix(const ConstraintParams* cParams,
                                       DataMatrixDeriv &c1, DataMatrixDeriv &c2,
                                       unsigned int &cIndex,
                                       const DataVecCoord &x1, const DataVecCoord &x2) override ;

    virtual void getConstraintViolation(const ConstraintParams* cParams ,
                                        BaseVector *v,
                                        const DataVecCoord &x1, const DataVecCoord &x2,
                                        const DataVecDeriv &v1, const DataVecDeriv &v2) override;

    virtual void getConstraintResolution(const ConstraintParams* cParams, vector<ConstraintResolution*>& resTab, unsigned int& offset) override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////

protected:
    AdaptiveBeamSlidingConstraint(TypedMechanicalState* object1, TypedMechanicalState* object2) ;
    AdaptiveBeamSlidingConstraint(TypedMechanicalState* object) ;
    AdaptiveBeamSlidingConstraint() ;
    virtual ~AdaptiveBeamSlidingConstraint(){}

private:
    void internalInit();

    unsigned int          m_cid;
    int                   m_nbConstraints; 		/*!< number of constraints created */
    vector<Real>     m_violations;
    vector<Real>     m_previousPositions;	/*!< the position on which each point was projected */
    vector<double>   m_displacements; 		/*!< displacement=double for compatibility with constraint resolution */
    vector<bool>     m_projected;

    /*! link to the interpolation component in the scene */
    SingleLink<AdaptiveBeamSlidingConstraint<DataTypes>,
               WireBeamInterpolation<DataTypes>,
               BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_interpolation;


    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using PairInteractionConstraint<DataTypes>::mstate1;
    using PairInteractionConstraint<DataTypes>::mstate2;
    ////////////////////////////////////////////////////////////////////////////
};

#if !defined(SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMSLIDINGCONSTRAINT_CPP)
extern template class SOFA_BEAMADAPTER_API AdaptiveBeamSlidingConstraint<defaulttype::Rigid3Types>;
#endif

} // namespace _adaptiveBeamSlidingConstraint_

/// 'Export' the objects defined in the private namespace into the 'public' one so that
/// we can use them instead as sofa::component::constraintset::AdaptiveBeamSlidingConstraint.
using _adaptiveBeamSlidingConstraint_::AdaptiveBeamSlidingConstraint ;

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINTSET_AdaptiveBeamSlidingConstraint_H
