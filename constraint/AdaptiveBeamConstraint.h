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

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/behavior/PairInteractionConstraint.h>
#include "WireBeamInterpolation.h"


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
namespace _adaptivebeamconstraint_
{

using sofa::core::behavior::PairInteractionConstraint ;
using sofa::core::behavior::ConstraintResolution ;
using sofa::core::visual::VisualParams ;
using sofa::core::ConstraintParams ;
using sofa::defaulttype::SolidTypes ;
using sofa::defaulttype::BaseVector ;
using sofa::core::objectmodel::Data ;
using sofa::defaulttype::Vec ;
using sofa::component::fem::WireBeamInterpolation ;

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////
/*!
 * \class AdaptiveBeamConstraint
 * @brief AdaptiveBeamConstraint Class constrain a rigid to be attached to a beam (only in position, not the orientation)
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class AdaptiveBeamConstraint : public PairInteractionConstraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamConstraint, DataTypes),
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

    virtual void getConstraintResolution(std::vector<ConstraintResolution*>& resTab, unsigned int& offset) override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////

protected:
    AdaptiveBeamConstraint(TypedMechanicalState* object1, TypedMechanicalState* object2) ;
    AdaptiveBeamConstraint(TypedMechanicalState* object) ;
    AdaptiveBeamConstraint() ;
    virtual ~AdaptiveBeamConstraint(){}

private:
    void internalInit();

    unsigned int          m_cid;
    int                   m_nbConstraints; 		/*!< number of constraints created */
    std::vector<Real>     m_violations;
    std::vector<Real>     m_previousPositions;	/*!< the position on which each point was projected */
    std::vector<double>   m_displacements; 		/*!< displacement=double for compatibility with constraint resolution */
    std::vector<bool>     m_projected;

    /*! link to the interpolation component in the scene */
    SingleLink<AdaptiveBeamConstraint<DataTypes>,
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

} // namespace _adaptivebeamconstraint_

/// 'Export' the objects defined in the private namespace into the 'public' one so that
/// we can use them instead as sofa::component::constraintset::AdaptiveBeamConstraint.
using _adaptivebeamconstraint_::AdaptiveBeamConstraint ;

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMCONSTRAINT_H
