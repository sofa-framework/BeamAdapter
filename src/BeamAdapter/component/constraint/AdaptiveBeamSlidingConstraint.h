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
#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <BeamAdapter/component/WireBeamInterpolation.h>


namespace sofa::component::constraintset
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
using sofa::linearalgebra::BaseVector ;
using sofa::core::objectmodel::Data ;
using sofa::type::Vec ;
using sofa::component::fem::WireBeamInterpolation ;

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

    virtual void getConstraintResolution(const ConstraintParams* cParams, std::vector<ConstraintResolution*>& resTab, unsigned int& offset) override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////

protected:
    AdaptiveBeamSlidingConstraint(TypedMechanicalState* object1, TypedMechanicalState* object2) ;
    AdaptiveBeamSlidingConstraint(TypedMechanicalState* object) ;
    AdaptiveBeamSlidingConstraint() ;
    virtual ~AdaptiveBeamSlidingConstraint() = default;

private:
    void internalInit();

    bool getCurvAbsOfProjection(const Vec3& x_input, const VecCoord& x, Real& x_output, const Real& tolerance);

    unsigned int          m_cid;
    int                   m_nbConstraints; 		/*!< number of constraints created */
    type::vector<Real>     m_violations;
    type::vector<Real>     m_previousPositions;	/*!< the position on which each point was projected */
    type::vector<double>   m_displacements; 		/*!< displacement=double for compatibility with constraint resolution */
    type::vector<bool>     m_projected;

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



/*! When the Newton method doesn't converge when searching for the closest projection on the curve,
 *  we will start a dichotomic search, for now a unoptimized brute force approach.
 */
template<class DataTypes>
class ProjectionSearch
{
public:
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::Coord Coord;
    typedef typename Coord::value_type Real;

    typedef Vec<3, Real> Vec3;

    /**
     * @brief Default Constructor.
     */
    ProjectionSearch(WireBeamInterpolation<DataTypes>* inter, const Vec3& x_input, const VecCoord& vecX, Real& xcurv_output, const Real& tol)
        : m_interpolation(inter),
        m_x(vecX),
        m_target(x_input),
        m_found(false),
        m_tolerance(tol),
        m_totalIterations(0),
        m_dichotomicIterations(0),
        m_searchDirection(0),
        m_e(xcurv_output)
    {}

    WireBeamInterpolation<DataTypes>* m_interpolation; 			///< The interpolation using this class
    const VecCoord& m_x;											///< The positions of the beams we are working on*/
    Vec3 m_target;												///< The point to be projected on the curve*/
    bool m_found;													///< True when the estimation is acceptable for the given tolerance*/
    Real m_tolerance;												///< Tolerance for the end of the search (projection of the estimation on the tangent is < than tolerance)*/
    unsigned int m_beamIndex;									///< The current beam for the search
    unsigned int m_totalIterations;											///< # of iterations (including beam changes)
    unsigned int m_dichotomicIterations;										///< # of iterations for the current beam
    int m_searchDirection;										///< When we change beam, which direction is it ?
    Real m_e;
    Real m_le;
    Real m_de;												///< Current estimation of the projection, its value in the current beam ([0,1]), and its distance to the target
    Real m_beamStart, m_beamEnd;									///< Abscissa of the current beam extremities
    Real m_segStart, m_segEnd;										///< Abscissa of the current search segment extremities
    Real range, rangeSampling;									    ///< The length of the current search segment, and the distance between 2 sampling points
    Vec3 P0, P1, P2, P3;											    ///< Current beam control points

    //TODO(dmarchal 2017-05-17) Je ne suis pas sur du tout de la robustesse de ce code.
    //Pourquoi s_sampling n'est pas une constexpr ?
    static const unsigned int s_sampling{ 10 };    				    ///< How many points do we consider each step ? (at least > 3 or it will never converge)
    Real m_distTab[s_sampling + 1];									///< Array of distances

    bool doSearch(Real& result);
    void initSearch(Real curvAbs);
    void newtonMethod();
    bool changeCurrentBeam(int index);
    Real computeDistAtCurvAbs(Real curvAbs);
    bool testForProjection(Real curvAbs);
};


#if !defined(SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEBEAMSLIDINGCONSTRAINT_CPP)
extern template class SOFA_BEAMADAPTER_API AdaptiveBeamSlidingConstraint<defaulttype::Rigid3Types>;
#endif

} // namespace _adaptiveBeamSlidingConstraint_

/// 'Export' the objects defined in the private namespace into the 'public' one so that
/// we can use them instead as sofa::component::constraintset::AdaptiveBeamSlidingConstraint.
using _adaptiveBeamSlidingConstraint_::AdaptiveBeamSlidingConstraint ;

} // namespace sofa::component::constraintset
