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
//
// C++ Implementation : WireBeamInterpolation / AdaptiveBeamForceFieldAndMass
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
#pragma once

#include <BeamAdapter/config.h>

#include <BeamAdapter/component/engine/WireRestShape.h>
#include <BeamAdapter/component/BeamInterpolation.h>

#include <sofa/core/behavior/ForceField.h>
#include <sofa/core/behavior/Mass.h>
#include <sofa/core/objectmodel/Data.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/topology/BaseMeshTopology.h>

#include <sofa/type/vector.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>

#include <sofa/core/objectmodel/BaseObject.h>


namespace sofa::component::fem
{

namespace _wirebeaminterpolation_
{
using sofa::component::fem::BeamInterpolation ;
using sofa::core::topology::BaseMeshTopology ;
using sofa::type::Quat ;
using sofa::type::Vec ;
using sofa::type::Vec3d ;
using sofa::type::vector;
/*!
 * \class WireBeamInterpolation
 * WireAdaptiveBeam Interpolation provides the implementation of a 1D parameteric beam model.
 *
 * Compute Finite Element elastic force and mass based on Adaptive 6D beam elements.
 * - Adaptive beam interpolation
 * - Adaptive Force and Mass computation
 * - Adaptive Mapping
 *
 * TODO(dmarchal 2017-05-17) Please specify who/when this will be done
 * TODO : put in a separate class what is specific to wire shape !
 */
template<class DataTypes>
class WireBeamInterpolation : public virtual BeamInterpolation<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(WireBeamInterpolation, DataTypes) ,
               SOFA_TEMPLATE(BeamInterpolation, DataTypes) );

    typedef BeamInterpolation<DataTypes> Inherited;

    typedef typename Inherited::VecCoord VecCoord;
    typedef typename Inherited::VecDeriv VecDeriv;
    typedef typename Inherited::VecReal VecReal;
    typedef typename Inherited::Coord Coord;
    typedef typename Inherited::Deriv Deriv;

    typedef typename Inherited::Real Real;
    typedef typename Inherited::Index Index;
    typedef typename Inherited::ElementID ElementID;
    typedef typename Inherited::VecElementID VecElementID;
    typedef typename Inherited::VecEdges VecEdges;
    typedef typename Inherited::VecIndex VecIndex;

    typedef typename Inherited::Transform Transform;
    typedef typename Inherited::SpatialVector SpatialVector;

    typedef typename Inherited::Vec2 Vec2;
    typedef typename Inherited::Vec3 Vec3;
    typedef typename Inherited::Vec6 Vec6;

    WireBeamInterpolation(sofa::component::engine::WireRestShape<DataTypes> *_restShape = NULL);

    virtual ~WireBeamInterpolation();

    void init() override;
    void bwdInit() override;
    void reinit() override { init(); bwdInit(); }

    using BeamInterpolation<DataTypes>::addBeam;

    void addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1,
                 const Transform &DOF0_H_Node0, const Transform &DOF1_H_Node1);

    virtual void getSamplingParameters(type::vector<Real>& xP_noticeable, type::vector< int>& nbP_density) override
    {
        this->m_restShape->getSamplingParameters(xP_noticeable, nbP_density);
    }

    virtual Real getRestTotalLength() override
    {
        return this->m_restShape->getLength();
    }

    virtual void getCollisionSampling(Real &dx, const Real& x_localcurv_abs) override
    {
        this->m_restShape->getCollisionSampling(dx,x_localcurv_abs);
    }

    virtual void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines) override
    {
        this->m_restShape->getNumberOfCollisionSegment(dx,numLines);
    }


    virtual void getYoungModulusAtX(int beamId,Real& x_curv, Real& youngModulus, Real& cPoisson) override
    {
        this->getAbsCurvXFromBeam(beamId, x_curv);
        this->m_restShape->getYoungModulusAtX(x_curv, youngModulus, cPoisson);
    }

    virtual void getRestTransform(unsigned int edgeInList, Transform &local0_H_local1_rest);
    virtual void getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest) override;
    virtual void getBeamAtCurvAbs(const Real& x_input, unsigned int &edgeInList_output, Real& baryCoord_output, unsigned int start=0) override;

    void getCurvAbsAtBeam(const unsigned int &edgeInList_input, const Real& baryCoord_input, Real& x_output);
    bool getApproximateCurvAbs(const Vec3& x_input, const VecCoord& x,  Real& x_output);	// Project a point on the segments, return false if cant project
    bool getCurvAbsOfProjection(const Vec3& x_input, const VecCoord& x, Real& x_output, const Real& tolerance);

    bool breaksInTwo(const Real &x_min_out,  Real &x_break, int &numBeamsNotUnderControlled );

    void setPathToRestShape(const std::string &o){m_restShape.setPath(o);}

    void getRestTransformOnX(Transform &global_H_local, const Real &x)
    {
        if(this->m_restShape)
        {
            this->m_restShape->getRestTransformOnX(global_H_local, x);
            return;
        }
        else
        {
            global_H_local.set(Vec3(x,0,0), Quat<Real>());

        }
    }

    SingleLink<WireBeamInterpolation<DataTypes>, sofa::component::engine::WireRestShape<DataTypes>,
    BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_restShape; /*! link on an external rest-shape*/


    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using  BeamInterpolation<DataTypes>::m_componentstate ;
    ////////////////////////////////////////////////////////////////////////////

public:

    template<class T>
    static bool canCreate(T* obj, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        return Inherited::canCreate(obj,context,arg);
    }

    template<class T>
    static typename T::SPtr  create(T* tObj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg) ;
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

    WireBeamInterpolation<DataTypes> *m_interpolation; 			///< The interpolation using this class
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
    Vec3 P0,P1,P2,P3;											    ///< Current beam control points

    //TODO(dmarchal 2017-05-17) Je ne suis pas sur du tout de la robustesse de ce code.
    //Pourquoi s_sampling n'est pas une constexpr ?
    static const unsigned int s_sampling {10} ;    				    ///< How many points do we consider each step ? (at least > 3 or it will never converge)
    Real m_distTab[s_sampling+1];									///< Array of distances

    bool doSearch(Real& result);
    void initSearch(Real curvAbs);
    void newtonMethod();
    bool changeCurrentBeam(int index);
    Real computeDistAtCurvAbs(Real curvAbs);
    bool testForProjection(Real curvAbs);
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_WIREBEAMINTERPOLATION_CPP)
extern template class SOFA_BEAMADAPTER_API WireBeamInterpolation<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace _wirebeaminterpolation_

/// Import the privately defined into the expected sofa namespace.
using _wirebeaminterpolation_::WireBeamInterpolation ;

} // namespace sofa::component::fem
