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
    typedef typename Inherited::Coord Coord;
    typedef typename Inherited::Deriv Deriv;

    typedef typename Inherited::Real Real;

    typedef typename Inherited::Transform Transform;
    typedef typename Inherited::SpatialVector SpatialVector;

    typedef typename Inherited::Vec2 Vec2;
    typedef typename Inherited::Vec3 Vec3;
    typedef typename Inherited::Quat Quat;

    WireBeamInterpolation(sofa::component::engine::WireRestShape<DataTypes> *_restShape = nullptr);

    virtual ~WireBeamInterpolation() = default;

    void init() override;
    void bwdInit() override;
    void reinit() override { init(); bwdInit(); }

    using BeamInterpolation<DataTypes>::addBeam;

    void addBeam(const BaseMeshTopology::EdgeID &eID  , const Real &length, const Real &x0, const Real &x1,
                 const Transform &DOF0_H_Node0, const Transform &DOF1_H_Node1);

    void getSamplingParameters(type::vector<Real>& xP_noticeable, type::vector< int>& nbP_density) override
    {
        this->m_restShape->getSamplingParameters(xP_noticeable, nbP_density);
    }

    Real getRestTotalLength() override
    {
        return this->m_restShape->getLength();
    }

    void getCollisionSampling(Real &dx, const Real& x_localcurv_abs) override
    {
        this->m_restShape->getCollisionSampling(dx,x_localcurv_abs);
    }

    void getNumberOfCollisionSegment(Real &dx, unsigned int &numLines) override
    {
        this->m_restShape->getNumberOfCollisionSegment(dx,numLines);
    }


    void getYoungModulusAtX(int beamId,Real& x_curv, Real& youngModulus, Real& cPoisson) override
    {
        this->getAbsCurvXFromBeam(beamId, x_curv);
        this->m_restShape->getYoungModulusAtX(x_curv, youngModulus, cPoisson);
    }

    virtual void getRestTransform(unsigned int edgeInList, Transform &local0_H_local1_rest);
    void getSplineRestTransform(unsigned int edgeInList, Transform &local_H_local0_rest, Transform &local_H_local1_rest) override;

    void getCurvAbsAtBeam(const unsigned int &edgeInList_input, const Real& baryCoord_input, Real& x_output);
    bool getApproximateCurvAbs(const Vec3& x_input, const VecCoord& x,  Real& x_output);	// Project a point on the segments, return false if cant project

    

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
            global_H_local.set(Vec3(x,0,0), Quat());
        }
    }

    SingleLink<WireBeamInterpolation<DataTypes>, sofa::component::engine::WireRestShape<DataTypes>,
    BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_restShape; /*! link on an external rest-shape*/


    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using  BeamInterpolation<DataTypes>::d_componentState ;
    ////////////////////////////////////////////////////////////////////////////

public:

    template<class T>
    static bool canCreate(T* obj, sofa::core::objectmodel::BaseContext* context, sofa::core::objectmodel::BaseObjectDescription* arg)
    {
        return Inherited::canCreate(obj,context,arg);
    }

    template<class T>
    static typename T::SPtr  create(T* tObj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg) ;

    /////////////////////////// Deprecated Methods  ////////////////////////////////////////// 
    /// For coils: a part of the coil instrument can be brokenIn2  (by default the point of release is the end of the straight length)
    [[deprecated("Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06")]]
    bool breaksInTwo(const Real& x_min_out, Real& x_break, int& numBeamsNotUnderControlled) {
        SOFA_UNUSED(x_min_out);
        SOFA_UNUSED(x_break);
        SOFA_UNUSED(numBeamsNotUnderControlled);
        msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
        return 0.0;
    }
};


#if !defined(SOFA_PLUGIN_BEAMADAPTER_WIREBEAMINTERPOLATION_CPP)
extern template class SOFA_BEAMADAPTER_API WireBeamInterpolation<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace _wirebeaminterpolation_

/// Import the privately defined into the expected sofa namespace.
using _wirebeaminterpolation_::WireBeamInterpolation ;

} // namespace sofa::component::fem
