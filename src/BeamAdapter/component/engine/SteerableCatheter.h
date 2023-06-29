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
// C++ Implementation : AdaptiveBeamController
//
// Description:
//
//
// Author: Hugo Talbot, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
//
#pragma once

#include <BeamAdapter/config.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <BeamAdapter/component/engine/WireRestShape.h>

// For the events
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>
#include <sofa/simulation/AnimateBeginEvent.h>

namespace sofa::component::engine
{

/*!
 * \class SteerableCatheter
 * \brief Component derivating from WireRestShape but having a steerable property
 */

template <class DataTypes>
class SteerableCatheter : public WireRestShape<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(SteerableCatheter,DataTypes),SOFA_TEMPLATE(WireRestShape,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    /*!
     * @brief Default Constructor.
    */
    SteerableCatheter();

    /*!
    * @brief Default Destructor.
    */
    virtual ~SteerableCatheter() = default;

    /// Initialization function
    void init() override;

    /// Function handling the event (listening needs to be true)
    void handleEvent(core::objectmodel::Event* ev) override;

    /// Boolean for bending
    Data<bool> d_activeBending;
    /// Boolean for unbending
    Data<bool> d_deactiveBending;
    /// Maximum angle that the catheter can reach
    Data<Real> d_angleMax;
    /// Minimum angle considering the catheter as flat
    Data<Real> d_flatAngle;
    /// Rate of bending
    Data<unsigned int> d_bendingRate;

    using Inherit1::f_listening;

protected:
    Real m_tipLength;
    Real m_currentAngleRadian;
    Real m_maxAngleRadian;
    Real m_maxUnbendingDiameter;
    Real m_incrementalAngleRadian;

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    Data<Real> d_length;
    Data<Real> d_straightLength;
    Data<Real> d_spireDiameter;
    Data<Real> d_spireHeight;
    ///////////////////////////////////////////////////////////////////////////

};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_STEERABLECATHETER_CPP)
extern template class SOFA_BEAMADAPTER_API SteerableCatheter<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace sofa::component::engine
