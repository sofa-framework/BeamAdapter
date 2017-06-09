/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
 *                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
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

#ifndef SOFA_COMPONENT_ENGINE_STEERABLECATHETER_H
#define SOFA_COMPONENT_ENGINE_STEERABLECATHETER_H

#include "initBeamAdapter.h"
#include <sofa/defaulttype/SolidTypes.h>
#include "WireRestShape.h"

// For the events
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/core/objectmodel/KeyreleasedEvent.h>
#include <sofa/simulation/AnimateBeginEvent.h>

//For angle computation
#define PI 3.14159265359

namespace sofa
{

namespace component
{

namespace engine
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

    typedef WireRestShape<DataTypes> Inherited;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    /*!
     * @brief Default Constructor.
    */
    SteerableCatheter():
        m_activeBending( initData(&m_activeBending,(bool)false, "activeBending","Boolean activating the bending of the steerable catheter") ),
        m_deactiveBending( initData(&m_deactiveBending,(bool)false, "deactiveBending","Boolean deactivating the bending of the steerable catheter") ),
        m_angleMax( initData(&m_angleMax,(Real) 180.0, "angleMax","Maximum angle that the catheter can reach \n (in degree [0-360])") ),
        m_flatAngle( initData(&m_flatAngle,(Real) 1.0, "flatAngle","Angle below which we consider the catheter as flat/n (Can't be zero)") ),
        m_bendingRate( initData(&m_bendingRate,(unsigned int) 10, "bendingRate","Nb of step needed to reach the maximum bending angle /n (the lower, the faster)") )
    {
     this->f_listening.setValue(true);
    }

    /*!
    * @brief Default Destructor.
    */
    ~SteerableCatheter(){}

    /// Initialization function
    void init();

    /// Function handling the event (listening needs to be true)
    void handleEvent(core::objectmodel::Event* ev);


    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T* obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
         return core::objectmodel::BaseObject::create(obj, context, arg);
    }

    virtual std::string getTemplateName() const
    {
         return templateName(this);
    }

    static std::string templateName(const SteerableCatheter<DataTypes>* = NULL)
    {
         return DataTypes::Name();
    }


    /// Boolean for bending
    Data<bool> m_activeBending;
    /// Boolean for unbending
    Data<bool> m_deactiveBending;
    /// Maximum angle that the catheter can reach
    Data<Real> m_angleMax;
    /// Minimum angle considering the catheter as flat
    Data<Real> m_flatAngle;
    /// Rate of bending
    Data<unsigned int> m_bendingRate;

protected:
    Real tipLength;
    unsigned int bendingRate;

    Real currentAngleRadian;
    Real maxAngleRadian;
    Real maxBendingDiameter;
    Real maxUnbendingDiameter;
    Real incrementalAngleRadian;
};


} // namespace engine

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_ENGINE_STEERABLECATHETER_H */
