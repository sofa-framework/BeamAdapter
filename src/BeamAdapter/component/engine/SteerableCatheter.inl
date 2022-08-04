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
//
// Description:
//
//
// Author: Hugo Talbot, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#pragma once

#include <BeamAdapter/component/engine/SteerableCatheter.h>

#include <math.h>

namespace sofa::component::engine
{

template <class DataTypes>
SteerableCatheter<DataTypes>::SteerableCatheter():
    d_activeBending( initData(&d_activeBending,(bool)false, "activeBending","Boolean activating the bending of the steerable catheter") ),
    d_deactiveBending( initData(&d_deactiveBending,(bool)false, "deactiveBending","Boolean deactivating the bending of the steerable catheter") ),
    d_angleMax( initData(&d_angleMax,(Real) 180.0, "angleMax","Maximum angle that the catheter can reach \n (in degree [0-360])") ),
    d_flatAngle( initData(&d_flatAngle,(Real) 1.0, "flatAngle","Angle below which we consider the catheter as flat/n (Can't be zero)") ),
    d_bendingRate( initData(&d_bendingRate,(unsigned int) 10, "bendingRate","Nb of step needed to reach the maximum bending angle /n (the lower, the faster)") )
{
    f_listening.setValue(true);
}

template <class DataTypes>
void SteerableCatheter<DataTypes>::init()
{
    Inherit1::init();

    Real flatAngle = d_flatAngle.getValue();
    Real angleMax = d_angleMax.getValue();

    if(angleMax<=0.0 || angleMax>360.0)
    {
        msg_info()<<"(SteerableCatheter) Wrong angleMax given: should be [0-360]. Bending deactivated.";
        return;
    }
    else if(d_bendingRate.getValue() <=0)
    {
        msg_info()<<"(SteerableCatheter) Wrong bendingRate given: should be non zero.";
        return;
    }
    else
    {
        if(flatAngle <= 0.0)
        {
            msg_info()<<"(SteerableCatheter) Wrong flatAngle given: must be higher than 0°.";
            flatAngle = 1.0;
            d_flatAngle.setValue(1.0);
        }

        // Initialize the increment for bending/unbending
        m_tipLength = d_length.getValue() - d_straightLength.getValue();
        m_maxAngleRadian = (angleMax*M_PI) / 360; // all angle are actually angle/2
        m_incrementalAngleRadian = m_maxAngleRadian / d_bendingRate.getValue();

        //Limiting the straight position avoiding going to infinity
        m_maxUnbendingDiameter = 360.0 * m_tipLength / (M_PI*flatAngle);

        // Reajust the initial spireDiameter (associated angle) to be multiple of incrementalAngleRadian
        Real spireDiameter = d_spireDiameter.getValue();

        if(spireDiameter>m_maxUnbendingDiameter || spireDiameter == 0.0)
        {
                if(spireDiameter>m_maxUnbendingDiameter)
                {
                    msg_info()<<"(SteerableCatheter) Wrong spireDiameter: must be below "<<m_maxUnbendingDiameter;
                }
                else
                {
                    msg_info()<<"(SteerableCatheter) Wrong spireDiameter: must be non-zero (==infinite curvature) ";
                }
                m_currentAngleRadian = flatAngle * M_PI / 360;
                d_spireDiameter.setValue( m_maxUnbendingDiameter );
        }
        else
        {
            if(m_tipLength != 0.0)
            {
                Real initialAngleRadian  = m_tipLength/spireDiameter;
                unsigned int initialAngleIncrement = (unsigned int)(initialAngleRadian/m_incrementalAngleRadian);
                if(initialAngleIncrement==0)
                    initialAngleIncrement = 1;
                m_currentAngleRadian = (Real)(initialAngleIncrement)*m_incrementalAngleRadian;
                d_spireDiameter.setValue( m_tipLength / m_currentAngleRadian );
            }
            else
            {
                d_spireDiameter.setValue( 0.0 );
                m_currentAngleRadian = 0.0;
                f_listening.setValue(false);
            }
        }
    }
}


template <class DataTypes>
void SteerableCatheter<DataTypes>::handleEvent(core::objectmodel::Event* event)
{
    // *****************************
    // Update bending at benginEvent
    if (simulation::AnimateBeginEvent::checkEventType(event))
    {
        //Bending activated
        if(d_activeBending.getValue())
        {
            if(m_currentAngleRadian < m_maxAngleRadian)
            {
                m_currentAngleRadian += m_incrementalAngleRadian;
                d_spireDiameter.setValue( m_tipLength / m_currentAngleRadian );
            }
        }
        //Bending activated
        if(d_deactiveBending.getValue())
        {
            if(m_currentAngleRadian > m_incrementalAngleRadian )
            {
                m_currentAngleRadian -= m_incrementalAngleRadian;
                d_spireDiameter.setValue( m_tipLength / m_currentAngleRadian );
            }
        }
    }

    // **********************
    // Keyboard events ...
    if (sofa::core::objectmodel::KeypressedEvent::checkEventType(event))
    {
        auto* ev = dynamic_cast<sofa::core::objectmodel::KeypressedEvent*>(event);
        assert(ev != nullptr);

        switch(ev->getKey())
        {
        case '6':
            {
                //********
                // Bending
                d_activeBending.setValue(true);
                break;
            }

        case '9':
            {
                //**********
                // Unbending
                d_deactiveBending.setValue(true);
                break;
            }
        }
    }

    if (sofa::core::objectmodel::KeyreleasedEvent::checkEventType(event))
    {
        auto* ev = dynamic_cast<sofa::core::objectmodel::KeyreleasedEvent*>(event);
        assert(ev != nullptr);

        switch(ev->getKey())
        {
        case '6':
            {
                //********
                // Bending
                d_activeBending.setValue(false);
                break;
            }

        case '9':
            {
                //**********
                // Unbending
                d_deactiveBending.setValue(false);
                break;
            }
        }
    }
}


} // namespace sofa::component::engine
