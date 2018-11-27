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
//
// Description:
//
//
// Author: Hugo Talbot, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_ENGINE_STEERABLECATHETER_INL
#define SOFA_COMPONENT_ENGINE_STEERABLECATHETER_INL

#include "SteerableCatheter.h"
#include "math.h"

namespace sofa
{

namespace component
{

namespace engine
{

template <class DataTypes>
void SteerableCatheter<DataTypes>::init()
{
    Inherited::init();

    bendingRate = m_bendingRate.getValue();
    Real angleMax = m_angleMax.getValue();
    Real flatAngle = m_flatAngle.getValue();

    if(angleMax<=0.0 || angleMax>360.0)
    {
        std::cout<<"(SteerableCatheter) Wrong angleMax given: should be [0-360]. Bending deactivated."<<std::endl;
        return;
    }
    else if(bendingRate <=0)
    {
        std::cout<<"(SteerableCatheter) Wrong bendingRate given: should be non zero."<<std::endl;
        return;
    }
    else
    {
        if( flatAngle <= 0.0)
        {
            std::cout<<"(SteerableCatheter) Wrong flatAngle given: must be higher than 0°."<<std::endl;
            flatAngle = 1.0;
            m_flatAngle.setValue(1.0);
        }

        // Initialize the increment for bending/unbending
        tipLength = this->d_length.getValue() - this->d_straightLength.getValue();
        maxAngleRadian = (angleMax*PI) / 360; // all angle are actually angle/2
        incrementalAngleRadian = maxAngleRadian / bendingRate;

        //Limiting the straight position avoiding going to infinity
        maxUnbendingDiameter = 360.0 * tipLength / (PI*flatAngle);

        // Reajust the initial spireDiameter (associated angle) to be multiple of incrementalAngleRadian
        Real _spireDiameter = this->d_spireDiameter.getValue();

        if(_spireDiameter>maxUnbendingDiameter || _spireDiameter == 0.0)
        {
                if(_spireDiameter>maxUnbendingDiameter)
                    std::cout<<"(SteerableCatheter) Wrong spireDiameter: must be below "<<maxUnbendingDiameter<<std::endl;
                else
                    std::cout<<"(SteerableCatheter) Wrong spireDiameter: must be non-zero (==infinite curvature) "<<std::endl;
                _spireDiameter = maxUnbendingDiameter;
                currentAngleRadian = flatAngle * PI / 360;
                this->d_spireDiameter.setValue( maxUnbendingDiameter );
        }
        else
        {
            if(tipLength != 0.0)
            {
                Real initialAngleRadian  = tipLength/_spireDiameter;
                unsigned int initialAngleIncrement = (unsigned int)(initialAngleRadian/incrementalAngleRadian);
                if(initialAngleIncrement==0)
                    initialAngleIncrement = 1;
                currentAngleRadian = (Real)(initialAngleIncrement)*incrementalAngleRadian;
                this->d_spireDiameter.setValue( tipLength / currentAngleRadian );
            }
            else
            {
                this->d_spireDiameter.setValue( 0.0 );
                currentAngleRadian = 0.0;
                this->f_listening.setValue(false);
            }
        }
    }
}


template <class DataTypes>
void SteerableCatheter<DataTypes>::handleEvent(core::objectmodel::Event* event)
{
    // *****************************
    // Update bending at benginEvent
    if (dynamic_cast<sofa::simulation::AnimateBeginEvent *>(event))
    {
        //Bending activated
        if(m_activeBending.getValue())
        {
            if(currentAngleRadian < maxAngleRadian)
            {
                currentAngleRadian += incrementalAngleRadian;
                this->d_spireDiameter.setValue( tipLength / currentAngleRadian );
            }
        }
        //Bending activated
        if(m_deactiveBending.getValue())
        {
            if(currentAngleRadian > incrementalAngleRadian )
            {
                currentAngleRadian -= incrementalAngleRadian;
                this->d_spireDiameter.setValue( tipLength / currentAngleRadian );
            }
        }
    }




    // **********************
    // Keyboard events ...
    if(sofa::core::objectmodel::KeypressedEvent* ev = dynamic_cast<sofa::core::objectmodel::KeypressedEvent*>(event))
    {
        switch(ev->getKey())
        {
        case '6':
            {
                //********
                // Bending
                m_activeBending.setValue(true);
                break;
            }

        case '9':
            {
                //**********
                // Unbending
                m_deactiveBending.setValue(true);
                break;
            }
        }
    }


    if(sofa::core::objectmodel::KeyreleasedEvent* ev = dynamic_cast<sofa::core::objectmodel::KeyreleasedEvent*>(event))
    {
        switch(ev->getKey())
        {
        case '6':
            {
                //********
                // Bending
                m_activeBending.setValue(false);
                break;
            }

        case '9':
            {
                //**********
                // Unbending
                m_deactiveBending.setValue(false);
                break;
            }
        }
    }
}


} // namespace engine
} // namespace component
} // namespace sofa

#endif /* SOFA_COMPONENT_ENGINE_STEERABLECATHETER_INL */
