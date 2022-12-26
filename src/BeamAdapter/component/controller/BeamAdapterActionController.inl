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
#include <BeamAdapter/component/controller/BeamAdapterActionController.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>

namespace sofa::component::controller
{

using namespace sofa::beamadapter;

template <class DataTypes>
BeamAdapterActionController<DataTypes>::BeamAdapterActionController()
    : d_writeMode(initData(&d_writeMode, false, "writeMode", "If true, will accumulate actions from keyboard and dump the actions and times when key 'E' is pressed."))
    , d_actions(initData(&d_actions, "actions", "List of actions to script the BeamAdapter"))
    , d_actionString(initData(&d_actionString, "actionString", "List of actions as string to script the BeamAdapter"))
    , d_timeSteps(initData(&d_timeSteps, "timeSteps", "List of key times corresponding to the actions"))
    , l_interventionController(initLink("interventionController", "Path to the InterventionalRadiologyController component on scene"))
{
    
}

template <class DataTypes>
void BeamAdapterActionController<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Loading);
    if (!l_interventionController.get())
    {
        msg_error() << "No l_interventionController given. Component will be set to Invalid.";
        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }

    // the controller must listen to the event (in particular BeginAnimationStep event)
    this->f_listening.setValue(true);

    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}


template <class DataTypes>
void BeamAdapterActionController<DataTypes>::onKeyPressedEvent(core::objectmodel::KeypressedEvent* kev)
{
    if (!d_writeMode.getValue())
        return;

    /// Control keys for interventonal Radiology simulations:
    switch (kev->getKey())
    {
    case 'E':
        m_currAction = BeamAdapterAction::NO_ACTION;
        m_exportActions = !m_exportActions;
        break;
    case 'D':
        m_currAction = BeamAdapterAction::DROP_TOOL;
        break;
    case '2':
        m_currAction = BeamAdapterAction::USE_TOOL_2;
        break;
    case '1':
        m_currAction = BeamAdapterAction::USE_TOOL_1;
        break;
    case '0':
        m_currAction = BeamAdapterAction::USE_TOOL_0;
        break;
    case 20: // droite = 20
        if (m_currAction == BeamAdapterAction::SPIN_RIGHT)
            m_currAction = BeamAdapterAction::NO_ACTION;
        else
            m_currAction = BeamAdapterAction::SPIN_RIGHT;

        break;
    case 18: // gauche = 18
        if (m_currAction == BeamAdapterAction::SPIN_LEFT)
            m_currAction = BeamAdapterAction::NO_ACTION;
        else
            m_currAction = BeamAdapterAction::SPIN_LEFT;

        break;
    case 19: // fleche haut = 19
        if (m_currAction == BeamAdapterAction::MOVE_FORWARD)
            m_currAction = BeamAdapterAction::NO_ACTION;
        else
            m_currAction = BeamAdapterAction::MOVE_FORWARD;
        
        break;
    case 21: // bas = 21
        if (m_currAction == BeamAdapterAction::MOVE_BACKWARD)
            m_currAction = BeamAdapterAction::NO_ACTION;
        else
            m_currAction = BeamAdapterAction::MOVE_BACKWARD;

        break;
    default:
        m_currAction = BeamAdapterAction::NO_ACTION;
    break;
    }
}


template <class DataTypes>
void BeamAdapterActionController<DataTypes>::onBeginAnimationStep(const double /*dt*/)
{
    const auto currentTime = this->getContext()->getTime();
    if (d_writeMode.getValue())
    {
        interventionCtrl* ctrl = l_interventionController.get();
        ctrl->applyAction(m_currAction);

        if (m_lastAction != m_currAction)
        {
            auto times = sofa::helper::WriteAccessor(d_timeSteps);
            auto actions = sofa::helper::WriteAccessor(d_actions);
            times.push_back(currentTime);
            actions.push_back(int(m_currAction));

            if (m_exportActions)
            {
                std::cout << "timeSteps='" << times.wref() << "'" << std::endl;
                std::cout << "actions='" << actions.wref() << "'" << std::endl;
            }

            m_lastAction = m_currAction;
        }

        if (m_currAction >= BeamAdapterAction::SWITCH_NEXT_TOOL) // action regarding tool needs only to be triggered once
        {
            m_currAction = BeamAdapterAction::NO_ACTION;
        }
    }
    else
    {
        const type::vector<Real>& times = d_timeSteps.getValue();
        if (!times.empty())
        {            
            if (m_readStep < times.size())
            {
                Real time = times[m_readStep];
                if (currentTime >= time) // check if another key time has been reached and change action
                {
                    m_currAction = BeamAdapterAction(d_actions.getValue()[m_readStep]);
                    m_readStep++;
                }
            }

            interventionCtrl* ctrl = l_interventionController.get();
            ctrl->applyAction(m_currAction);

            if (m_currAction >= BeamAdapterAction::SWITCH_NEXT_TOOL) // action regarding tool needs only to be triggered once
            {
                m_currAction = BeamAdapterAction::NO_ACTION;
            }
        }
    }
}


} // namespace sofa::component::controller
