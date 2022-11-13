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
    : l_interventionController(initLink("interventionController", "Path to the Interpolation component on scene"))
    , d_writeMode(initData(&d_writeMode, false, "writeMode", "List of actions to script the intervention"))
    , d_actions(initData(&d_actions, "actions", "List of actions to script the intervention"))
    , d_actionString(initData(&d_actionString, "actionString", "List of actions as string to script the intervention"))
    , d_timeSteps(initData(&d_timeSteps, "timeSteps", "List of time to change the action"))
{
    
}

template <class DataTypes>
void BeamAdapterActionController<DataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Loading);
    if (!l_interventionController.get())
    {
        msg_error() << "No l_interventionController given";
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
        currAction = BeamAdapterAction::NO_ACTION;
        m_exportActions = !m_exportActions;
        break;
    case 'D':
        currAction = BeamAdapterAction::DROP_TOOL;
        break;
    case '2':
        currAction = BeamAdapterAction::USE_TOOL_2;
        break;
    case '1':
        currAction = BeamAdapterAction::USE_TOOL_1;
        break;
    case '0':
        currAction = BeamAdapterAction::USE_TOOL_0;
        break;
    case 20: // droite = 20
        if (currAction == BeamAdapterAction::SPIN_RIGHT)
            currAction = BeamAdapterAction::NO_ACTION;
        else
            currAction = BeamAdapterAction::SPIN_RIGHT;

        break;
    case 18: // gauche = 18
        if (currAction == BeamAdapterAction::SPIN_LEFT)
            currAction = BeamAdapterAction::NO_ACTION;
        else
            currAction = BeamAdapterAction::SPIN_LEFT;

        break;
    case 19: // fleche haut = 19
        if (currAction == BeamAdapterAction::MOVE_FORWARD)
            currAction = BeamAdapterAction::NO_ACTION;
        else
            currAction = BeamAdapterAction::MOVE_FORWARD;
        
        break;
    case 21: // bas = 21
        if (currAction == BeamAdapterAction::MOVE_BACKWARD)
            currAction = BeamAdapterAction::NO_ACTION;
        else
            currAction = BeamAdapterAction::MOVE_BACKWARD;

        break;
    default:
        currAction = BeamAdapterAction::NO_ACTION;
    break;
    }
}


template <class DataTypes>
void BeamAdapterActionController<DataTypes>::onBeginAnimationStep(const double /*dt*/)
{
    if (d_writeMode.getValue())
    {
        interventionCtrl* ctrl = l_interventionController.get();
        ctrl->applyAction(currAction);

        if (lastAction != currAction)
        {
            auto times = sofa::helper::WriteAccessor(d_timeSteps);
            auto actions = sofa::helper::WriteAccessor(d_actions);
            times.push_back(getContext()->getTime());
            actions.push_back(int(currAction));

            if (m_exportActions)
            {
                std::cout << "timeSteps='" << times.wref() << "'" << std::endl;
                std::cout << "actions='" << actions.wref() << "'" << std::endl;
            }

            lastAction = currAction;
        }

        if (currAction >= BeamAdapterAction::SWITCH_NEXT_TOOL) // action regarding tool needs only to be triggered once
        {
            currAction = BeamAdapterAction::NO_ACTION;
        }
    }
    else
    {
        const type::vector<Real>& times = d_timeSteps.getValue();
        if (!times.empty())
        {
            const auto currentTime = this->getContext()->getTime();
            if (m_readStep < times.size())
            {
                Real time = times[m_readStep];
                if (currentTime >= time) // check if another key time has been reached and change action
                {
                    currAction = BeamAdapterAction(d_actions.getValue()[m_readStep]);
                    m_readStep++;
                }
            }

            interventionCtrl* ctrl = l_interventionController.get();
            ctrl->applyAction(currAction);

            if (currAction >= BeamAdapterAction::SWITCH_NEXT_TOOL) // action regarding tool needs only to be triggered once
            {
                currAction = BeamAdapterAction::NO_ACTION;
            }
        }
    }
}

template <class DataTypes>
void BeamAdapterActionController<DataTypes>::onMouseEvent(core::objectmodel::MouseEvent* mev)
{

}


} // namespace sofa::component::controller
