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
    : d_writeMode(initData(&d_writeMode, true, "writeMode", "If true, will accumulate actions from keyboard and dump the actions and times when key 'E' is pressed."))
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
        return;
    }

    // the controller must listen to the event (in particular BeginAnimationStep event)
    this->f_listening.setValue(true);

    if (d_writeMode.getValue() && (d_actions.isSet() || d_actionString.isSet()))
    {
        msg_warning() << "WriteMode is set to on but a list of actions has been set as input. The list will be overwritten.";
    }

    interventionCtrl* ctrl = l_interventionController.get();
    ctrl->useBeamAction(true);

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
        m_currAction = BeamAdapterAction::SPIN_RIGHT;
        break;
    case 18: // gauche = 18
        m_currAction = BeamAdapterAction::SPIN_LEFT;
        break;
    case 19: // fleche haut = 19
        m_currAction = BeamAdapterAction::MOVE_FORWARD;
        break;
    case 21: // bas = 21
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
    interventionCtrl* ctrl = l_interventionController.get();

    if (d_writeMode.getValue())
    {
        if (m_currAction == BeamAdapterAction::NO_ACTION)
            return ctrl->applyInterventionalRadiologyController();

        ctrl->applyAction(m_currAction);

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
        m_currAction = BeamAdapterAction::NO_ACTION;
    }
    else
    {
        const type::vector<Real>& times = d_timeSteps.getValue();
        if (!times.empty() && m_readStep < int(times.size()))
        {
            const Real& time = times[m_readStep];

            if (currentTime >= time) // check if another key time has been reached and change action
            {                
                m_currAction = BeamAdapterAction(d_actions.getValue()[m_readStep]);
                m_readStep++;

                ctrl->applyAction(m_currAction);
            }
        }
    }

    ctrl->applyInterventionalRadiologyController();
}


} // namespace sofa::component::controller
