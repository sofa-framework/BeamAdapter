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
    , d_actions(initData(&d_actions, "actions", "List of actions to script the intervention"))
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
void BeamAdapterActionController<DataTypes>::onBeginAnimationStep(const double /*dt*/)
{
    const type::vector<Real>& times = d_timeSteps.getValue();
    if (!times.empty())
    {
        const auto currentTime = this->getContext()->getTime();
        if (readStep < times.size())
        {
            Real time = times[readStep];
            if (currentTime >= time) // check if another key time has been reached and change action
            {
                currAction = BeamAdapterAction(d_actions.getValue()[readStep]);
                readStep++;
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

template <class DataTypes>
void BeamAdapterActionController<DataTypes>::onMouseEvent(core::objectmodel::MouseEvent* mev)
{

}


} // namespace sofa::component::controller
