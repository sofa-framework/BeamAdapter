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

#include <sofa/component/controller/MechanicalStateController.h>
#include <BeamAdapter/utils/BeamActions.h>
#include <BeamAdapter/component/controller/InterventionalRadiologyController.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>

namespace sofa::component::controller
{

/***
* This class is a SOFA component inheriting from MechanicalStateController.
* It can be used to script the InterventionalRadiologyController using @sa BeamActions 
* If @sa d_writeMode mode is on, each keyboard interaction used to control the Beam will be exported with key times.
* Otherwise, it will load a list of @sa BeamActions with their corresponding key times to replay a navigation. 
*/
template<class DataTypes>
class BeamAdapterActionController : public MechanicalStateController<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BeamAdapterActionController, DataTypes), SOFA_TEMPLATE(MechanicalStateController, DataTypes));

    using BeamAdapterAction = sofa::beamadapter::BeamAdapterAction;
    using interventionCtrl = InterventionalRadiologyController<DataTypes>;
    using Real = typename DataTypes::Real;

    BeamAdapterActionController();
    ~BeamAdapterActionController() override = default;

    /// Method to init the component. Will search for a InterventionalRadiologyController using link @sa l_interventionController
    void init() override;

    /// Method called at each timestep. Will check if an action has to be read or write
    void onBeginAnimationStep(const double dt) override;

    /// Method to control the Beam using keyboard and save the actions in @sa d_actions
    void onKeyPressedEvent(core::objectmodel::KeypressedEvent* kev) override;

    /// Unused metho for mouse event
    void onMouseEvent(core::objectmodel::MouseEvent* ev) override { SOFA_UNUSED(ev);}

    Data <bool> d_writeMode; ///< If true, will accumulate actions in @sa d_actions for export. Press key E for export.
    Data <type::vector<int> > d_actions; ///< List of actions to import or export.
    Data <type::vector<std::string> > d_actionString; ///< List of actions to import or export as string.
    Data <type::vector<Real> > d_timeSteps; ///< List of key times corresponding to BeamActions in @sa d_actions or @sa d_actionString

    /// Link to the InterventionalRadiologyController, controlling the Beam, to script.
    SingleLink<BeamAdapterActionController<DataTypes>, InterventionalRadiologyController<DataTypes>, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_interventionController;

private:
    int m_readStep = 0; ///< counter to the current action to read in @sa d_actions
    BeamAdapterAction m_currAction = BeamAdapterAction::NO_ACTION; ///< Current action imported or to export
    BeamAdapterAction m_lastAction = BeamAdapterAction::NO_ACTION; ///< Previous action imported or to export

    bool m_exportActions = false; ///< Bool to dump actions, will be set to true if key 'E' is pressed
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_ACTIONCONTROLLER_CPP)
extern template class SOFA_BEAMADAPTER_API BeamAdapterActionController<sofa::defaulttype::Rigid3Types>;
#endif

} /// namespace sofa::component::controller
