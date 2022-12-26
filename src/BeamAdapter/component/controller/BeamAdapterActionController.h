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

template<class DataTypes>
class BeamAdapterActionController : public MechanicalStateController<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BeamAdapterActionController, DataTypes), SOFA_TEMPLATE(MechanicalStateController, DataTypes));

    using interventionCtrl = InterventionalRadiologyController<DataTypes>;
    using Real = typename DataTypes::Real;

    BeamAdapterActionController();
    virtual ~BeamAdapterActionController() = default;

    void init() override;
    void onBeginAnimationStep(const double dt) override;
    void onKeyPressedEvent(core::objectmodel::KeypressedEvent* kev) override;
    void onMouseEvent(core::objectmodel::MouseEvent*) override;

    Data <bool> d_writeMode;
    Data <type::vector<int> > d_actions;
    Data <type::vector<std::string> > d_actionString;
    Data <type::vector<Real> > d_timeSteps;

    SingleLink<BeamAdapterActionController<DataTypes>, InterventionalRadiologyController<DataTypes>, BaseLink::FLAG_STOREPATH | BaseLink::FLAG_STRONGLINK> l_interventionController;

private:
    int m_readStep = 0;
    sofa::beamadapter::BeamAdapterAction currAction = BeamAdapterAction::NO_ACTION;
    sofa::beamadapter::BeamAdapterAction lastAction = BeamAdapterAction::NO_ACTION;

    bool m_exportActions = false;
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_ACTIONCONTROLLER_CPP)
extern template class SOFA_BEAMADAPTER_API BeamAdapterActionController<sofa::defaulttype::Rigid3Types>;
#endif

} /// namespace sofa::component::controller
