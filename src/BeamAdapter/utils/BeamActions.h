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

namespace sofa::beamadapter
{
    
    /// \brief Enum class listing all possible actions during Beam navigation
    enum class BeamAdapterAction
    {
        NO_ACTION = 0,
        MOVE_FORWARD,
        MOVE_BACKWARD,
        SPIN_RIGHT,
        SPIN_LEFT,
        SWITCH_NEXT_TOOL,
        SWITCH_PREVIOUS_TOOL,
        DROP_TOOL,
        USE_TOOL_0,
        USE_TOOL_1,
        USE_TOOL_2,
    };

    /// \brief map of action as string keyword instead of int for better clarity in scene scripting
    const static std::map<std::string, BeamAdapterAction> beamActionNames = {
        {"stop", BeamAdapterAction::NO_ACTION},
        {"fwd", BeamAdapterAction::MOVE_FORWARD},
        {"bwd", BeamAdapterAction::MOVE_BACKWARD},
        {"right", BeamAdapterAction::SPIN_RIGHT},
        {"nextT", BeamAdapterAction::SWITCH_NEXT_TOOL},
        {"prevT", BeamAdapterAction::SWITCH_PREVIOUS_TOOL},
        {"dropT", BeamAdapterAction::DROP_TOOL},
        {"tool0", BeamAdapterAction::USE_TOOL_0},
        {"tool1", BeamAdapterAction::USE_TOOL_1},
        {"tool2", BeamAdapterAction::USE_TOOL_2}
    };

    /// static method to convert an action as string into enum class using @sa beamActionNames
    static BeamAdapterAction convertBeamAdapterAction(const std::string& sAction)
    {
        auto bAction = beamActionNames.find(sAction);

        if (bAction != beamActionNames.end())
            return bAction->second;
        else
            return BeamAdapterAction::NO_ACTION;
    }

} // namespace sofa::beamAdapter
