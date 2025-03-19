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
#include <BeamAdapter/initBeamAdapter.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>

namespace beamadapter
{

extern void registerBeamInterpolation(sofa::core::ObjectFactory* factory);
extern void registerWireBeamInterpolation(sofa::core::ObjectFactory* factory);
extern void registerAdaptiveBeamLengthConstraint(sofa::core::ObjectFactory* factory);
extern void registerAdaptiveBeamSlidingConstraint(sofa::core::ObjectFactory* factory);
extern void registerAdaptiveBeamController(sofa::core::ObjectFactory* factory);
extern void registerBeamAdapterActionController(sofa::core::ObjectFactory* factory);
extern void registerInterventionalRadiologyController(sofa::core::ObjectFactory* factory);
extern void registerSutureController(sofa::core::ObjectFactory* factory);
extern void registerSteerableCatheter(sofa::core::ObjectFactory* factory);
extern void registerWireRestShape(sofa::core::ObjectFactory* factory);
extern void registerAdaptiveBeamForceFieldAndMass(sofa::core::ObjectFactory* factory);
extern void registerAdaptiveInflatableBeamForceField(sofa::core::ObjectFactory* factory);
extern void registerAdaptiveBeamMapping(sofa::core::ObjectFactory* factory);
extern void registerBeamLengthMapping(sofa::core::ObjectFactory* factory);
extern void registerMultiAdaptiveBeamMapping(sofa::core::ObjectFactory* factory);
extern void registerRodMeshSection(sofa::core::ObjectFactory* factory);
extern void registerRodSpireSection(sofa::core::ObjectFactory* factory);
extern void registerRodStraightSection(sofa::core::ObjectFactory* factory);
extern void registerBeamProjectionDifferenceMultiMapping(sofa::core::ObjectFactory* factory);

extern "C" {
    SOFA_BEAMADAPTER_API void initExternalModule();
    SOFA_BEAMADAPTER_API const char* getModuleLicense();
    SOFA_BEAMADAPTER_API const char* getModuleName();
    SOFA_BEAMADAPTER_API const char* getModuleVersion();
    SOFA_BEAMADAPTER_API const char* getModuleDescription();
    SOFA_BEAMADAPTER_API void registerObjects(sofa::core::ObjectFactory* factory);
}

void initBeamAdapter()
{
    static bool first = true;
    if (first)
    {
        // make sure that this plugin is registered into the PluginManager
        sofa::helper::system::PluginManager::getInstance().registerPlugin(MODULE_NAME);
        
        first = false;
    }
}

//Here are just several convenient functions to help user to know what contains the plugin

void initExternalModule()
{
    initBeamAdapter();
}

const char* getModuleLicense()
{
    return "INRIA and Digital-Trainers";
}

const char* getModuleName()
{
    return MODULE_NAME;
}

const char* getModuleVersion()
{
    return MODULE_VERSION;
}

const char* getModuleDescription()
{
    return "A dynamic adapter that modulates the DOF repartition of a beam model according to its radius of curvature.";
}

void registerObjects(sofa::core::ObjectFactory* factory)
{
    registerBeamInterpolation(factory);
    registerWireBeamInterpolation(factory);
    registerAdaptiveBeamLengthConstraint(factory);
    registerAdaptiveBeamSlidingConstraint(factory);
    registerAdaptiveBeamController(factory);
    registerBeamAdapterActionController(factory);
    registerInterventionalRadiologyController(factory);
    registerSutureController(factory);
    registerSteerableCatheter(factory);
    registerWireRestShape(factory);
    registerAdaptiveBeamForceFieldAndMass(factory);
    registerAdaptiveInflatableBeamForceField(factory);
    registerAdaptiveBeamMapping(factory);
    registerBeamLengthMapping(factory);
    registerMultiAdaptiveBeamMapping(factory);
    registerRodMeshSection(factory);
    registerRodSpireSection(factory);
    registerRodStraightSection(factory);
    registerBeamProjectionDifferenceMultiMapping(factory);
}

} // namespace beamadapter
