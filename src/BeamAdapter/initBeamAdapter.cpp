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
using sofa::core::ObjectFactory;

#ifdef SOFA_DEV
#include <BeamAdapter/component/AdaptiveBeamContactMapper.h>
#include <BeamAdapter/component/MultiAdaptiveBeamContactMapper.h>
#endif // SOFA_DEV

namespace sofa::component
{
	extern "C" {
		SOFA_BEAMADAPTER_API void initExternalModule();
		SOFA_BEAMADAPTER_API const char* getModuleLicense();
		SOFA_BEAMADAPTER_API const char* getModuleName();
		SOFA_BEAMADAPTER_API const char* getModuleVersion();
		SOFA_BEAMADAPTER_API const char* getModuleDescription();
		SOFA_BEAMADAPTER_API const char* getModuleComponentList();
	}

	void initBeamAdapter()
	{
		static bool first = true;
		if (first)
		{
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
		return sofa_tostring(SOFA_TARGET);
	}

	const char* getModuleVersion()
	{
		return sofa_tostring(BEAMADAPTER_VERSION);
	}

	const char* getModuleDescription()
	{
		return "A dynamic adapter that modulates the DOF repartition of a beam model according to its radius of curvature.";
	}

	const char* getModuleComponentList()
	{
		/// string containing the names of the classes provided by the plugin
		static std::string classes = ObjectFactory::getInstance()->listClassesFromTarget(sofa_tostring(SOFA_TARGET));
		return classes.c_str();
	}

} // namespace sofa::component
