/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: M. Adam, J. Allard, B. Andre, P-J. Bensoussan, S. Cotin, C. Duriez,*
* H. Delingette, F. Falipou, F. Faure, S. Fonteneau, L. Heigeas, C. Mendoza,  *
* M. Nesme, P. Neumann, J-P. de la Plata Alcade, F. Poyer and F. Roy          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <BeamAdapter/component/initBeamAdapter.h>

#include <sofa/core/ObjectFactory.h>
using sofa::core::ObjectFactory;

#ifdef SOFA_DEV
#include <BeamAdapter/component/AdaptiveBeamContactMapper.h>
#include <BeamAdapter/component/MultiAdaptiveBeamContactMapper.h>
#endif // SOFA_DEV

namespace sofa
{

namespace component
{
	extern "C" {
		SOFA_BEAMADAPTER_API void initExternalModule();
		SOFA_BEAMADAPTER_API const char* getModuleLicense();
		SOFA_BEAMADAPTER_API const char* getModuleName();
		SOFA_BEAMADAPTER_API const char* getModuleVersion();
		SOFA_BEAMADAPTER_API const char* getModuleDescription();
		SOFA_BEAMADAPTER_API const char* getModuleComponentList();
	}

	//Here are just several convenient functions to help user to know what contains the plugin
	
	void initExternalModule()
	{
		static bool first = true;
		if (first)
		{
			first = false;
		}
	}
	const char* getModuleLicense()
	{
		return "INRIA and Digital-Trainers";
	}

	const char* getModuleName()
	{
		return "BeamAdapter";
	}

	const char* getModuleVersion()
	{
		return "0.1";
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
} 
} 
