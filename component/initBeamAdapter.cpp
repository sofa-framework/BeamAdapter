/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
#include "initBeamAdapter.h"
#include <cstring>
#include <string>

#ifdef SOFA_DEV
#include "AdaptiveBeamContactMapper.h"
#include "MultiAdaptiveBeamContactMapper.h"
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
                return "Copyright © Inria, All rights reserved";
	}

	const char* getModuleName()
	{
		return "BeamAdapter";
	}

	const char* getModuleVersion()
	{
                return "1.0";
	}

	const char* getModuleDescription()
	{
		return "A dynamic adapter that modulates the DOF repartition of a beam model according to its radius of curvature.";
	}

	const char* getModuleComponentList()
	{
        return	"AdaptiveBeamConstraint"
                "AdaptiveBeamController"
                "AdaptiveBeamForceFieldAndMass"
                "AdaptiveBeamLengthConstraint"
                "AdaptiveBeamMapping"
                "BeamInterpolation"
                "InterventionalRadiologyController"
                "MultiAdaptiveBeamMapping"
                "SteerableCatheter"
                "SutureController"
                "UnbuiltGenericConstraintSolver"
                "WireBeamInterpolation"
                "WireRestShape";
	}
} 
} 


SOFA_LINK_CLASS(AdaptiveBeamController)
SOFA_LINK_CLASS(AdaptiveBeamForceFieldAndMass)
SOFA_LINK_CLASS(AdaptiveBeamMapping)
SOFA_LINK_CLASS(BeamInterpolation)
SOFA_LINK_CLASS(InterventionalRadiologyController)
SOFA_LINK_CLASS(MultiAdaptiveBeamMapping)
SOFA_LINK_CLASS(SteerableCatheter)
SOFA_LINK_CLASS(WireBeamInterpolation)
SOFA_LINK_CLASS(WireRestShape)
