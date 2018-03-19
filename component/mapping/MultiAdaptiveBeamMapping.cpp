/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
//
// C++ Implementation : UnifiedMultiMultiAdaptiveBeamMapping
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

#include "MultiAdaptiveBeamMapping.inl"
#include <sofa/core/ObjectFactory.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/Mapping.inl>

namespace sofa
{

namespace component
{

namespace mapping
{
using namespace defaulttype;
using namespace core;
using namespace core::behavior;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(MultiAdaptiveBeamMapping)

// Register in the Factory
int MultiAdaptiveBeamMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs")
#ifdef SOFA_WITH_FLOAT
.add< MultiAdaptiveBeamMapping<Rigid3fTypes, Vec3fTypes   > >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< MultiAdaptiveBeamMapping< Rigid3dTypes, Vec3dTypes > >()
#endif
#ifdef SOFA_WITH_FLOAT
#ifdef SOFA_WITH_DOUBLE
.add< MultiAdaptiveBeamMapping< Rigid3dTypes, Vec3fTypes > >()
.add< MultiAdaptiveBeamMapping< Rigid3fTypes, Vec3dTypes > >()
#endif
#endif
;


#ifdef SOFA_WITH_FLOAT
    template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamMapping<Rigid3fTypes, Vec3fTypes   >;
#endif
#ifdef SOFA_WITH_DOUBLE
    template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamMapping< Rigid3dTypes, Vec3dTypes >;
#endif

#ifdef SOFA_WITH_FLOAT
#ifdef SOFA_WITH_DOUBLE
    template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamMapping< Rigid3dTypes, Vec3fTypes >;
    template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamMapping< Rigid3fTypes, Vec3dTypes >;
#endif
#endif


} // namespace mapping

} // namespace component

} // namespace sofa

