/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
//
// C++ Implementation : AdaptiveBeamMapping
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/ObjectFactory.h>

#include "AdaptiveBeamMapping.inl"

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
SOFA_DECL_CLASS(AdaptiveBeamMapping)

// Register in the Factory
int AdaptiveBeamMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs")
#ifdef SOFA_WITH_DOUBLE
.add< AdaptiveBeamMapping<Rigid3dTypes, Vec3dTypes   > >(true) //default template
.add< AdaptiveBeamMapping<Rigid3dTypes, Rigid3dTypes > >()
#endif
#ifdef SOFA_WITH_FLOAT
.add< AdaptiveBeamMapping< Rigid3fTypes, Vec3fTypes > >()
#endif

#ifdef SOFA_WITH_FLOAT
#ifdef SOFA_WITH_DOUBLE
.add< AdaptiveBeamMapping< Rigid3dTypes, Vec3fTypes > >()
.add< AdaptiveBeamMapping< Rigid3fTypes, Vec3dTypes > >()
#endif
#endif
;

} // namespace mapping

} // namespace component

} // namespace sofa

