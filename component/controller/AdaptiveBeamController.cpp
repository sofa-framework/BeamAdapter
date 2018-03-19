/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
//
// C++ Implementation : AdaptiveBeamController
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
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

#include "../initBeamAdapter.h"
#include "../controller/AdaptiveBeamController.inl"

namespace sofa
{

namespace component
{

namespace controller
{

namespace _adaptivebeamcontroller_
{

using sofa::defaulttype::Rigid3dTypes;
using sofa::defaulttype::Rigid3fTypes;
using core::RegisterObject;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(AdaptiveBeamController)

//TODO(dmarchal 2017-06-01): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int AdaptiveBeamControllerClass = RegisterObject("")
#ifdef SOFA_WITH_FLOAT
.add< AdaptiveBeamController<Rigid3fTypes> >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< AdaptiveBeamController<Rigid3dTypes> >()
#endif
;

#ifdef SOFA_WITH_FLOAT
template class SOFA_BEAMADAPTER_API AdaptiveBeamController<Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class SOFA_BEAMADAPTER_API AdaptiveBeamController<Rigid3dTypes>;
#endif

} // namespace _adaptivebeamcontroller_

} // namespace controller

} // namespace component

} // namespace sofa
