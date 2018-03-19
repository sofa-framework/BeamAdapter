/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
//
// C++ Implementation : InterventionalRadiologyController
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "../initBeamAdapter.h"
#include "InterventionalRadiologyController.inl"

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>


namespace sofa
{

namespace component
{


namespace controller
{

using namespace sofa::defaulttype;


/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(InterventionalRadiologyController)

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int InterventionalRadiologyControllerClass = core::RegisterObject("")
#ifdef SOFA_WITH_FLOAT
.add< InterventionalRadiologyController<Rigid3fTypes> >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< InterventionalRadiologyController<Rigid3dTypes> >()
#endif
;

#ifdef SOFA_WITH_FLOAT
template class SOFA_BEAMADAPTER_API InterventionalRadiologyController<Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class SOFA_BEAMADAPTER_API InterventionalRadiologyController<Rigid3dTypes>;
#endif



} // namespace controller

} // namespace component

} // namespace sofa
