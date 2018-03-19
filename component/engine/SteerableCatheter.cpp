/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
//
//
// Description:
//
//
// Author: Hugo Talbot, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "SteerableCatheter.inl"

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>


namespace sofa
{

namespace component
{

namespace engine
{

using namespace sofa::defaulttype;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(SteerableCatheter)

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int SteerableCatheterClass = core::RegisterObject("")
#ifdef SOFA_WITH_FLOAT
.add< SteerableCatheter<Rigid3fTypes> >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< SteerableCatheter<Rigid3dTypes> >()
#endif
;

#ifdef SOFA_WITH_FLOAT
template class SOFA_BEAMADAPTER_API SteerableCatheter<Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class SOFA_BEAMADAPTER_API SteerableCatheter<Rigid3dTypes>;
#endif


}// namespace engine

} // namespace component

} // namespace sofa
