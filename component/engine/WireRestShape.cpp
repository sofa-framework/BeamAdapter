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
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "WireRestShape.inl"

#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>

namespace sofa
{

namespace component
{

namespace engine
{

namespace _wirerestshape_
{
using namespace sofa::defaulttype;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(WireRestShape)

//TODO(damien 2017-05-17): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int WireRestShapeClass = core::RegisterObject("")
#ifdef SOFA_WITH_FLOAT
.add< WireRestShape<Rigid3fTypes> >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< WireRestShape<Rigid3dTypes> >()
#endif
;

#ifdef SOFA_WITH_FLOAT
template class SOFA_BEAMADAPTER_API WireRestShape<Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class SOFA_BEAMADAPTER_API WireRestShape<Rigid3dTypes>;
#endif

} // namespace _wirerestshape_

}// namespace engine

} // namespace component

} // namespace sofa
