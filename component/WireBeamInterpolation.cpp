/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
//
// C++ Implementation : WireBeamInterpolation / AdaptiveBeamForceFieldAndMass
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//


#include "initBeamAdapter.h"
#include "WireBeamInterpolation.inl"
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>


namespace sofa
{

namespace component
{

namespace fem
{

namespace _wirebeaminterpolation_
{
using namespace sofa::defaulttype;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(WireBeamInterpolation)

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int WireBeamInterpolationClass = core::RegisterObject("Adaptive Beam Interpolation on Wire rest Shape")
#ifdef SOFA_WITH_FLOAT
.add< WireBeamInterpolation<Rigid3fTypes> >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< WireBeamInterpolation<Rigid3dTypes> >()
#endif
;

#ifdef SOFA_WITH_FLOAT
template class SOFA_BEAMADAPTER_API WireBeamInterpolation<Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class SOFA_BEAMADAPTER_API WireBeamInterpolation<Rigid3dTypes>;
#endif

} // namespace _wirebeaminterpolation_

} // namespace fem

} // namespace component

} // namespace sofa

