/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
//
// C++ Implementation : BeamInterpolation / AdaptiveBeamForceFieldAndMass
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

/// This define is here to prevent the declaration of the template instances as "extern".
/// Have a look a the end of BeamInterpolation.h
#define SOFA_BEAMINTERPOLATION_CPP
#include "BeamInterpolation.inl"

namespace sofa
{

namespace component
{

namespace fem
{

namespace _beaminterpolation_
{

using namespace sofa::defaulttype;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(BeamInterpolation)

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int BeamInterpolationClass = core::RegisterObject("Adaptive Beam Interpolation")
#ifdef SOFA_WITH_FLOAT
.add< BeamInterpolation<Rigid3fTypes> >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< BeamInterpolation<Rigid3dTypes> >()
#endif
;



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Explicit template instanciation of extern template.
////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef SOFA_WITH_FLOAT
template class SOFA_BEAMADAPTER_API BeamInterpolation<Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class SOFA_BEAMADAPTER_API BeamInterpolation<Rigid3dTypes>;
#endif

} /// namespace _beaminterpolation_

} /// namespace fem

} /// namespace component

} /// namespace sofa

