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

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

#define SOFA_PLUGIN_BEAMADAPTER_ADAPTVEBEAMFORCEFIELD_CPP
#include "AdaptiveBeamForceFieldAndMass.inl"

namespace sofa
{

namespace component
{

namespace forcefield
{

namespace _adaptivebeamforcefieldandmass_
{

using sofa::core::RegisterObject ;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(AdaptiveBeamForceFieldAndMass)

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int AdaptiveBeamForceFieldAndMassClass = RegisterObject("Adaptive Beam finite elements")
#ifdef SOFA_WITH_FLOAT
.add< AdaptiveBeamForceFieldAndMass<sofa::defaulttype::Rigid3fTypes> >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< AdaptiveBeamForceFieldAndMass<sofa::defaulttype::Rigid3dTypes> >()
#endif
;

#ifdef SOFA_WITH_FLOAT
template class SOFA_BEAMADAPTER_API AdaptiveBeamForceFieldAndMass<sofa::defaulttype::Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class SOFA_BEAMADAPTER_API AdaptiveBeamForceFieldAndMass<sofa::defaulttype::Rigid3dTypes>;
#endif
////////////////////////////////////////////////////////////////////////////////////////////////////

} /// _adaptivebeamforcefiedlandmass_

} /// namespace forcefield

} /// namespace component

} /// namespace sofa

