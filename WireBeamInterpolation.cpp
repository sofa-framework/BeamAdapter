/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
 *                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
 *                                                                             *
 * This library is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU Lesser General Public License as published by    *
 * the Free Software Foundation; either version 2.1 of the License, or (at     *
 * your option) any later version.                                             *
 *                                                                             *
 * This library is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
 * for more details.                                                           *
 *                                                                             *
 * You should have received a copy of the GNU Lesser General Public License    *
 * along with this library; if not, write to the Free Software Foundation,     *
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
 *******************************************************************************
 *                               SOFA :: Modules                               *
 *                                                                             *
 * Authors: The SOFA Team and external contributors (see Authors.txt)          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
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

