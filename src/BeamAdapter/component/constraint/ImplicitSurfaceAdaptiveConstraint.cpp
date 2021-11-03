/******************************************************************************
*                              BeamAdapter plugin                             *
*                  (c) 2006 Inria, University of Lille, CNRS                  *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Authors: see Authors.md                                                     *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "ImplicitSurfaceAdaptiveConstraint.inl"

#include <sofa/defaulttype/VecTypes.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace constraint
{

using namespace sofa::defaulttype;
using namespace sofa::helper;
using core::RegisterObject;

#ifdef SOFAEVE
SOFA_DECL_CLASS(ImplicitSurfaceAdaptiveConstraint)

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int ImplicitSurfaceAdaptiveConstraintClass = RegisterObject("PROUT TODO-ImplicitSurfaceAdaptiveConstraint")
.add< ImplicitSurfaceAdaptiveConstraint<Rigid3Types> >(true)

;

template class ImplicitSurfaceAdaptiveConstraint<Rigid3Types>;


#endif  // SOFAEVE

} // namespace constraint

} // namespace component

} // namespace sofa
