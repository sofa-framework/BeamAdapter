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
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "BaseRestShape.inl"

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

SOFA_DECL_CLASS(BaseRestShape)


//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int BaseRestShapeClass = core::RegisterObject("Straignt Rest Shape")
#ifndef SOFA_FLOAT
.add< BaseRestShape<Rigid3dTypes> >()
#endif
#ifndef SOFA_DOUBLE
.add< BaseRestShape<Rigid3fTypes> >()
#endif
.addAlias("StraightRestShape")
;

#ifndef SOFA_FLOAT
template class SOFA_BEAMADAPTER_API BaseRestShape<Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
template class SOFA_BEAMADAPTER_API BaseRestShape<Rigid3fTypes>;
#endif


}// namespace engine

} // namespace component

} // namespace sofa
