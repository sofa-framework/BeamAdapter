/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
#include "AdaptiveBeamLengthConstraint.inl"

#include <sofa/defaulttype/Vec3Types.h>
//#include <BaseMechanics/MechanicalObject.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

namespace _adaptivebeamlengthconstraint_
{

using namespace sofa::defaulttype;
using namespace sofa::helper;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(AdaptiveBeamLengthConstraint)

//TODO(dmarchal 2017-05-17): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int AdaptiveBeamLengthConstraintClass = core::RegisterObject("TODO")
        #ifdef SOFA_WITH_FLOAT
        .add< AdaptiveBeamLengthConstraint<Rigid3fTypes> >()
        #endif
        #ifdef SOFA_WITH_DOUBLE
        .add< AdaptiveBeamLengthConstraint<Rigid3dTypes> >()
        #endif
        ;

#ifdef SOFA_WITH_FLOAT
template class AdaptiveBeamLengthConstraint<Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class AdaptiveBeamLengthConstraint<Rigid3dTypes>;
#endif

} /// namespace _adaptivebeamlengthconstraint_

} /// namespace constraintset

} /// namespace component

} /// namespace sofa

