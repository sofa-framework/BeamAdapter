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
// C++ Implementation : BeamLengthMapping
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
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/ObjectFactory.h>

#include "BeamLengthMapping.inl"

namespace sofa
{

namespace component
{

namespace mapping
{

//using namespace defaulttype;
using namespace core;
using namespace core::behavior;
using namespace sofa::defaulttype;


/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
//SOFA_DECL_CLASS(BeamLengthMapping)

// Register in the Factory
int BeamLengthMappingClass = core::RegisterObject("computes the lengths of the beams")
        .add< BeamLengthMapping<Rigid3Types, Vec1dTypes   > >(true) //default template
        //.add< BeamLengthMapping<Rigid3Types, Rigid3Types > >()
;


//int BeamLengthMappingClass = core::RegisterObject("computes the lengths of the beams")
//#ifdef SOFA_WITH_DOUBLE
//.add< BeamLengthMapping<Rigid3dTypes, Vec1dTypes   > >(true) //default template
//#endif

//#ifdef SOFA_WITH_FLOAT
//.add< BeamLengthMapping< Rigid3fTypes, Vec1fTypes > >()
//#endif


//#ifdef SOFA_WITH_FLOAT
//#ifdef SOFA_WITH_DOUBLE
//.add< BeamLengthMapping< Rigid3dTypes, Vec1fTypes > >()
//.add< BeamLengthMapping< Rigid3fTypes, Vec1dTypes > >()
//#endif
//#endif

//;

} // namespace mapping

} // namespace component

} // namespace sofa
