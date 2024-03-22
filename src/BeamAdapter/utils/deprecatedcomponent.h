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
#pragma once

#include <sofa/core/objectmodel/BaseObject.h>


namespace sofa::component
{
using sofa::core::objectmodel::BaseObject ;
using sofa::core::objectmodel::BaseContext ;
using sofa::core::objectmodel::BaseObjectDescription ;

class DeprecatedComponent : public BaseObject
{
public:
      SOFA_CLASS(DeprecatedComponent, BaseObject) ;

      /// Pre-construction check method called by ObjectFactory.
      template<class T>
      static bool canCreate(T* obj, BaseContext* /*context*/, BaseObjectDescription* /*arg*/)
      {
          obj->serr << "[Deprecated component]: " << obj->getName() ;
          return false;
      }
private:
} ;

} // namespace sofa::component


