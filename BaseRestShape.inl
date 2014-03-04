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

#ifndef SOFA_COMPONENT_ENGINE_BASERESTSHAPE_INL
#define SOFA_COMPONENT_ENGINE_BASERESTSHAPE_INL
#define VERIF 1

#include "BaseRestShape.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <BaseTopology/EdgeSetGeometryAlgorithms.h>
#include <BaseTopology/QuadSetTopologyModifier.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/TopologyChangeVisitor.h>

#define PI 3.14159265
#define EPSILON 0.0000000001

namespace sofa
{

namespace component
{

namespace engine
{


template<class DataTypes>
void BaseRestShape<DataTypes>::init()
{

	sofa::helper::vector<int> & _density   =*density.beginEdit();
	sofa::helper::vector<Real>& _keyPoints =*keyPoints.beginEdit();

	_density.resize(1)  ; _density[0]  =       10 ;
	_keyPoints.resize(2); _keyPoints[0]= (Real)0.; _keyPoints[1]= length.getValue();

	density.endEdit();
	keyPoints.endEdit();

}


template <class DataTypes>
void BaseRestShape<DataTypes>::releaseWirePart()
{

}


template <class DataTypes>
void BaseRestShape<DataTypes>::getSamplingParameters(helper::vector<Real>& xP_noticeable, helper::vector<int>& nbP_density)
{
	xP_noticeable.resize(keyPoints.getValue().size());
	nbP_density.resize(density.getValue().size());

	nbP_density[0]   = density.getValue()[0];
	xP_noticeable[0] = keyPoints.getValue()[0];
	xP_noticeable[1] = keyPoints.getValue()[1];
}



template <class DataTypes>
void BaseRestShape<DataTypes>::getRestTransformOnX(Transform &global_H_local, const Real &x)
{
	Real x_used = x ;

	if(x_used>length.getValue())
		x_used=length.getValue();

	if(x_used<0.0)
		x_used=0.0;

	if( x_used < length.getValue())
	{
		global_H_local.set(Vec3(x_used, 0.0, 0.0 ), sofa::defaulttype::Quat());
		return;
	}
}


template <class DataTypes>
void BaseRestShape<DataTypes>::getYoungModulusAtX(Real& /*x_curv*/, Real& youngModulus, Real& cPoisson)
{
	youngModulus = _youngModulus1.getValue() ;
	cPoisson     = _poissonRatio.getValue();
}







}// namespace engine


} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_ENGINE_BASERESTSHAPE_INL */
