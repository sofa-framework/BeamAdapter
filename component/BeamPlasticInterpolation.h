/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
//
// C++ Implementation : BeamPlasticInterpolation
//
// Description: Implementation of plasticity over the AdaptiveBeamForceFieldAndMass
// force field interface.
//
//
// Author: Camille Krewcun, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#pragma once

#include "BeamInterpolation.h"

namespace sofa::plugin::beamadapter::component::forcefield
{

namespace _beamplasticinterpolation_
{

struct BeamGeometry
{
    std::string sectionShape;   ///<cross-section shape (rectangular, square, elliptic, circular)
    double _L;                  ///<length of the beam
    double _Ly;                 ///<length of the beam section along Y (if rectangular)
    double _Lz;                 ///<length of the beam section along Z (if rectangular)
    double _r; 			        ///<radius of the section (if circular)
    double _rInner;		        ///<inner radius of the section if beam is hollow
    double _rSmall;             ///<small radius of the section (if elliptic)
    double _rLarge;             ///<large radius of the section (if elliptic)
};

using sofa::helper::vector;
using sofa::component::fem::BeamInterpolation;
using sofa::plugin::beamadapter::component::forcefield::_beamplasticinterpolation_::BeamGeometry;

template<class DataTypes>
class BeamPlasticInterpolation : public virtual BeamInterpolation<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(BeamPlasticInterpolation, DataTypes),
        SOFA_TEMPLATE(BeamInterpolation, DataTypes));

    BeamPlasticInterpolation();
    virtual ~BeamPlasticInterpolation() {}

    const vector<BeamGeometry> getBeamGeometryParameters();

};

} // namespace _beamplasticinterpolation_

} // namespace sofa::plugin::beamadapter::component::forcefield