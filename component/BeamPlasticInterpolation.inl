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

#include "BeamPlasticInterpolation.h"

namespace sofa::plugin::beamadapter::component::forcefield
{

namespace _beamplasticinterpolation_
{

template <class DataTypes>
BeamPlasticInterpolation<DataTypes>::BeamPlasticInterpolation()
{

}

template <class DataTypes>
const vector<BeamGeometry> BeamPlasticInterpolation<DataTypes>::getBeamGeometryParameters()
{
    const vector<double>& lengthList = d_lengthList.getValue();
    vector<BeamGeometry> beamGeometryParams;

    for (unsigned int i = 0; i < lengthList.size(); i++)
    {
        double L = lengthList[i];
        if (crossSectionShape.getValue().getSelectedItem() == "rectangular")
        {
            BeamGeometry newBeamGeometry = {
                "rectangular",          //cross-section shape
                L,                      //length
                d_lengthY.getValue(),   //section length along Y (if rectangular)
                d_lengthZ.getValue(),   //section length along Z (if rectangular)
                0,                      //radius (if circular)
                0,                      //inner radius (if hollow)
                0,                      //small radius (if elliptic)
                0                       //large radius (if elliptic)
            };
            beamGeometryParams.push_back(newBeamGeometry);
        }
        else if (crossSectionShape.getValue().getSelectedItem() == "square")
        {
            BeamGeometry newBeamGeometry = {
                "square",                //cross-section shape
                L,                       //length
                d_sideLength.getValue(), //section length along Y (if rectangular)
                d_sideLength.getValue(), //section length along Z (if rectangular)
                0,                       //radius (if circular)
                0,                       //inner radius (if hollow)
                0,                       //small radius (if elliptic)
                0                        //large radius (if elliptic)
            };
            beamGeometryParams.push_back(newBeamGeometry);
        }
        else if (crossSectionShape.getValue().getSelectedItem() == "elliptic")
        {
            BeamGeometry newBeamGeometry = {
                "circular",                 //cross-section shape
                L,                          //length
                0,                          //section length along Y (if rectangular)
                0,                          //section length along Z (if rectangular)
                d_radius.getValue(),        //radius (if circular)
                d_innerRadius.getValue(),   //inner radius (if hollow)
                0,                          //small radius (if elliptic)
                0                           //large radius (if elliptic)
            };
            beamGeometryParams.push_back(newBeamGeometry);
        }
        else if (crossSectionShape.getValue().getSelectedItem() == "elliptic")
        {
            BeamGeometry newBeamGeometry = {
                "elliptic",                 //cross-section shape
                L,                          //length
                0,                          //section length along Y (if rectangular)
                0,                          //section length along Z (if rectangular)
                0,                          //radius (if circular)
                0,                          //inner radius (if hollow)
                d_smallRadius.getValue(),   //small radius (if elliptic)
                d_largeRadius.getValue()    //large radius (if elliptic)
            };
            beamGeometryParams.push_back(newBeamGeometry);
        }
        else
        {
            //TO DO: implement quadrature method for a disc and a hollow-disc cross section
            msg_error() << "Quadrature method for " << crossSectionShape.getValue().getSelectedItem()
                << " shape cross section has not been implemented yet. Methods for rectangular cross sections are available";
        }
    }

    return beamGeometryParams;
}

} // namespace _beamplasticinterpolation_

} // namespace sofa::plugin::beamadapter::component::forcefield