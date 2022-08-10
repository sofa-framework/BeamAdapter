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


namespace sofa::beamadapter
{
    struct BeamSection {
        double _r; 			///< Radius of the beam section
        double _rInner; 	///< Inner radius of the section if beam is hollow
        double _Iy;         ///< Iy is the cross-section moment of inertia (assuming mass ratio = 1) about the y axis, see https ://en.wikipedia.org/wiki/Second_moment_of_area
        double _Iz; 		///< Iz is the cross-section moment of inertia (assuming mass ratio = 1) about the z axis, see https ://en.wikipedia.org/wiki/Second_moment_of_area
        double _J;  		///< Polar moment of inertia (J = Iy + Iz)
        double _A; 			///< A is the cross-sectional area
        double _Asy; 		///< _Asy is the y-direction effective shear area =  10/9 (for solid circular section) or 0 for a non-Timoshenko beam
        double _Asz; 		///< _Asz is the z-direction effective shear area
    };

} // namespace sofa::beamAdapter
