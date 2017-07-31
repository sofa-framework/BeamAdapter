/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU General Public License as published by the Free  *
* Software Foundation; either version 2 of the License, or (at your option)   *
* any later version.                                                          *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
* more details.                                                               *
*                                                                             *
* You should have received a copy of the GNU General Public License along     *
* with this program; if not, write to the Free Software Foundation, Inc., 51  *
* Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
*******************************************************************************
*                            SOFA :: Applications                             *
*                                                                             *
* Authors: M. Adam, J. Allard, B. Andre, P-J. Bensoussan, S. Cotin, C. Duriez,*
* H. Delingette, F. Falipou, F. Faure, S. Fonteneau, L. Heigeas, C. Mendoza,  *
* M. Nesme, P. Neumann, J-P. de la Plata Alcade, F. Poyer and F. Roy          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include "GeometricalUtils.h"

namespace sofa
{

namespace helper
{


template<class Real>
void computeSplineLength(Real &length, const Vec<3, Real>& P0, const Vec<3, Real>& P1, const Vec<3, Real>& P2, const Vec<3, Real> &P3)
{

        // the computation of integral Int[0,1] ||dP(x)||  dx = length
        // is done using Gauss Points

        // definition of the Gauss points
        Real x1, x2, x3, x4;
        Real A = 2*sqrt(6.0/5.0);
        x1 = -sqrt((3.0 - A)/7.0 )/2.0+ 0.5;
        x2 = sqrt((3.0 - A) /7.0 )/2.0+ 0.5;
        x3 = -sqrt((3.0 + A)/7.0 )/2.0+ 0.5;
        x4 = sqrt((3.0 + A) /7.0 )/2.0+ 0.5;

        Vec<3, Real> dP1, dP2, dP3, dP4;



        dP1 = P0*(-3*(1-x1)*(1-x1)) + P1*(3-12*x1+9*x1*x1) + P2*(6*x1-9*x1*x1) + P3*(3*x1*x1);
        dP2 = P0*(-3*(1-x2)*(1-x2)) + P1*(3-12*x2+9*x2*x2) + P2*(6*x2-9*x2*x2) + P3*(3*x2*x2);
        dP3 = P0*(-3*(1-x3)*(1-x3)) + P1*(3-12*x3+9*x3*x3) + P2*(6*x3-9*x3*x3) + P3*(3*x3*x3);
        dP4 = P0*(-3*(1-x4)*(1-x4)) + P1*(3-12*x4+9*x4*x4) + P2*(6*x4-9*x4*x4) + P3*(3*x4*x4);

        // formula with 4 Gauss Points
        Real B= sqrt(30.0);
        length = ((18.0 + B) /72.0 )*dP1.norm() + ((18.0 + B) /72.0 )*dP2.norm() + ((18.0 - B) /72.0 )*dP3.norm() + ((18.0 - B) /72.0 )*dP4.norm();





}


} 
} 




