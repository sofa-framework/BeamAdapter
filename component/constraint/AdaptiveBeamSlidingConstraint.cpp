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
//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/Vec3Types.h>

#include "AdaptiveBeamSlidingConstraint.inl"

namespace sofa
{

namespace component
{

namespace constraintset
{

namespace _AdaptiveBeamSlidingConstraint_
{

AdaptiveBeamSlidingConstraintResolution::AdaptiveBeamSlidingConstraintResolution(double* sliding)
    : m_slidingDisp(sliding)
{
    nbLines = 3;
}

void AdaptiveBeamSlidingConstraintResolution::resolution(int line, double** w, double* d, double* force, double* dfree)
{
    SOFA_UNUSED(dfree);
    resolution(line,w,d,force);
}

void AdaptiveBeamSlidingConstraintResolution::resolution(int line, double** w, double* d, double* force)
{
    double f[2] = {force[line], force[line+1]};

    force[line] -= d[line] / w[line][line];
    d[line+1] += w[line+1][line] * (force[line]-f[0]);
    force[line+1] -= d[line+1] / w[line+1][line+1];
    d[line+2] += w[line+2][line] * (force[line]-f[0]) + w[line+2][line+1] * (force[line+1]-f[1]);
}

void AdaptiveBeamSlidingConstraintResolution::init(int line, double** w, double* force)
{
    SOFA_UNUSED(force);
    m_slidingW = w[line+2][line+2];
}

void AdaptiveBeamSlidingConstraintResolution::store(int line, double* force, bool convergence)
{
    SOFA_UNUSED(convergence);
    if(m_slidingDisp)
        *m_slidingDisp = force[line+2] * m_slidingW;
}


/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
SOFA_DECL_CLASS(AdaptiveBeamSlidingConstraint)

int AdaptiveBeamSlidingConstraintClass = core::RegisterObject("Constrain a rigid to be attached to a beam (only in position, not the orientation)")
        #ifdef SOFA_WITH_FLOAT
        .add< AdaptiveBeamSlidingConstraint<sofa::defaulttype::Rigid3fTypes> >()
        #endif
        #ifdef SOFA_WITH_DOUBLE
        .add< AdaptiveBeamSlidingConstraint<sofa::defaulttype::Rigid3dTypes> >()
        #endif
        ;

#ifdef SOFA_WITH_FLOAT
template class AdaptiveBeamSlidingConstraint<sofa::defaulttype::Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class AdaptiveBeamSlidingConstraint<sofa::defaulttype::Rigid3dTypes>;
#endif
///////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace _AdaptiveBeamSlidingConstraint_

} // namespace constraintset

} // namespace component

} // namespace sofa

