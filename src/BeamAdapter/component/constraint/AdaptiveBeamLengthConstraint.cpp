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
#include <BeamAdapter/component/constraint/AdaptiveBeamLengthConstraint.inl>

#include <sofa/defaulttype/VecTypes.h>
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace constraintset
{

namespace _adaptivebeamlengthconstraint_
{


void AdaptiveBeamLengthConstraintResolution::init(int line, double** w, double* force)
{
    SOFA_UNUSED(w);

    if(m_initF)
        force[line] = *m_initF;
}
void AdaptiveBeamLengthConstraintResolution::resolution(int line, double** w, double* d, double* force)
{
    force[line] -= d[line] / w[line][line];
    if(force[line] < 0)
        force[line] = 0;
}

void AdaptiveBeamLengthConstraintResolution::store(int line, double* force, bool convergence)
{
    SOFA_UNUSED(convergence) ;

    if(m_initF)
        *m_initF = force[line];
    if(m_active)
        *m_active = (force[line] != 0);
}


using namespace sofa::defaulttype;
using namespace sofa::helper;
using core::RegisterObject;

/////////////////////////////////////////// FACTORY ////////////////////////////////////////////////
///
/// Register the component into the sofa factory.
/// For more details:
/// https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/the-objectfactory/
///
////////////////////////////////////////////////////////////////////////////////////////////////////
static int AdaptiveBeamLengthConstraintClass = RegisterObject("Constrain the length of a beam.")
                .add< AdaptiveBeamLengthConstraint<Rigid3Types> >(true) // default template

        ;

template class AdaptiveBeamLengthConstraint<Rigid3Types>;


} /// namespace _adaptivebeamlengthconstraint_

} /// namespace constraintset

} /// namespace component

} /// namespace sofa

