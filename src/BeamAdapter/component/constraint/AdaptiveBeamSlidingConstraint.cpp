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
#define SOFA_PLUGIN_BEAMADAPTER_ADAPTIVEBEAMSLIDINGCONSTRAINT_CPP

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/component/statecontainer/MechanicalObject.h>
#include <BeamAdapter/config.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/defaulttype/VecTypes.h>

#include <BeamAdapter/component/constraint/AdaptiveBeamSlidingConstraint.inl>


using sofa::core::objectmodel::BaseObject ;
using sofa::core::objectmodel::BaseContext ;
using sofa::core::objectmodel::BaseObjectDescription ;
#include <BeamAdapter/utils/deprecatedcomponent.h>
using sofa::component::DeprecatedComponent;
using sofa::defaulttype::Rigid3Types;
using sofa::defaulttype::Rigid3Types;


namespace sofa::component::constraintset::_adaptiveBeamSlidingConstraint_
{

AdaptiveBeamSlidingConstraintResolution::AdaptiveBeamSlidingConstraintResolution(double* sliding)
    : ConstraintResolution(3)
     ,m_slidingDisp(sliding)
{
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

template class SOFA_BEAMADAPTER_API AdaptiveBeamSlidingConstraint<Rigid3Types>;


///////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace sofa::component::constraintset::_adaptiveBeamSlidingConstraint_

namespace beamadapter
{

void registerAdaptiveBeamSlidingConstraint(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(sofa::core::ObjectRegistrationData("Constrain a rigid to be attached to a beam (only in position, not the orientation).")
                             .add< sofa::component::constraintset::_adaptiveBeamSlidingConstraint_::AdaptiveBeamSlidingConstraint<Rigid3Types> >());
}

} // namespace beamadapter
