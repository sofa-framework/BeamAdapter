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
#ifndef SOFA_COMPONENT_COLLISION_MULTIADAPTIVEBEAMCONTACTMAPPER_INL
#define SOFA_COMPONENT_COLLISION_MULTIADAPTIVEBEAMCONTACTMAPPER_INL

#include "MultiAdaptiveBeamContactMapper.h"
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/Simulation.h>
#include <sofa/simulation/common/DeleteVisitor.h>
#include <iostream>

namespace sofa
{

namespace component
{

namespace collision
{


template < class TCollisionModel, class DataTypes >
void MultiAdaptiveBeamContactMapper<TCollisionModel,DataTypes>::cleanup()
{
    if (mapping!=NULL)
    {
    simulation::Node* parent = dynamic_cast<simulation::Node*>(model->getContext());
        if (parent!=NULL)
        {
            simulation::Node::SPtr child = dynamic_cast<simulation::Node*>(mapping->getContext());
            child->detachFromGraph();
            child->execute<simulation::DeleteVisitor>(sofa::core::ExecParams::defaultInstance());
            child.reset();
            mapping.reset();
        }
    }
}

template < class TCollisionModel, class DataTypes >
typename MultiAdaptiveBeamContactMapper<TCollisionModel,DataTypes>::MMechanicalState* MultiAdaptiveBeamContactMapper<TCollisionModel,DataTypes>::createMapping(const char* name)
{
    if (model==NULL) return NULL;
    InMechanicalState* instate = model->getMechanicalState();
    if (instate!=NULL)
    {
        simulation::Node* parent = dynamic_cast<simulation::Node*>(instate->getContext());
        if (parent==NULL )
        {
            std::cerr << "ERROR: MultiAdaptiveBeamContactMapper only works for scenegraph scenes.\n";
            return NULL;
        }
        instate->getContext()->get(m_ircontroller);
        if (m_ircontroller==NULL )
        {
            std::cerr << "ERROR: MultiAdaptiveBeamContactMapper only works if having InterventionalRadiologyController.\n";
            return NULL;
        }

        simulation::Node::SPtr childPtr=  parent->createChild(name);
        child =childPtr.get();
        parent->addChild(child); child->updateSimulationContext();
        typename MMechanicalObject::SPtr outmodel = sofa::core::objectmodel::New<MMechanicalObject>();
        child->addObject(outmodel);
        outmodel->useMask.setValue(true);

        mapping =  sofa::core::objectmodel::New<MMapping>(instate, outmodel.get(),m_ircontroller);
        mapping->setName(name);
        mapping->setBarycentricMapping();
        child->addObject(mapping);

    }
    else
    {
        simulation::Node* parent = dynamic_cast<simulation::Node*>(model->getContext());
        if (parent==NULL)
        {
            std::cerr << "ERROR: MultiAdaptiveBeamContactMapper only works for scenegraph scenes.\n";
            return NULL;
        }
        simulation::Node::SPtr childPtr=  parent->createChild(name);
        child =childPtr.get();
        parent->addChild(child); child->updateSimulationContext();
        typename MMechanicalObject::SPtr outmodel = sofa::core::objectmodel::New<MMechanicalObject>();
        outmodel->useMask.setValue(true);
        mapping = NULL;
    }
    return outmodel;
}


} // namespace collision

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_COLLISION_MULTIADAPTIVEBEAMCONTACTMAPPER_INL */
