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
#ifndef SOFA_COMPONENT_COLLISION_MULTIADAPTIVEBEAMCONTACTMAPPER_H
#define SOFA_COMPONENT_COLLISION_MULTIADAPTIVEBEAMCONTACTMAPPER_H

#include <sofa/helper/system/config.h>
#include <sofa/helper/Factory.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/simulation/common/Node.h>
#include <sofa/simulation/common/Simulation.h>
#include <SofaBaseCollision/BaseContactMapper.h>
#include <SofaAdvancedConstraint/BSplineModel.h> //ctn_DEV
#include <sofa/core/VecId.h>
#include <iostream>

#include "MultiAdaptiveBeamMapping.h"
#include "WireBeamInterpolation.h"
#include "InterventionalRadiologyController.h"

namespace sofa
{

namespace component
{

namespace collision
{

using namespace sofa::defaulttype;


/*!
 * \class MultiAdaptiveBeamContactMapper
 * @brief MultiAdaptiveBeamContactMapper Class
 */
template < class TCollisionModel, class DataTypes >
class MultiAdaptiveBeamContactMapper : public BaseContactMapper<DataTypes>
{
public:
    typedef typename DataTypes::Real Real;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef TCollisionModel MCollisionModel;
    typedef typename MCollisionModel::InDataTypes InDataTypes;
    typedef core::behavior::MechanicalState<InDataTypes> InMechanicalState;
    typedef core::behavior::MechanicalState<typename MultiAdaptiveBeamContactMapper::DataTypes> MMechanicalState;
    typedef component::container::MechanicalObject<typename MultiAdaptiveBeamContactMapper::DataTypes> MMechanicalObject;
    typedef component::mapping::MultiAdaptiveBeamMapping< InDataTypes, typename MultiAdaptiveBeamContactMapper::DataTypes > MMapping;
    typedef component::controller::InterventionalRadiologyController<InDataTypes> InterventionalRadiologyController;

    MCollisionModel* model;
    //    simulation::Node* child;
    simulation::Node::SPtr child;
    typename MMapping::SPtr mapping;
    MMechanicalState* outmodel;
    InterventionalRadiologyController* m_ircontroller;
    int nbp;

    MultiAdaptiveBeamContactMapper()
    : model(NULL), child(NULL), mapping(NULL), outmodel(NULL), nbp(0)
    {
    }

    void setCollisionModel(MCollisionModel* model)
    {
	this->model = model;
    }

    void cleanup();

    MMechanicalState* createMapping(const char* name="contactPoints");

    void resize(int size)
    {
        if (mapping!=NULL)
            mapping->clear(size);
        nbp = 0;
    }

    //A barycentric coordinate of a spline is define by its curvilinear coordinate and its radius
    //in baryP argument here, the first is curvilinear, second is radius, third is not use for instance
    int addBaryPoint(const Vector3& baryP, int splineId, Real& r)
    {
        int i = nbp++;
        if ((int)outmodel->getX()->size() <= i)
            outmodel->resize(i+1);
        if (mapping)
        {
        	CubicBezierCurve<2> cb(this->model, splineId);
        	int localsplineId;
        	int edgeId = model->comuteEdgeIdFromSplineId(splineId,localsplineId);
        	if ( cb.isCubicBezier() )
        		i = mapping->addBaryPoint(edgeId,baryP,true);
        	else
        		i = mapping->addBaryPoint(edgeId,baryP,false);
            r=cb.r();
        }
        else
        {
            std::cout << "WARNING[MultiAdaptiveBeamContactMapper] AdaptiveBeamMapping is not created" << std::endl;
        }
        return i;
    }



	void update()
    {
        if (mapping!=NULL)
        {
            core::BaseMapping* map = mapping.get();
            map->apply(core::MechanicalParams::defaultInstance(), core::VecCoordId::position(), core::ConstVecCoordId::position());
            map->applyJ(core::MechanicalParams::defaultInstance(), core::VecDerivId::velocity(), core::ConstVecDerivId::velocity());
        }
    }

    void updateXfree()
    {
        if (mapping!=NULL)
        {
            core::BaseMapping* map = mapping.get();
            map->apply(core::MechanicalParams::defaultInstance(), core::VecCoordId::freePosition(), core::ConstVecCoordId::freePosition());
            map->applyJ(core::MechanicalParams::defaultInstance(), core::VecDerivId::freeVelocity(), core::ConstVecDerivId::freeVelocity());
        }
    }
};


/// Mapper for BSplineModel
template<class DataTypes>
class ContactMapper<BSplineModel<2>, DataTypes> : public MultiAdaptiveBeamContactMapper<BSplineModel<2>, DataTypes>
{
public:
	typedef typename DataTypes::Real Real;
	typedef typename DataTypes::Coord Coord;
};




#if defined(WIN32) && !defined(SOFA_BUILD_COMPONENT_COLLISION)

		extern template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamContactMapper<BSplineModel<2>,Vec3Types>;
	

	extern template class SOFA_BEAMADAPTER_API ContactMapper<BSplineModel<2>>;
#endif


} // namespace collision

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_COLLISION_MULTIADAPTIVEBEAMCONTACTMAPPER_H */
