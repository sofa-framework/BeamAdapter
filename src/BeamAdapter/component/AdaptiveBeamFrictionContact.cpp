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
/*
 *  Specialization of Friction Contact Reponse for BSplineModel collision
 *
 * */


#ifdef SOFA_DEV


#include <sofa/component/collision/FrictionContact.inl>

#include "AdaptiveBeamContactMapper.h"
#include "MultiAdaptiveBeamContactMapper.h"



namespace sofa
{

namespace component
{

namespace collision
{
//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
template < >
FrictionContact<BSplineModel<1> , PointModel>::FrictionContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod)
: model1(model1)
, model2(model2)
, intersectionMethod(intersectionMethod)
, m_constraint(NULL)
, parent(NULL)
, mu (initData(&mu, 0.8, "mu", "friction coefficient (0 for frictionless contacts)"))
{
    selfCollision = ((core::CollisionModel*)model1 == (core::CollisionModel*)model2);
    mapper1.setCollisionModel(model1);
    if (!selfCollision) mapper2.setCollisionModel(model2);
    contacts.clear();
    mappedContacts.clear();

}
template < >
FrictionContact<BSplineModel<1> , SphereModel>::FrictionContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod)
: model1(model1)
, model2(model2)
, intersectionMethod(intersectionMethod)
, m_constraint(NULL)
, parent(NULL)
, mu (initData(&mu, 0.8, "mu", "friction coefficient (0 for frictionless contacts)"))
{
    selfCollision = ((core::CollisionModel*)model1 == (core::CollisionModel*)model2);
    mapper1.setCollisionModel(model1);
    if (!selfCollision) mapper2.setCollisionModel(model2);
    contacts.clear();
    mappedContacts.clear();

}





template <>
void FrictionContact<BSplineModel<1> , PointModel>::activateMappers()
{
    if (!m_constraint)
	{
        // Get the mechanical model from mapper1 to fill the constraint vector
        MechanicalState1* mmodel1 = mapper1.createMapping();
        // Get the mechanical model from mapper2 to fill the constraints vector
        MechanicalState2* mmodel2 = selfCollision ? mmodel1 : mapper2.createMapping();
        m_constraint = sofa::core::objectmodel::New<constraintset::UnilateralInteractionConstraint<Vec3Types> >(mmodel1, mmodel2);
        m_constraint->setName( getName() );
    }

    int size = contacts.size();
    m_constraint->clear(size);
    if (selfCollision)
		mapper1.resize(2*size);
    else
    {
        mapper1.resize(size);
        mapper2.resize(size);
    }
    int i = 0;
    const double d0 = intersectionMethod->getContactDistance() + model1->getProximity() + model2->getProximity(); // - 0.001;

    //std::cout<<" d0 = "<<d0<<std::endl;

    mappedContacts.resize(contacts.size());
    for (std::vector<DetectionOutput*>::const_iterator it = contacts.begin(); it!=contacts.end(); it++, i++)
    {
        DetectionOutput* o = *it;
        //std::cout<<" collisionElements :"<<o->elem.first<<" - "<<o->elem.second<<std::endl;
        CollisionElement1 elem1(o->elem.first);
        CollisionElement2 elem2(o->elem.second);
        int index1 = elem1.getIndex();
        int index2 = elem2.getIndex();
        //std::cout<<" indices :"<<index1<<" - "<<index2<<std::endl;

        DataTypes1::Real r1 = 0.;
        DataTypes2::Real r2 = 0.;
        //double constraintValue = ((o->point[1] - o->point[0]) * o->normal) - intersectionMethod->getContactDistance();

        // Create mapping for first point
        index1 = mapper1.addBaryPoint(o->baryCoords[0], index1, r1);
        // Create mapping for second point
        index2 = selfCollision ? mapper1.addPoint(o->point[1], index2, r2) : mapper2.addPoint(o->point[1], index2, r2);
        double distance = d0 + r1 + r2;

        mappedContacts[i].first.first = index1;
        mappedContacts[i].first.second = index2;
        mappedContacts[i].second = distance;
    }

    // Update mappings
    mapper1.update();
    mapper1.updateXfree();
    if (!selfCollision) mapper2.update();
    if (!selfCollision) mapper2.updateXfree();
	//std::cerr<<" end activateMappers call"<<std::endl;

}

template <>
void FrictionContact<BSplineModel<1> , SphereModel>::activateMappers()
{
    if (!m_constraint)
	{
        // Get the mechanical model from mapper1 to fill the constraint vector
        MechanicalState1* mmodel1 = mapper1.createMapping();
        // Get the mechanical model from mapper2 to fill the constraints vector
        MechanicalState2* mmodel2 = selfCollision ? mmodel1 : mapper2.createMapping();
        m_constraint = sofa::core::objectmodel::New<constraintset::UnilateralInteractionConstraint<Vec3Types> >(mmodel1, mmodel2);
        m_constraint->setName( getName() );
    }

    int size = contacts.size();
    m_constraint->clear(size);
    if (selfCollision)
		mapper1.resize(2*size);
    else
    {
        mapper1.resize(size);
        mapper2.resize(size);
    }
    int i = 0;
    const double d0 = intersectionMethod->getContactDistance() + model1->getProximity() + model2->getProximity(); // - 0.001;

    //std::cout<<" d0 = "<<d0<<std::endl;

    mappedContacts.resize(contacts.size());
    for (std::vector<DetectionOutput*>::const_iterator it = contacts.begin(); it!=contacts.end(); it++, i++)
    {
        DetectionOutput* o = *it;
        //std::cout<<" collisionElements :"<<o->elem.first<<" - "<<o->elem.second<<std::endl;
        CollisionElement1 elem1(o->elem.first);
        CollisionElement2 elem2(o->elem.second);
        int index1 = elem1.getIndex();
        int index2 = elem2.getIndex();
        //std::cout<<" indices :"<<index1<<" - "<<index2<<std::endl;

        DataTypes1::Real r1 = 0.;
        DataTypes2::Real r2 = 0.;
        //double constraintValue = ((o->point[1] - o->point[0]) * o->normal) - intersectionMethod->getContactDistance();

        // Create mapping for first point
        index1 = mapper1.addBaryPoint(o->baryCoords[0], index1, r1);
        // Create mapping for second point
        index2 = selfCollision ? mapper1.addPoint(o->point[1], index2, r2) : mapper2.addPoint(o->point[1], index2, r2);
        double distance = d0 + r1 + r2;

        mappedContacts[i].first.first = index1;
        mappedContacts[i].first.second = index2;
        mappedContacts[i].second = distance;
    }

    // Update mappings
    mapper1.update();
    mapper1.updateXfree();
    if (!selfCollision) mapper2.update();
    if (!selfCollision) mapper2.updateXfree();
	//std::cerr<<" end activateMappers call"<<std::endl;

}

Creator<Contact::Factory, FrictionContact<BSplineModel<1> , PointModel> > AdaptiveBSplinePointFrictionContactClass("FrictionContact", true);
Creator<Contact::Factory, FrictionContact<BSplineModel<1> , SphereModel> > AdaptiveBSplineSphereFrictionContactClass("FrictionContact", true);




//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
template < >
FrictionContact<BSplineModel<2> , PointModel>::FrictionContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod)
: model1(model1)
, model2(model2)
, intersectionMethod(intersectionMethod)
, m_constraint(NULL)
, parent(NULL)
, mu (initData(&mu, 0.8, "mu", "friction coefficient (0 for frictionless contacts)"))
{
    selfCollision = ((core::CollisionModel*)model1 == (core::CollisionModel*)model2);
    mapper1.setCollisionModel(model1);
    if (!selfCollision) mapper2.setCollisionModel(model2);
    contacts.clear();
    mappedContacts.clear();

}
template < >
FrictionContact<BSplineModel<2> , SphereModel>::FrictionContact(CollisionModel1* model1, CollisionModel2* model2, Intersection* intersectionMethod)
: model1(model1)
, model2(model2)
, intersectionMethod(intersectionMethod)
, m_constraint(NULL)
, parent(NULL)
, mu (initData(&mu, 0.8, "mu", "friction coefficient (0 for frictionless contacts)"))
{
    selfCollision = ((core::CollisionModel*)model1 == (core::CollisionModel*)model2);
    mapper1.setCollisionModel(model1);
    if (!selfCollision) mapper2.setCollisionModel(model2);
    contacts.clear();
    mappedContacts.clear();

}





template <>
void FrictionContact<BSplineModel<2> , PointModel>::activateMappers()
{
    if (!m_constraint)
	{
        // Get the mechanical model from mapper1 to fill the constraint vector
        MechanicalState1* mmodel1 = mapper1.createMapping();
        // Get the mechanical model from mapper2 to fill the constraints vector
        MechanicalState2* mmodel2 = selfCollision ? mmodel1 : mapper2.createMapping();
        m_constraint = sofa::core::objectmodel::New<constraintset::UnilateralInteractionConstraint<Vec3Types> >(mmodel1, mmodel2);
        m_constraint->setName( getName() );
    }

    int size = contacts.size();
    m_constraint->clear(size);
    if (selfCollision)
		mapper1.resize(2*size);
    else
    {
        mapper1.resize(size);
        mapper2.resize(size);
    }
    int i = 0;
    const double d0 = intersectionMethod->getContactDistance() + model1->getProximity() + model2->getProximity(); // - 0.001;

    //std::cout<<" d0 = "<<d0<<std::endl;

    mappedContacts.resize(contacts.size());
    for (std::vector<DetectionOutput*>::const_iterator it = contacts.begin(); it!=contacts.end(); it++, i++)
    {
        DetectionOutput* o = *it;
        //std::cout<<" collisionElements :"<<o->elem.first<<" - "<<o->elem.second<<std::endl;
        CollisionElement1 elem1(o->elem.first);
        CollisionElement2 elem2(o->elem.second);
        int index1 = elem1.getIndex();
        int index2 = elem2.getIndex();
        //std::cout<<" indices :"<<index1<<" - "<<index2<<std::endl;

        DataTypes1::Real r1 = 0.;
        DataTypes2::Real r2 = 0.;
        //double constraintValue = ((o->point[1] - o->point[0]) * o->normal) - intersectionMethod->getContactDistance();

        // Create mapping for first point
        index1 = mapper1.addBaryPoint(o->baryCoords[0], index1, r1);
        // Create mapping for second point
        index2 = selfCollision ? mapper1.addPoint(o->point[1], index2, r2) : mapper2.addPoint(o->point[1], index2, r2);
        double distance = d0 + r1 + r2;

        mappedContacts[i].first.first = index1;
        mappedContacts[i].first.second = index2;
        mappedContacts[i].second = distance;
    }

    // Update mappings
    mapper1.update();
    mapper1.updateXfree();
    if (!selfCollision) mapper2.update();
    if (!selfCollision) mapper2.updateXfree();
	//std::cerr<<" end activateMappers call"<<std::endl;

}

template <>
void FrictionContact<BSplineModel<2> , SphereModel>::activateMappers()
{
    if (!m_constraint)
	{
        // Get the mechanical model from mapper1 to fill the constraint vector
        MechanicalState1* mmodel1 = mapper1.createMapping();
        // Get the mechanical model from mapper2 to fill the constraints vector
        MechanicalState2* mmodel2 = selfCollision ? mmodel1 : mapper2.createMapping();
        m_constraint = sofa::core::objectmodel::New<constraintset::UnilateralInteractionConstraint<Vec3Types> >(mmodel1, mmodel2);
        m_constraint->setName( getName() );
    }

    int size = contacts.size();
    m_constraint->clear(size);
    if (selfCollision)
		mapper1.resize(2*size);
    else
    {
        mapper1.resize(size);
        mapper2.resize(size);
    }
    int i = 0;
    const double d0 = intersectionMethod->getContactDistance() + model1->getProximity() + model2->getProximity(); // - 0.001;

    //std::cout<<" d0 = "<<d0<<std::endl;

    mappedContacts.resize(contacts.size());
    for (std::vector<DetectionOutput*>::const_iterator it = contacts.begin(); it!=contacts.end(); it++, i++)
    {
        DetectionOutput* o = *it;
        //std::cout<<" collisionElements :"<<o->elem.first<<" - "<<o->elem.second<<std::endl;
        CollisionElement1 elem1(o->elem.first);
        CollisionElement2 elem2(o->elem.second);
        int index1 = elem1.getIndex();
        int index2 = elem2.getIndex();
        //std::cout<<" indices :"<<index1<<" - "<<index2<<std::endl;

        DataTypes1::Real r1 = 0.;
        DataTypes2::Real r2 = 0.;
        //double constraintValue = ((o->point[1] - o->point[0]) * o->normal) - intersectionMethod->getContactDistance();

        // Create mapping for first point
        index1 = mapper1.addBaryPoint(o->baryCoords[0], index1, r1);
        // Create mapping for second point
        index2 = selfCollision ? mapper1.addPoint(o->point[1], index2, r2) : mapper2.addPoint(o->point[1], index2, r2);
        double distance = d0 + r1 + r2;

        mappedContacts[i].first.first = index1;
        mappedContacts[i].first.second = index2;
        mappedContacts[i].second = distance;
    }

    // Update mappings
    mapper1.update();
    mapper1.updateXfree();
    if (!selfCollision) mapper2.update();
    if (!selfCollision) mapper2.updateXfree();
	//std::cerr<<" end activateMappers call"<<std::endl;

}

Creator<Contact::Factory, FrictionContact<BSplineModel<2> , PointModel> > MultiAdaptiveBSplinePointFrictionContactClass("FrictionContact", true);
Creator<Contact::Factory, FrictionContact<BSplineModel<2> , SphereModel> > MultiAdaptiveBSplineSphereFrictionContactClass("FrictionContact", true);


} // namespace collision

} // namespace component

} // namespace sofa

#endif  /* SOFA_DEV */
