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
//
// C++ Implementation : SutureController
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_SutureController_INL
#define SOFA_COMPONENT_CONTROLLER_SutureController_INL

#include "SutureController.h"
#include "WireBeamInterpolation.h"

//#define DEBUG

namespace sofa
{

namespace component
{

namespace controller
{

template <class DataTypes>
SutureController<DataTypes>::SutureController(fem::WireBeamInterpolation<DataTypes>* _adaptiveinterpolation)
: startingPos(initData(&startingPos,Coord(),"startingPos","starting pos for inserting the instrument"))
, threshold(initData(&threshold, (Real)0.000001, "threshold", "threshold for controller precision which is homogeneous to the unit of length"))
, maxBendingAngle(initData(&maxBendingAngle, (Real)0.1, "maxBendingAngle", "max bending criterion (in rad) for one beam"))
, useDummyController(initData(&useDummyController, false, "useDummyController"," use a very simple controller of adaptativity (use for debug)" ))
, fixRigidTransforms(initData(&fixRigidTransforms, false, "fixRigidTransforms", "fix the sampling and transformations of rigid segments"))
, m_rigidCurvAbs(initData(&m_rigidCurvAbs, "rigidCurvAbs", "pairs of curv abs for beams we want to rigidify"))
, m_adaptiveinterpolation(initLink("interpolation", "Path to the Interpolation component on scene"), _adaptiveinterpolation)
, m_nodeCurvAbs(initData(&m_nodeCurvAbs, "nodeCurvAbs", ""))
, m_curvatureList(initData(&m_curvatureList, "curvatureList", "List of the beams curvature (abscissa - curvature)"))
, m_controlPoints(initData(&m_controlPoints, "controlPoints", "List of the spline control points positions"))
, m_topology(0)
{
}

template <class DataTypes>
void SutureController<DataTypes>::init()
{
	this->f_listening.setValue(true);
	Inherit::init();

	if(!m_adaptiveinterpolation)
           {
	      core::objectmodel::BaseContext *c=this->getContext();
	      m_adaptiveinterpolation.set(c->get<WInterpolation>(core::objectmodel::BaseContext::Local));
	   }

	if(!m_adaptiveinterpolation)
		serr << "No Beam Interpolation found, the component can not work!" << sendl;

	m_adaptiveinterpolation->setControlled(true);

	xAbs_collisionPoints_buf.clear();

	/// @TODO : verifier que la topologie est celle d'un wire
	this->getContext()->get(m_topology);

	if (!wireIsAlreadyInitialized())
		initWireModel();
}


template <class DataTypes>
void SutureController<DataTypes>::reinit()
{
    this->getMechanicalState()->cleanup();
    init();
    applyController();
}


template <class DataTypes>
void SutureController<DataTypes>::initWireModel()
{
	std::cout << "SutureController initWireModel\n";

	m_topology->clear();
	m_topology->cleanup();
	
	// on initialise le "wire" en prenant la position de départ + la forme au repos + la discretisation proposée...
	Transform global_T_init;
	const Coord startPos = startingPos.getValue();
	global_T_init.setOrigin(startPos.getCenter());
	global_T_init.setOrientation(startPos.getOrientation());


	helper::vector< Real > xP_noticeable;
	helper::vector< int > nbP_density;
	m_adaptiveinterpolation->getSamplingParameters(xP_noticeable, nbP_density);

	// computation of the number of node on the structure:
	unsigned int numNodes = 1;

	for (unsigned int i=0; i<nbP_density.size(); i++)
	{
		numNodes += nbP_density[i]; // numBeams between each noticeable point
	}

	sofa::helper::vector<Real> &nodeCurvAbs = *m_nodeCurvAbs.beginEdit();

	// Initial position of the nodes:
	nodeCurvAbs.clear();

	Real x_curv = 0.0;

	Data<VecCoord>* datax = this->getMechanicalState()->write(sofa::core::VecCoordId::position());
	Data<VecDeriv>* datav = this->getMechanicalState()->write(sofa::core::VecDerivId::velocity());
	VecCoord& x = *datax->beginEdit();
	VecDeriv& v = *datav->beginEdit();

	this->getMechanicalState()->resize(numNodes);

	x.resize(numNodes);
	v.resize(numNodes);

	x[0].getCenter()= global_T_init.getOrigin();
	x[0].getOrientation()= global_T_init.getOrientation();

	unsigned int id_node=0;
	m_topology->addPoint(x[0][0], x[0][1], x[0][2]);
	nodeCurvAbs.push_back(0.0);

#ifdef DEBUG
	std::cout<< "Init the positions of the beam "<<std::endl;
#endif
	m_adaptiveinterpolation->clear();

	for (unsigned int i=0; i<nbP_density.size(); i++)
	{
		if (nbP_density[i]<=0)
			continue;

		Real length = xP_noticeable[i+1] - xP_noticeable[i];
		Real dx = length/nbP_density[i];

		for(int p=0; p<nbP_density[i]; p++)
		{
			x_curv+=dx;
			id_node++;

			Transform global_H_localP;
			m_adaptiveinterpolation->getRestTransformOnX(global_H_localP, x_curv);
#ifdef DEBUG
			std::cout<<" getRestTransformOnX ="<<global_H_localP<<" at "<< x_curv <<std::endl;
#endif


			Transform global_H_P = global_T_init * global_H_localP;

			// Todo => Put in Mstate !!
			x[id_node].getCenter()= global_H_P.getOrigin();
			x[id_node].getOrientation()= global_H_P.getOrientation();


			// modif the topology

			m_topology->addEdge( (int)(id_node-1), (int)(id_node) );
			m_topology->addPoint( x[id_node][0], x[id_node][1], x[id_node][2] );

			nodeCurvAbs.push_back(x_curv);


			// add the beam to the m_adaptiveinterpolation
			m_adaptiveinterpolation->addBeam(id_node-1,dx,x_curv-dx,x_curv,0.0  );

#ifdef DEBUG
			std::cout<<" Pos at x_curv ="<<x_curv<<" : "<<x[id_node]<<std::endl;
#endif
		}

	}

	datax->endEdit();
	datav->endEdit();
	m_nodeCurvAbs.endEdit();
}


template <class DataTypes>
bool SutureController<DataTypes>::wireIsAlreadyInitialized()
{
	unsigned int numDofs = 0;

	if (this->getMechanicalState() != NULL)
	{
		numDofs = this->getMechanicalState()->getX()->size();
		if (numDofs == 0)
			return false;
	}
	else
	{
		std::cerr << "SutureController should have a MechanicalState in its context\n";
		return false;
	}

	// When there are rigid segments, # of dofs is different than # of edges and beams
	unsigned int numRigidPts = 0;
	helper::ReadAccessor< Data< helper::set< Real > > > rigidCurvAbs = m_rigidCurvAbs;
	int nbRigidAbs = rigidCurvAbs->size();
	if (nbRigidAbs>0 && (nbRigidAbs%2)==0)
	{
		const helper::vector<Real>& curvAbs = m_nodeCurvAbs.getValue();
		RealConstIterator it;
		unsigned int i = 0, nbCurvAbs = curvAbs.size();
		for(it=rigidCurvAbs->begin(); it!=rigidCurvAbs->end();)
		{
			Real start, end;
			start = *it++;
			end = *it++;

			// Look for the start of the rigid segment
			while(i<nbCurvAbs && curvAbs[i] < start)
				i++;
			
			// Count the # of points in this rigid segment
			unsigned int tmpNumRigidPts = 0;
			while(i<nbCurvAbs && curvAbs[i] < end)
			{
				i++;
				tmpNumRigidPts++;
			}

			if(!tmpNumRigidPts)
				tmpNumRigidPts = 1;	// At least one beam for this rigid segment
			numRigidPts += tmpNumRigidPts;
		}
	}

	if (m_topology != NULL)
	{
		if ((unsigned int)m_topology->getNbPoints() != numDofs || (unsigned int)m_topology->getNbEdges() != (numDofs+numRigidPts-1))
			return false;
	}
	else
	{
		std::cerr << "SutureController should have a topology container in its context\n";
		return false;
	}

	if (m_adaptiveinterpolation != NULL)
	{
		if (m_adaptiveinterpolation->getNumBeams() != (numDofs + numRigidPts - 1))
			return false;
	}
	else
	{
		std::cerr << "SutureController should have a WireBeamInterpolation in its context\n";
		return false;
	}

	if (m_nodeCurvAbs.getValue().size() != (numDofs + numRigidPts))
		return false;

	return true;
}

///////////// TODO: Faire une fonction qui recalcule entièrement la topologie (pour remplacer addNodesAndEdge et removeNodesAndEdge


template <class DataTypes>
void SutureController<DataTypes>::recreateTopology()
{
    /* A chaque pas de temps, ON REMET A JOUR LA TOPO ENTIEREMENT */
    m_topology->cleanup();
    m_topology->clear();
}



template <class DataTypes>
void SutureController<DataTypes>::addNodesAndEdge(unsigned int num, Real &xend)
{
    unsigned int numNodes = this->getMechanicalState()->getSize();
    this->getMechanicalState()->resize( numNodes + num);
    unsigned int numBeams = m_adaptiveinterpolation->getNumBeams();
#ifdef DEBUG
    std::cout<<" ####  addNodesAndEdge -- Size Mstate ="<<this->getMechanicalState()->getSize()<<std::endl;
#endif
    for (unsigned int i=0; i<num; i++)
    {
        m_topology->addPoint( 0.0, 0.0, 0.0 );
        m_topology->addEdge( (int)(numNodes-1+i), (int)(numNodes+i) );
        m_adaptiveinterpolation->addBeam(numBeams-1+i,xend,0.0,xend,0.0  ); // the parameters will be coerrected in applyNewSampling
    }
}


template <class DataTypes>
void SutureController<DataTypes>::removeNodesAndEdge(unsigned int num)
{
    unsigned int numNodes = this->getMechanicalState()->getSize();
    if (numNodes - num <= 0)
    {
        serr<<"  in removeNodesAndEdge : no more nodes !!"<<sendl;
        return;
    }

    unsigned int numBeams = m_adaptiveinterpolation->getNumBeams();


    /* ON REMET A JOUR LA TOPO ENTIEREMENT */
    m_topology->cleanup();
    m_topology->clear();

    m_topology->addPoint(0.0,0.0,0.0);
    for (int i=0; i<(int) (numBeams-num); i++)
    {
        m_topology->addPoint((double)(i+1),0.0,0.0);
        m_topology->addEdge(i,i+1);
    }

    this->getMechanicalState()->resize( numNodes - num);
    std::cout<<" #### removeNodesAndEdge -- Size Mstate ="<<this->getMechanicalState()->getSize()<<std::endl;
}


template <class DataTypes>
void SutureController<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
	applyController();
}


template <class DataTypes>
void SutureController<DataTypes>::dummyController(sofa::helper::vector<Real> &newCurvAbs)
{
    static int compteur = 0;
    Real length = m_nodeCurvAbs.getValue()[m_nodeCurvAbs.getValue().size() - 1];
    newCurvAbs.clear();

    for (unsigned int i=0; i<=10; i++)
    {
        Real xtest = length * (Real)i / 10.0;
        newCurvAbs.push_back( xtest );
    }

    Real decalage;
    if (compteur<50)
        decalage=-compteur*0.1;
    else
        decalage=(compteur -100)*0.1;


    for (unsigned int i=1; i<newCurvAbs.size()-1; i++)
    {
        newCurvAbs[i] = newCurvAbs[i] + decalage;
    }


    // between compteur 10 and 90 A node is added !!
    if( compteur > 10 && compteur < 90)
    {
        Real totalLength = newCurvAbs[newCurvAbs.size()-1];
        newCurvAbs[newCurvAbs.size()-1] = totalLength + decalage;
        newCurvAbs.push_back(totalLength );
    }


    compteur++;
    if (compteur==100)
        compteur=0;


    ////////// rigidify a part of the beam //////
    this->rigidCurveSegments.clear();

    std::pair<Real, Real> rigidSegment;

    rigidSegment.first=3.2*length/10.0;
    rigidSegment.second = 5.5*length/10.0;
    this->rigidCurveSegments.push_back(rigidSegment);



    rigidSegment.first = 7.5*length/10.0;
    rigidSegment.second = 10*length/10.0;
    this->rigidCurveSegments.push_back(rigidSegment);



    addRigidCurvAbs(newCurvAbs, 0.0001);



    Real totalLength = newCurvAbs[newCurvAbs.size()-1];

    listOfImposedNodesOnXcurv.push_back(0.23*totalLength);

    listOfImposedNodesOnXcurv.push_back(0.57*totalLength);



   // newCurvAbs.push_back(length);
#ifdef DEBUG
    std::cout<<"newCurvAbs  = "<<newCurvAbs<<std::endl;
#endif
}


template <class DataTypes>
void SutureController<DataTypes>::addRigidCurvAbs(sofa::helper::vector<Real> &newCurvAbs, const Real &tol)
{


    sofa::helper::vector<Real> newCurvAbsBuf=newCurvAbs;

    newCurvAbs.clear();

    unsigned int iterator=1;
    newCurvAbs.push_back(newCurvAbsBuf[0]);

    for (unsigned int i=0; i<rigidCurveSegments.size(); i++)
    {

        while(rigidCurveSegments[i].first > newCurvAbsBuf[iterator])
        {
            newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
            iterator++;
        }

        // newCurvAbsBuf[iterator-1]  <  rigidCurveSegments[i].first <  newCurvAbsBuf[iterator]
        //=>  if rigidCurveSegments[i].first ~ newCurvAbsBuf[iterator-1]
        // or if rigidCurveSegments[i].first ~ newCurvAbsBuf[iterator]  => do not add the point
        if (rigidCurveSegments[i].first - newCurvAbsBuf[iterator-1] > tol && newCurvAbsBuf[iterator] - rigidCurveSegments[i].first > tol)
        {
            newCurvAbs.push_back(rigidCurveSegments[i].first);
        }



        while(rigidCurveSegments[i].second > newCurvAbsBuf[iterator])
        {
            newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
            iterator++;
        }

        // newCurvAbsBuf[iterator-1]  <  rigidCurveSegments[i].second <  newCurvAbsBuf[iterator]
        //=>  if rigidCurveSegments[i].second ~ newCurvAbsBuf[iterator-1]
        // or if rigidCurveSegments[i].second ~ newCurvAbsBuf[iterator]  => do not add the point
        if (rigidCurveSegments[i].second - newCurvAbsBuf[iterator-1] > tol && newCurvAbsBuf[iterator] - rigidCurveSegments[i].second > tol)
        {
            newCurvAbs.push_back(rigidCurveSegments[i].second);
        }

    }


    while(iterator< newCurvAbsBuf.size())
    {
        newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
        iterator++;
    }
#ifdef DEBUG
    std::cout<<"newCurvAbsBuf = "<<newCurvAbsBuf<<std::endl;
    std::cout<<"newCurvAbs    = "<<newCurvAbs<<std::endl;
#endif


}



template <class DataTypes>
void SutureController<DataTypes>::addImposedCurvAbs(sofa::helper::vector<Real> &newCurvAbs, const Real &tol)
{

#ifdef DEBUG
    std::cout<<" --------- addImposedCurvAbs  called with tolerance= "<< tol<<"   ----------"<<std::endl;
    std::cout<<" newCurvAbs = "<<newCurvAbs<<std::endl;
#endif

    listOfImposedNodesOnXcurv.sort();
    listOfImposedNodesOnXcurv.unique();

    ListRealIterator it_xcurv_imposed;

    sofa::helper::vector<Real> newCurvAbsBuf=newCurvAbs;
    newCurvAbs.clear();
    unsigned int iterator=1;
    newCurvAbs.push_back(newCurvAbsBuf[0]);

#ifdef DEBUG
    std::cout<<" listOfImposedNodesOnXcurv = ";
#endif
    for (it_xcurv_imposed=listOfImposedNodesOnXcurv.begin(); it_xcurv_imposed!=listOfImposedNodesOnXcurv.end(); it_xcurv_imposed++)
    {

#ifdef DEBUG
        std::cout<<" "<<(*it_xcurv_imposed);
#endif
        while( (*it_xcurv_imposed) > newCurvAbsBuf[iterator])   // newCurvAbsBuf=  [ 0 2 5 7]    // xcurvImposed= [ 2.1  3  4.9 ]
        {
            newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
            iterator++;
        }

        if ( newCurvAbsBuf[iterator] - (*it_xcurv_imposed)  > tol  &&  (*it_xcurv_imposed) - newCurvAbsBuf[iterator-1] > tol)
            // cas 2 par exemple =>  xcurvImposed=3  newCurvAbs[iterator]= 5
        {
            newCurvAbs.push_back( (*it_xcurv_imposed) );            // newCurvAbsBuf=  [ 0 2 3 ... ]
        }
        else if( (*it_xcurv_imposed) - newCurvAbsBuf[iterator-1] < tol ) //  cas 1 par exemple => xcurvImposed=2.1  newCurvAbs
        {
            newCurvAbs[iterator-1] = (*it_xcurv_imposed);   // newCurvAbsBuf= [0 2.1 ...]
        }
        else
        {
            newCurvAbs.push_back( (*it_xcurv_imposed) );
            iterator++;
        }


    }

    while(iterator < newCurvAbsBuf.size())
    {
        newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
        iterator++;
    }

#ifdef DEBUG
    std::cout<<" "<<std::endl;

    std::cout<<" result : newCurvAbs  ="<<newCurvAbs<<std::endl;
#endif




}




template <class DataTypes>
void SutureController<DataTypes>::applyController()
{
//    std::cout<<"applyController: numBeams =" <<m_adaptiveinterpolation->getNumBeams()<<std::endl;
    Data<VecCoord>* datax = this->getMechanicalState()->write(sofa::core::VecCoordId::position());
    Data<VecDeriv>* datav = this->getMechanicalState()->write(sofa::core::VecDerivId::velocity());
    VecCoord& x = *datax->beginEdit();
    VecDeriv& v = *datav->beginEdit();

    sofa::helper::vector<Real> newCurvAbs;

	this->storeRigidSegmentsTransformations();

    if (useDummyController.getValue())
        this->dummyController(newCurvAbs);
    else
        this->computeSampling(newCurvAbs, x);

    this->addImposedCurvAbs(newCurvAbs, 0.0001);

	this->verifyRigidSegmentsSampling(newCurvAbs);

#ifdef DEBUG
    std::cout<<" newCurvAbs = "<<newCurvAbs<<"  - nodeCurvAbs = "<<m_nodeCurvAbs.getValue()<<std::endl;

    for (unsigned int i=0; i<rigidCurveSegments.size(); i++)
    {
        std::cout<<" rigidCurveSegments["<<i<<"] "<<rigidCurveSegments[0].first<<" - "<<rigidCurveSegments[0].second<<std::endl;
    }

    std::cout<<" Before... Beam  lengths = ";
    for (unsigned int b=0; b<m_adaptiveinterpolation->getNumBeams(); b++)
        std::cout<<" "<<m_adaptiveinterpolation->getLength(b);
    std::cout<<" "<<std::endl;

#endif

    this->applyNewSampling(newCurvAbs, m_nodeCurvAbs.getValue(), x, v);

 #ifdef DEBUG
    std::cout<<"After.... Beam  lengths = ";

    for (unsigned int b=0; b<m_adaptiveinterpolation->getNumBeams(); b++)
        std::cout<<" "<<m_adaptiveinterpolation->getLength(b);
    std::cout<<" "<<std::endl;
#endif

	this->verifyRigidSegmentsTransformations();
	prevRigidCurvSegments = rigidCurveSegments;

	sofa::helper::vector<Real> &nodeCurvAbs = *m_nodeCurvAbs.beginEdit();
	nodeCurvAbs.assign(newCurvAbs.begin(), newCurvAbs.end());
	m_nodeCurvAbs.endEdit();

	VecCoord& ctrlPts = *m_controlPoints.beginEdit();
	ctrlPts.clear();

	unsigned int numBeams = m_adaptiveinterpolation->getNumBeams();
	Transform global_H0_local,  global_H1_local;
	for (unsigned int b=0; b<numBeams; b++)
	{
		m_adaptiveinterpolation->computeTransform2(b, global_H0_local, global_H1_local, x);
		Coord pt;
		pt.getCenter() = global_H0_local.getOrigin();
		pt.getOrientation() = global_H0_local.getOrientation();
		ctrlPts.push_back(pt);
	}
	Coord pt;
	pt.getCenter() = global_H1_local.getOrigin();
	pt.getOrientation() = global_H1_local.getOrientation();
	ctrlPts.push_back(pt);
	
	m_controlPoints.endEdit();
	
	datax->endEdit();
    datav->endEdit();
}


//////*************** PRIVATE FUNCTIONS ****************//

/// this function calls ComputeTotalBendingRotationAngle on the beams between xmin and xmax
template <class DataTypes>
typename SutureController<DataTypes>::Real SutureController<DataTypes>::computeBendingAngle(const Real& xmin, const Real& xmax, const Real& dx_comput, const VecCoord& Pos)
{

    //test for verification
    if (xmax < xmin || dx_comput==0.0){
        serr<<"wrong parameters in computeBendingAngle function"<<sendl;
        return 0.0;
    }

    if (xmax > m_adaptiveinterpolation->getRestTotalLength()){
        serr<<" in computeBendingAngle : max > getRestTotalLength"<<sendl;
        return 0.0;
    }


    unsigned int idBeamMin, idBeamMax;
    Real baryCoordMin, baryCoordMax;
    Transform Tnode0, Tnode1;

    m_adaptiveinterpolation->getBeamAtCurvAbs(xmin, idBeamMin, baryCoordMin);
    m_adaptiveinterpolation->getBeamAtCurvAbs(xmax, idBeamMax, baryCoordMax);
    m_adaptiveinterpolation->computeTransform2(idBeamMin,Tnode0,Tnode1, Pos);

    if (idBeamMin==idBeamMax)
    {
         //std::cout<<" idBeamMin==idBeamMax, angle ="<<angle<<std::endl;
        return m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(dx_comput, Tnode0, Tnode1, m_adaptiveinterpolation->getLength(idBeamMin), baryCoordMin, baryCoordMax);
    }

    if (idBeamMin>idBeamMax)
    {
        serr<<"ERROR in computeBendingAngle"<<sendl;
        return 0.0;
    }

    //////////


    // compute the angle for the first beam

    Real angle = m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(dx_comput, Tnode0, Tnode1, m_adaptiveinterpolation->getLength(idBeamMin), baryCoordMin, 1.0);

    unsigned int b=idBeamMin+1;
    while(b<idBeamMax)
    {
        m_adaptiveinterpolation->computeTransform2(b,Tnode0,Tnode1,Pos);
        angle += m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(dx_comput, Tnode0, Tnode1, m_adaptiveinterpolation->getLength(b), 0.0, 1.0);
        b++;
    }

    m_adaptiveinterpolation->computeTransform2(idBeamMax,Tnode0,Tnode1,Pos);
    angle += m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(dx_comput, Tnode0, Tnode1, m_adaptiveinterpolation->getLength(idBeamMax), 0.0, baryCoordMax);

	return angle;
}


/// this function computes the tangent value on a series of discrete points (store also the curv_abs of these discrete points)
template <class DataTypes>
void SutureController<DataTypes>::computeTangentOnDiscretePoints(sofa::helper::vector<Vec3> TangTable, sofa::helper::vector<Real> xTable,  unsigned int numDiscretePoints, const VecCoord& Pos)
{

    TangTable.clear();
    Transform Tnode0, Tnode1;
    Real baryCoord;
    unsigned int beam;

    Real dx = m_adaptiveinterpolation->getRestTotalLength()/(numDiscretePoints-1);
    Real x=dx;


    // compute intial tang for the beginning of the wire:
    m_adaptiveinterpolation->computeTransform2(0,Tnode0,Tnode1, Pos);
    Vec3 t = Tnode0.getOrientation().rotate(Vec3(1.0,0.0,0.0));
    TangTable.push_back(t);
    xTable.push_back(0.0);

    for (unsigned int p=0; p<(numDiscretePoints-1) ; p++)
    {
        m_adaptiveinterpolation->getBeamAtCurvAbs(x, beam, baryCoord);
        m_adaptiveinterpolation->computeTransform2(beam,Tnode0,Tnode1, Pos);
        m_adaptiveinterpolation->getTangent(t, baryCoord, Tnode0,Tnode1,m_adaptiveinterpolation->getLength(beam) );

        TangTable.push_back(t);
        xTable.push_back(x);

        x+=dx;

    }
}



template <class DataTypes>
void SutureController<DataTypes>::detectRigidBeams(const sofa::helper::vector<Real> &newCurvAbs)
{

    unsigned int seg=0;

    rigidBeamList.clear();

    for (unsigned int i=1; i<newCurvAbs.size(); i++)
    {
        ////////////////////////// Rigidification ////////////////////////
        //(1) look if there is a rigidCurveSegment...
        if(seg < rigidCurveSegments.size() )
        {
            // (2) look if [newCurvAbs[i-1] newCurvAbs[i] are contained in the next rigidified segment....
            if(newCurvAbs[i-1]+threshold.getValue() > rigidCurveSegments[seg].first && newCurvAbs[i]-threshold.getValue() < rigidCurveSegments[seg].second )
                rigidBeamList.push_back(true);

            else
                rigidBeamList.push_back(false);


            if( newCurvAbs[i]+threshold.getValue() > rigidCurveSegments[seg].second)
                seg++;

        }
        else
            rigidBeamList.push_back(false);

    }
#ifdef DEBUG
    std::cout<<" detectRigidBeams result : "<<rigidBeamList<<std::endl;
#endif
}



// When a new sampling is defined in "newCurvAbs", the position and the velocity needs to be "re-interpolated"

template <class DataTypes>
void SutureController<DataTypes>::applyNewSampling(const sofa::helper::vector<Real> &newCurvAbs, const sofa::helper::vector<Real> &oldCurvAbs, VecCoord &x, VecDeriv &v)
{

    VecCoord x_buf=x;
    VecDeriv v_buf=v;



    x.clear();
    v.clear();

#ifdef DEBUG
    std::cout<<" Begin applyNewSampling : newCurvAbs="<<newCurvAbs<<std::endl;
#endif
    detectRigidBeams(newCurvAbs);



    ////////////////////// interpolation of the position and velocities //////////

    Transform global_H_interpol;
    Deriv v_interpol;
    unsigned int j=0;
    Vec3 null(0,0,0);

    vec_global_H_node.clear();
    vec_global_Vel_node.clear();

    m_adaptiveinterpolation->InterpolateTransformAndVelUsingSpline(0,0.0,null,x_buf, v_buf, global_H_interpol, v_interpol);
    vec_global_H_node.push_back(global_H_interpol);
    vec_global_Vel_node.push_back(v_interpol);


    for (unsigned int i=1; i<newCurvAbs.size()-1; i++)
    {

        while(newCurvAbs[i]>oldCurvAbs[j])
        {
                j++;
                if (j>=oldCurvAbs.size()) // DEBUG //
                {
                        std::cerr<<"**************** > WARNING j ="<<j<<">=oldCurvAbs.size()"<<std::endl;
                        return;
                }
        }
        Real L = m_adaptiveinterpolation->getLength(j-1);
        Real L0 = newCurvAbs[i] - oldCurvAbs[j-1];
        Real ratio=L0/L;

        m_adaptiveinterpolation->InterpolateTransformAndVelUsingSpline(j-1,ratio,null,x_buf, v_buf, global_H_interpol, v_interpol);
        vec_global_H_node.push_back(global_H_interpol);
        vec_global_Vel_node.push_back(v_interpol);

    }
    m_adaptiveinterpolation->InterpolateTransformAndVelUsingSpline(oldCurvAbs.size()-2,1.0,null,x_buf, v_buf, global_H_interpol, v_interpol);
    vec_global_H_node.push_back(global_H_interpol);
    vec_global_Vel_node.push_back(v_interpol);


#ifdef DEBUG
    std::cout<<"vec_global_H_node = "<<vec_global_H_node<<std::endl;
    std::cout<<"vec_global_Vel_node = "<<vec_global_Vel_node<<std::endl;
#endif


    ////////////////////// compute a gravity center for the rigid segments //////////
    vec_global_H_gravityCenter.clear();
    sofa::helper::vector<Deriv> vec_Vel_gravityCenter;


    Transform global_H_gravityC;
    Deriv vel_gravityC;
    bool rigidification;
    rigidification=false;
    unsigned int Rseg=0;
    Real length_of_rigidSegment ;


    for (unsigned int s=0; s<rigidBeamList.size();s++)
    {
        if(rigidBeamList[s]) // the beam is rigidified
        {
            if (!rigidification) //begining of the rigidification
            {

                length_of_rigidSegment = rigidCurveSegments[Rseg].second  - rigidCurveSegments[Rseg].first;
                Rseg++;
                rigidification=true;
                global_H_gravityC.clear();

            }

            Real length_of_beam = newCurvAbs[s+1]-newCurvAbs[s];
            Real weight = (length_of_beam/(length_of_rigidSegment*2));

            global_H_gravityC.setOrigin( vec_global_H_node[s+1].getOrigin()*weight
                                         + vec_global_H_node[s].getOrigin()*weight
                                         + global_H_gravityC.getOrigin() );


            vel_gravityC.getVCenter() += vec_global_Vel_node[s+1].getVCenter()*weight+ vec_global_Vel_node[s].getVCenter()*weight;
            vel_gravityC.getVOrientation() += vec_global_Vel_node[s+1].getVOrientation()*weight+ vec_global_Vel_node[s].getVOrientation()*weight;

        }
        else
        {
            if (rigidification)
            {

                Vec3 pos_G = global_H_gravityC.getOrigin();
                global_H_gravityC.setOrientation(vec_global_H_node[s].getOrientation() );
                global_H_gravityC.setOrigin(pos_G);

 #ifdef DEBUG
                std::cout<<" global_H_gravityC.origin = "<<global_H_gravityC.getOrigin()<<"  global_H_gravityC ="<<global_H_gravityC <<std::endl;


                /////////////////// Computation of the velocities (can be put in an other function) /////////////////

                std::cout<<"rigidification using vec_global_Vel_node["<<s+1<<"] = "<<vec_global_Vel_node[s+1]<<std::endl;
#endif
                SpatialVector VelNode_in_global, VelNode_in_Node, VelGravityC_inGravityC, VelGravity_inGlobal;
                VelNode_in_global.setAngularVelocity( vec_global_Vel_node[s].getVOrientation());
                VelNode_in_global.setLinearVelocity(  vec_global_Vel_node[s].getVCenter());


                // projection of the velocity of the node in the frame of the node
                VelNode_in_Node.setLinearVelocity( vec_global_H_node[s].backProjectVector(VelNode_in_global.getLinearVelocity()) );
                VelNode_in_Node.setAngularVelocity( vec_global_H_node[s].backProjectVector(VelNode_in_global.getAngularVelocity() ) );



                // TRANSPORT of the velocity of the node to the gravity center (rigid link)
                Transform gravityC_H_Node = global_H_gravityC.inversed()*vec_global_H_node[s];
                VelGravityC_inGravityC = gravityC_H_Node*VelNode_in_Node;


                // projection of the velocity of the gravity center in the global frame
                VelGravity_inGlobal.setLinearVelocity(global_H_gravityC.projectVector( VelGravityC_inGravityC.getLinearVelocity() ) );
                VelGravity_inGlobal.setAngularVelocity(  global_H_gravityC.projectVector( VelGravityC_inGravityC.getAngularVelocity() ) );

                vel_gravityC.getVCenter()       = VelGravity_inGlobal.getLinearVelocity();
                vel_gravityC.getVOrientation()  = VelGravity_inGlobal.getAngularVelocity();


                /////////////////////////////////////  ///////////////////////////////////  //////////////////////////////////



                vec_global_H_gravityCenter.push_back(global_H_gravityC);
                vec_Vel_gravityCenter.push_back(vel_gravityC);
                rigidification=false; // end of the rigidification

                global_H_gravityC.clear();
            }

        }

    }
    if (rigidification) // rigidification at the end tip
    {

        unsigned s=rigidBeamList.size()-1;
#ifdef DEBUG
        std::cout<<"!!!! Rigidification at the end tip !!!! "<<std::endl;
#endif
        Vec3 pos_G = global_H_gravityC.getOrigin();
        global_H_gravityC.setOrientation(vec_global_H_node[s+1].getOrientation() );
        global_H_gravityC.setOrigin(pos_G);


        /////////////////// Computation of the velocities (can be put in an other function) /////////////////


        SpatialVector VelNode_in_global, VelNode_in_Node, VelGravityC_inGravityC, VelGravity_inGlobal;
        VelNode_in_global.setAngularVelocity( vec_global_Vel_node[s].getVOrientation());
        VelNode_in_global.setLinearVelocity(  vec_global_Vel_node[s].getVCenter());

        // projection of the velocity of the node in the frame of the node
        VelNode_in_Node.setLinearVelocity( vec_global_H_node[s].backProjectVector(VelNode_in_global.getLinearVelocity()) );
        VelNode_in_Node.setAngularVelocity( vec_global_H_node[s].backProjectVector(VelNode_in_global.getAngularVelocity() ) );

        // TRANSPORT of the velocity of the node to the gravity center (rigid link)
        Transform gravityC_H_Node = global_H_gravityC.inversed()*vec_global_H_node[s];
        VelGravityC_inGravityC = gravityC_H_Node*VelNode_in_Node;

        // projection of the velocity of the gravity center in the global frame
        VelGravity_inGlobal.setLinearVelocity(global_H_gravityC.projectVector( VelGravityC_inGravityC.getLinearVelocity() ) );
        VelGravity_inGlobal.setAngularVelocity(  global_H_gravityC.projectVector( VelGravityC_inGravityC.getAngularVelocity() ) );

        vel_gravityC.getVCenter()       = VelGravity_inGlobal.getLinearVelocity();
        vel_gravityC.getVOrientation()  = VelGravity_inGlobal.getAngularVelocity();


        /////////////////////////////////////  ///////////////////////////////////  //////////////////////////////////

        vec_global_H_gravityCenter.push_back(global_H_gravityC);
        vec_Vel_gravityCenter.push_back(vel_gravityC);
        rigidification=false; // end of the rigidification
    }




    //////////////// Set the beam and the topology /////////////////

    recreateTopology();
    unsigned int numNodes=0;
    m_adaptiveinterpolation->clear();

    numNodes++;
    Coord xDof;
    xDof.getCenter()     = vec_global_H_node[0].getOrigin();
    xDof.getOrientation()= vec_global_H_node[0].getOrientation();
    x.push_back(xDof);
    v.push_back(vec_global_Vel_node[0]);
    m_topology->addPoint( xDof[0], xDof[1], xDof[2] );
    Real L;

    Rseg=0;
    rigidification=false;
    for (unsigned int s=0; s<rigidBeamList.size();s++)
    {
        if(rigidBeamList[s]) // the beam is rigidified
        {

            if (!rigidification) //begining of the rigidification
            {
                // the last element of vector x is replaced by the gravity center of the rigid zone
                x.pop_back();
                xDof.getCenter()     = vec_global_H_gravityCenter[Rseg].getOrigin();
                xDof.getOrientation()= vec_global_H_gravityCenter[Rseg].getOrientation();
                x.push_back(xDof);
                v.pop_back();
                v.push_back( vec_Vel_gravityCenter[Rseg] );
                Rseg++;
                rigidification=true;

                // add a transformation between the node 1 and the gravity center on previous beam
                Transform GravityCenter_H_node1 = vec_global_H_gravityCenter[Rseg-1].inversed()*vec_global_H_node[s];
                m_adaptiveinterpolation->setTransformBetweenDofAndNode(s-1,GravityCenter_H_node1,1);

            }

            Transform GravityCenter_H_interpol0 = vec_global_H_gravityCenter[Rseg-1].inversed()*vec_global_H_node[s];
            Transform GravityCenter_H_interpol1 = vec_global_H_gravityCenter[Rseg-1].inversed()*vec_global_H_node[s+1];




             //////// ADD A BEAM THAT IS ON A SEGMENT THAT LINKS THE SAME DOF ///////
            m_topology->addEdge( (int)(numNodes-1), (int)(numNodes-1) );

            L =  newCurvAbs[s+1] - newCurvAbs[s];
            m_adaptiveinterpolation->addBeam(s, L, newCurvAbs[s], newCurvAbs[s+1] ,GravityCenter_H_interpol0, GravityCenter_H_interpol1 );

#ifdef DEBUG
            std::cout<<" add rigid segment on segment "<<s<<" that links the same dof: [ "<<m_topology->getEdge(s)[0]<<" "<<m_topology->getEdge(s)[1]
                    <<"]   x curv = ["<< newCurvAbs[s]<<" " <<newCurvAbs[s+1]<<"]"<<std::endl;
            std::cout<<"GravityCenter_H_interpol0 ="<<GravityCenter_H_interpol0<<"   - GravityCenter_H_interpol1"<<GravityCenter_H_interpol1<<std::endl;
#endif
        }

        else // the beam is deformable
        {


            L =  newCurvAbs[s+1] - newCurvAbs[s];

            // ADD the beam in the topology and in the interpolation

            m_topology->addEdge( (int)(numNodes-1), (int)(numNodes) );
            m_adaptiveinterpolation->addBeam(s, L, newCurvAbs[s], newCurvAbs[s+1] ,0.0 );
#ifdef DEBUG
            std::cout<<" add deformable segment on segment "<<s<<" that links the same dof: [ "<<m_topology->getEdge(s)[0]<<" "<<m_topology->getEdge(s)[1]
                    <<"]   x curv = ["<< newCurvAbs[s]<<" " <<newCurvAbs[s+1]<<"]"<<std::endl;
#endif
            // ADD a DOF for the second node of the beam
            numNodes++;
            xDof.getCenter()     = vec_global_H_node[s+1].getOrigin();
            xDof.getOrientation()= vec_global_H_node[s+1].getOrientation();
            x.push_back(xDof);

            v.push_back(vec_global_Vel_node[s+1]);

            m_topology->addPoint( xDof[0], xDof[1], xDof[2] );


            if (rigidification) // end of the rigidification
            {
                rigidification=false;

                // add a transformation between the node 0 and the gravity center on current beam
                Transform GravityCenter_H_node0 = vec_global_H_gravityCenter[Rseg-1].inversed()*vec_global_H_node[s];
                m_adaptiveinterpolation->setTransformBetweenDofAndNode(s,GravityCenter_H_node0,0);
#ifdef DEBUG
                std::cout<<" Add transformation"<< GravityCenter_H_node0<<"  between the node 0 and the gravity center on beam "<<s<<std::endl;
#endif
            }

        }
    }

#ifdef DEBUG
    std::cout<<"xbuf = "<<x_buf<<std::endl;
    std::cout<<" x   = "<<x<<std::endl;

    std::cout<<"vbuf = "<<v_buf<<std::endl;
    std::cout<<" v   = "<<v<<std::endl;
#endif

    this->getMechanicalState()->resize(x.size());


}

template <class DataTypes>
bool SutureController<DataTypes>::verifyRigidCurveSegmentSort()
{

    if(this->rigidCurveSegments.size()==0)
        return true;

    for (unsigned int seg=0; seg<this->rigidCurveSegments.size()-1; seg++)
    {
        if(rigidCurveSegments[seg].second > rigidCurveSegments[seg+1].first)
            return false;
    }
    return true;

}


template <class DataTypes>
void SutureController<DataTypes>::computeSampling(sofa::helper::vector<Real> &newCurvAbs, VecCoord &x)
{
    helper::vector<Real> xP_noticeable;
    helper::vector<int> nbP_density;

    m_adaptiveinterpolation->getSamplingParameters(xP_noticeable, nbP_density);

    if(xP_noticeable.size()<2){
        serr<<" xP_noticeable_buf.size()= "<<xP_noticeable.size()<<sendl;
        return;
    }

	std::vector<Real> beamsCurvature;
	unsigned int nbBeams = m_adaptiveinterpolation->getNumBeams();
	beamsCurvature.resize(nbBeams);

	helper::WriteAccessor< Data< sofa::helper::vector<Vec2> > > curvatureList = m_curvatureList;
	curvatureList.clear();
	curvatureList.resize(nbBeams);
	// Computing the curvature of each beam (from the previous timestep)
	for(unsigned int b=0; b<nbBeams; ++b)
	{
		Real beamLength = m_adaptiveinterpolation->getLength(b);
		m_adaptiveinterpolation->getAbsCurvXFromBeam(b, curvatureList[b][0]);
		Transform Tnode0, Tnode1;
		m_adaptiveinterpolation->computeTransform2(b, Tnode0, Tnode1, x);

		curvatureList[b][1] = beamsCurvature[b] = m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(beamLength / 5, Tnode0, Tnode1, beamLength, 0.0, 1.0);
	}

	sofa::helper::vector<Real> newCurvAbs_notSecure;
	newCurvAbs_notSecure.clear();
	Real currentCurvAbs = 0.0, currentAngle = 0.0, maxAngle = maxBendingAngle.getValue();
	unsigned int currentBeam = 0;
	for(unsigned int part=0; part<nbP_density.size(); part++)
	{
		Real maxBeamLength = (xP_noticeable[part+1] - xP_noticeable[part]) / nbP_density[part];

		while(currentCurvAbs < xP_noticeable[part+1]-threshold.getValue())
		{
			newCurvAbs_notSecure.push_back(currentCurvAbs);
			Real maxCurvAbs = std::min(currentCurvAbs + maxBeamLength, xP_noticeable[part+1]);

			while(currentBeam < nbBeams)
			{
				Real beamStart, beamEnd, beamLength;
				m_adaptiveinterpolation->getAbsCurvXFromBeam(currentBeam, beamStart, beamEnd);
				beamLength = beamEnd - beamStart;

				Real beamAngle = beamsCurvature[currentBeam];				// Curvature of the whole beam
				Real baryStart = (beamEnd - currentCurvAbs) / beamLength;	// Where we are currently on the beam [0-1]
				Real remainingbeamAngle = beamAngle * baryStart;			// This is the curvature from the current abs to the end of the beam
				if(currentAngle + remainingbeamAngle > maxAngle)			// The new beam will end somewhere on this beam
				{
					Real baryEnd = (maxAngle - currentAngle) / beamAngle;
					currentCurvAbs += beamLength * baryEnd;
					if(currentCurvAbs > maxCurvAbs)
						currentCurvAbs = maxCurvAbs;
					currentAngle = 0.0;
					break;
				}
				else if(beamEnd > maxCurvAbs)								// We got to a limit
				{
					currentCurvAbs = maxCurvAbs;
					currentAngle = 0.0;
					break;
				}
				else														// Continue to next beam
				{
					currentAngle += remainingbeamAngle;
					currentCurvAbs += beamLength * baryStart;
					++currentBeam;
				}
			}
		}

		newCurvAbs_notSecure.push_back(xP_noticeable[part+1]);
	}

/*
    Real sutureLength= m_adaptiveinterpolation->getRestTotalLength();
    Real dx=sutureLength/100.0;	// TODO : pourquoi 100 et pas un Data ou un nombre dépendant de la densité ?
    Real BendingAngle=0.0;
    this->computeBendingAngle(BendingAngle, 0.0, sutureLength, dx, x);

#ifdef DEBUG
    std::cout<<" BendingAngle ="<<BendingAngle<<std::endl;
#endif


    Real L_add=0.0;

    // in density it is defined the number of beams that are basically asked between
    sofa::helper::vector<Real> newCurvAbs_notSecure;
    newCurvAbs_notSecure.clear();
    for (unsigned int part=0; part<nbP_density.size(); part++)
    {
        unsigned int numBeams=nbP_density[part];
        Real L_beam_straight = (xP_noticeable[part+1] - xP_noticeable[part])/numBeams;

        while(L_add< xP_noticeable[part+1]-threshold.getValue())
        {
            newCurvAbs_notSecure.push_back(L_add);

            if (L_add+L_beam_straight < xP_noticeable[part+1])
                BendingAngle = this->computeBendingAngle(L_add, L_add+L_beam_straight, dx, x);
            else
                BendingAngle = this->computeBendingAngle(L_add, xP_noticeable[part+1], dx, x);

            Real L_beam = L_beam_straight;
            if (BendingAngle>maxBendingAngle.getValue())
                L_beam*=maxBendingAngle.getValue()/BendingAngle;
            L_add+=L_beam;
        }
        newCurvAbs_notSecure.push_back(xP_noticeable[part+1]);
    }
	*/
    /*

	// TEMP TEST : remove aligned dofs
    if(newCurvAbs_notSecure.size() >= 3)
	{
		dx = sutureLength / 20;
        Real prevAbs = newCurvAbs_notSecure[0], unused=0.0;
        unsigned int prevBeam = 0;
        m_adaptiveinterpolation->getBeamAtCurvAbs(prevAbs, prevBeam, unused);
		int size = newCurvAbs_notSecure.size();
		for(int i=1; i<size-1; )
		{
			Real curvAbs = newCurvAbs_notSecure[i];
			Real angle = 0;
			this->computeBendingAngle(angle, prevAbs, curvAbs, dx, x);

            unsigned int beamIndex = 0;
            m_adaptiveinterpolation->getBeamAtCurvAbs(curvAbs, beamIndex, unused);

            bool stretched = false;
            for(unsigned int j=beamIndex; j<=beamIndex; j++)
            {
                Vec3 P0,P1,P2,P3;
                Real rest_length = m_adaptiveinterpolation->getLength(j);
                m_adaptiveinterpolation->getSplinePoints(j,x,P0,P1,P2,P3);
                Real length = 0.0;
                m_adaptiveinterpolation->computeActualLength(length, P0,P1,P2,P3);

           //     std::cout<<" test on curvAbs"<<curvAbs<<" - prevBeam="<<prevBeam<<" - beamIndex ="<<beamIndex<<"- length="<<length<<"- rest_length"<<rest_length<<std::endl;

                if( length > 1.01 * rest_length)
                {
                    stretched = true;
                    break;
                }
            }

            if(angle<0.1 && stretched)	// less than 2 degrees is considered plan
			{
				newCurvAbs_notSecure.erase(newCurvAbs_notSecure.begin() + i);
				size--;
			}
			else
			{
				i++;
                prevAbs = curvAbs;

			}
            //prevBeam = beamIndex;
		}
	}
    */

	if(!verifyRigidCurveSegmentSort())
		serr<<" WARNING : rigidCurveSegments are not correctly sorted !"<<sendl;

    this->rigidCurveSegments.clear();
	helper::ReadAccessor< Data< helper::set< Real > > > rigidCurvAbs = m_rigidCurvAbs;
	int nb = rigidCurvAbs->size();
	if(nb>0 && (nb%2)==0)	// Make sure we have pairs of curv abs
	{
        RealConstIterator it;
		for(it=rigidCurvAbs->begin(); it!=rigidCurvAbs->end();)
		{
			Real start, end;
			start = *it++;
			end = *it++;
			this->rigidCurveSegments.push_back(std::make_pair(start, end));
		}
		addRigidCurvAbs(newCurvAbs_notSecure, 0.0001);
	}

    ///// Verify that there is no beams with null length ///
    newCurvAbs.clear();
    newCurvAbs.push_back(newCurvAbs_notSecure[0]);
    for (unsigned int i=1; i<newCurvAbs_notSecure.size(); i++)
    {
        if (newCurvAbs_notSecure[i] > newCurvAbs_notSecure[i-1]+threshold.getValue())
            newCurvAbs.push_back(newCurvAbs_notSecure[i]);
    }

#ifdef DEBUG
    std::cout<<" compute Sampling: newCurvAbs="<<newCurvAbs<<"   xP_noticeable="<<xP_noticeable<<"   nbP_density="<<nbP_density<<std::endl;

#endif


}

template <class DataTypes>
void SutureController<DataTypes>::verifyRigidSegmentsSampling(sofa::helper::vector<Real> &newCurvAbs)
{	// Making sure we keep the same sampling in the rigid segments from one timestep to the next
	if(!fixRigidTransforms.getValue())
		return;

	const sofa::helper::vector<Real> &oldCurvAbs = m_nodeCurvAbs.getValue();
	typename sofa::helper::vector<Real>::iterator newIter, newIter2;
	typename sofa::helper::vector<Real>::const_iterator oldIter, oldIter2;
	newIter = newCurvAbs.begin();
	oldIter = oldCurvAbs.begin();

	typename sofa::helper::vector< std::pair<Real, Real> >::const_iterator rigidIter;
	// For each segment
	for(rigidIter = rigidCurveSegments.begin(); rigidIter != rigidCurveSegments.end(); ++rigidIter)
	{
		if(std::find(prevRigidCurvSegments.begin(), prevRigidCurvSegments.end(), *rigidIter) == prevRigidCurvSegments.end())
			continue;	// If this is a new segment, don't modify it

		Real start = rigidIter->first, end = rigidIter->second;

		// Find indices in the curvAbs lists corresponding to the start of the rigid segment
		newIter = std::upper_bound(newIter, newCurvAbs.end(), start);
		oldIter = std::upper_bound(oldIter, oldCurvAbs.end(), start);

		newIter2 = newIter;
		--newIter;
		oldIter2 = oldIter;
		--oldIter;

		// Find indices in the curvAbs lists corresponding to the end of the rigid segment
		newIter2 = std::upper_bound(newIter, newCurvAbs.end(), end);
		oldIter2 = std::upper_bound(oldIter, oldCurvAbs.end(), end);

		// Removing what was computed in this timestep
		newIter = newCurvAbs.erase(newIter, newIter2);
		// And replacing by the data of the previous timestep
		newCurvAbs.insert(newIter, oldIter, oldIter2);
	}
}

template <class DataTypes>
void SutureController<DataTypes>::storeRigidSegmentsTransformations()
{
	if(!fixRigidTransforms.getValue())
		return;

	prevRigidTransforms.clear();

	typename sofa::helper::vector< std::pair<Real, Real> >::const_iterator rigidIter;
	// For each rigid segment
	for(rigidIter = prevRigidCurvSegments.begin(); rigidIter != prevRigidCurvSegments.end(); ++rigidIter)
	{
		double rigidStart = rigidIter->first, rigidEnd = rigidIter->second;

		// For all curv abs in this segment
		unsigned int beamId, lastBeamId;
		Real bary;
		m_adaptiveinterpolation->getBeamAtCurvAbs(rigidStart, beamId, bary);
		m_adaptiveinterpolation->getBeamAtCurvAbs(rigidEnd, lastBeamId, bary);

		while(beamId < lastBeamId)
		{
			// Save the transformation
			//  (we consider that we don't cut rigid beams so transformations are the same for 2 beams in the same curv abs)
			Real curvAbs0, curvAbs1;
			m_adaptiveinterpolation->getAbsCurvXFromBeam(beamId, curvAbs0, curvAbs1);

			Transform t0, t1;
			m_adaptiveinterpolation->getDOFtoLocalTransform(beamId, t0, t1);
			prevRigidTransforms[curvAbs0] = t0;
			prevRigidTransforms[curvAbs1] = t1;

			++beamId;
		}
	}

/*	int nb = m_adaptiveinterpolation->getNumBeams();
	for(int i=0;i<nb; ++i)
	{
		Real curvAbs0, curvAbs1;
		m_adaptiveinterpolation->getAbsCurvXFromBeam(i, curvAbs0, curvAbs1);

		Transform t0, t1;
		m_adaptiveinterpolation->getDOFtoLocalTransform(i, t0, t1);
		prevRigidTransforms[curvAbs0] = t0;
		prevRigidTransforms[curvAbs1] = t1;
	}	*/
}

template <class DataTypes>
void SutureController<DataTypes>::verifyRigidSegmentsTransformations()
{
	if(!fixRigidTransforms.getValue())
		return;

	typename sofa::helper::vector< std::pair<Real, Real> >::const_iterator rigidIter;
	// For each rigid segment
	for(rigidIter = rigidCurveSegments.begin(); rigidIter != rigidCurveSegments.end(); ++rigidIter)
	{
		if(std::find(prevRigidCurvSegments.begin(), prevRigidCurvSegments.end(), *rigidIter) == prevRigidCurvSegments.end())
			continue;	// If this is a new segment, don't modify it

		Real rigidStart = rigidIter->first, rigidEnd = rigidIter->second;

		// For all curv abs in this segment
		unsigned int beamId, lastBeamId;
		Real bary;
		m_adaptiveinterpolation->getBeamAtCurvAbs(rigidStart, beamId, bary);
		m_adaptiveinterpolation->getBeamAtCurvAbs(rigidEnd, lastBeamId, bary);

		while(beamId < lastBeamId)
		{
			// Load the transformations and replace what was computed this timestep
			Real curvAbs0, curvAbs1;
			m_adaptiveinterpolation->getAbsCurvXFromBeam(beamId, curvAbs0, curvAbs1);

			typename std::map<Real, Transform>::const_iterator iter;
			iter = prevRigidTransforms.find(curvAbs0);
			if(iter != prevRigidTransforms.end())
				m_adaptiveinterpolation->setTransformBetweenDofAndNode(beamId, iter->second, false);

			iter = prevRigidTransforms.find(curvAbs1);
			if(iter != prevRigidTransforms.end())
				m_adaptiveinterpolation->setTransformBetweenDofAndNode(beamId, iter->second, true);

			++beamId;
		}
	}	
	/*
	// TODO : adapt this simpler method to rigid segments that have been modified (the transformation need to change)
	int nb = m_adaptiveinterpolation->getNumBeams();
	for(int i=0;i<nb; ++i)
	{
		Real curvAbs0, curvAbs1;
		m_adaptiveinterpolation->getAbsCurvXFromBeam(i, curvAbs0, curvAbs1);

		std::map<Real, Transform>::const_iterator iter;
		iter = prevRigidTransforms.find(curvAbs0);
		if(iter != prevRigidTransforms.end())
			m_adaptiveinterpolation->setTransformBetweenDofAndNode(i, iter->second, false);

		iter = prevRigidTransforms.find(curvAbs1);
		if(iter != prevRigidTransforms.end())
			m_adaptiveinterpolation->setTransformBetweenDofAndNode(i, iter->second, true);
	}*/
}

template <class DataTypes>
void SutureController<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowBehaviorModels()) return;

    if (rigidCurveSegments.size() != vec_global_H_gravityCenter.size())
    {
        serr<<"in draw function rigidCurveSegments.size() ="<< rigidCurveSegments.size() <<" != vec_global_H_gravityCenter.size() = "<<vec_global_H_gravityCenter.size()<<sendl;
    }

    for (unsigned int i=0; i<vec_global_H_gravityCenter.size(); i++)
    {

        Real Length = rigidCurveSegments[i].second - rigidCurveSegments[i].first;
        Vec3 sizeArrows (Length/4, Length/8, Length/8);

        vparams->drawTool()->drawFrame(vec_global_H_gravityCenter[i].getOrigin(), vec_global_H_gravityCenter[i].getOrientation(), sizeArrows );

    }



}


} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_SutureController_INL */
