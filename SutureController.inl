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
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/component/topology/EdgeSetGeometryAlgorithms.h>
#include <sofa/core/behavior/MechanicalState.h>

//#define DEBUG


namespace sofa
{

namespace component
{

namespace controller
{

template <class DataTypes>
SutureController<DataTypes>::SutureController(WireBeamInterpolation<DataTypes>* _adaptiveinterpolation)
: m_adaptiveinterpolation(_adaptiveinterpolation)
, startingPos(initData(&startingPos,Coord(),"startingPos","starting pos for inserting the instrument"))
, threshold(initData(&threshold, (Real)0.000001, "threshold", "threshold for controller precision which is homogeneous to the unit of length"))
, maxBendingAngle(initData(&maxBendingAngle, (Real)0.1, "maxBendingAngle", "max bending criterion (in rad) for one beam"))
, m_interpolationPath(initData(&m_interpolationPath,"interpolation", "Path to the Interpolation component on scene"))
, useDummyController(initData(&useDummyController, false, "useDummyController"," use a very simple controller of adaptativity (use for debug)" ))
{


}

template <class DataTypes>
SutureController<DataTypes>::SutureController()
: m_adaptiveinterpolation(NULL)
, startingPos(initData(&startingPos,Coord(),"startingPos","starting pos for inserting the instrument"))
, threshold(initData(&threshold, (Real)0.000001, "threshold", "threshold for controller precision which is homogeneous to the unit of length"))
, maxBendingAngle(initData(&maxBendingAngle, (Real)0.1, "maxBendingAngle", "max bending criterion (in rad) for one beam"))
, m_interpolationPath(initData(&m_interpolationPath,"interpolation", "Path to the Interpolation component on scene"))
, useDummyController(initData(&useDummyController, false, "useDummyController"," use a very simple controller of adaptativity (use for debug)" ))
{


}


template <class DataTypes>
void SutureController<DataTypes>::init()
{
        if (m_adaptiveinterpolation==NULL) {
            ///////// get the Adaptive Interpolation component ///////
            //std::vector<sofa::core::behavior::LinearSolver*> solvers;
            core::objectmodel::BaseContext * c = this->getContext();

            const helper::vector<std::string>& interpolName = m_interpolationPath.getValue();
            if (interpolName.empty()) {
                m_adaptiveinterpolation = c->get<WInterpolation>(core::objectmodel::BaseContext::Local);
            } else {
                m_adaptiveinterpolation = c->get<WInterpolation>(m_interpolationPath.getValue()[0]);
            }

            if(m_adaptiveinterpolation==NULL)
                serr<<" no Beam Interpolation found !!! the component can not work"<<sendl;
            else
                sout<<" interpolation named"<<m_adaptiveinterpolation->getName()<<" found (for "<<this->getName()<<")"<<sendl;
        }


	xAbs_collisionPoints_buf.clear();


	this->f_listening.setValue(true);

	Inherit::init();
        m_adaptiveinterpolation->setControlled(true);


	// TODO : VERIFIER QUE LA TOPOLOGIE EST CELLE D'UN WIRE

        this->getContext()->get(_topology);
        _topology->clear();
        _topology->cleanup();


        // on initialise le "wire" en prenant la position de départ + la forme au repos + la discretisation proposée...
        Transform global_T_init;
        Coord startPos=startingPos.getValue();
        global_T_init.setOrigin( startPos.getCenter() );
        global_T_init.setOrientation(startPos.getOrientation());


        helper::vector<Real> xP_noticeable;
        helper::vector< int> nbP_density;
        m_adaptiveinterpolation->getSamplingParameters( xP_noticeable, nbP_density);

        // computation of the number of node on the structure:
        unsigned int numNodes;
        numNodes=1; //initial point
        for (unsigned int i=0; i<nbP_density.size(); i++)
        {
            numNodes+=nbP_density[i]; // numBeams between each noticeable point
        }




        // Initial position of the nodes:
        nodeCurvAbs.clear();

        Real x_curv=0.0;



        Data<VecCoord>* datax = this->getMechanicalState()->write(sofa::core::VecCoordId::position());
        Data<VecDeriv>* datav = this->getMechanicalState()->write(sofa::core::VecDerivId::velocity());
        VecCoord& x = *datax->beginEdit();
        VecDeriv& v = *datav->beginEdit();

        this->getMechanicalState()->resize( numNodes);

        x.resize( numNodes);
        v.resize( numNodes );



        x[0].getCenter()= global_T_init.getOrigin();
        x[0].getOrientation()= global_T_init.getOrientation();

        unsigned int id_node=0;
        _topology->addPoint( x[0][0], x[0][1], x[0][2] );
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



            Real dx= length/nbP_density[i];



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

                _topology->addEdge( (int)(id_node-1), (int)(id_node) );
                _topology->addPoint( x[id_node][0], x[id_node][1], x[id_node][2] );

                nodeCurvAbs.push_back(x_curv);


                // add the beam to the m_adaptiveinterpolation
                m_adaptiveinterpolation->addBeam(id_node-1,dx,x_curv-dx,x_curv,0.0  );

#ifdef DEBUG
                std::cout<<" Pos at x_curv ="<<x_curv<<" : "<<x[id_node]<<std::endl;
#endif
            }

        }


        this->getContext()->get(edgeMod);

        if (edgeMod == NULL)
                serr << "EdgeSetController has no binding EdgeSetTopologyModifier." << sendl;



        datax->endEdit();
        datav->endEdit();

}

template <class DataTypes>
void SutureController<DataTypes>::reinit()
{

    this->getMechanicalState()->cleanup();
    init();
    applyController();

}


///////////// TODO: Faire une fonction qui recalcule entièrement la topologie (pour remplacer addNodesAndEdge et removeNodesAndEdge


template <class DataTypes>
void SutureController<DataTypes>::recreateTopology()
{

    /* A chaque pas de temps, ON REMET A JOUR LA TOPO ENTIEREMENT */
    _topology->cleanup();
    _topology->clear();





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
        _topology->addPoint( 0.0, 0.0, 0.0 );
        _topology->addEdge( (int)(numNodes-1+i), (int)(numNodes+i) );
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
    _topology->cleanup();
    _topology->clear();

    _topology->addPoint(0.0,0.0,0.0);
    for (int i=0; i<(int) (numBeams-num); i++)
    {
        _topology->addPoint((double)(i+1),0.0,0.0);
        _topology->addEdge(i,i+1);
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


    static int compteur=0;
    Real length = nodeCurvAbs[nodeCurvAbs.size()-1];
    newCurvAbs.clear();

    for (unsigned int i=0; i<=10; i++)
    {
        Real xtest=length * (Real)i /10.0;
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
void SutureController<DataTypes>::applyController()
{


//    std::cout<<"applyController: numBeams =" <<m_adaptiveinterpolation->getNumBeams()<<std::endl;
    Data<VecCoord>* datax = this->getMechanicalState()->write(sofa::core::VecCoordId::position());
    Data<VecDeriv>* datav = this->getMechanicalState()->write(sofa::core::VecDerivId::velocity());
    VecCoord& x = *datax->beginEdit();
    VecDeriv& v = *datav->beginEdit();


    sofa::helper::vector<Real> newCurvAbs;

    if (useDummyController.getValue())
        this->dummyController(newCurvAbs);
    else
        this->computeSampling(newCurvAbs, x);


#ifdef DEBUG
    std::cout<<" newCurvAbs = "<<newCurvAbs<<"  - nodeCurvAbs = "<<nodeCurvAbs<<std::endl;

    for (unsigned int i=0; i<rigidCurveSegments.size(); i++)
    {
        std::cout<<" rigidCurveSegments["<<i<<"] "<<rigidCurveSegments[0].first<<" - "<<rigidCurveSegments[0].second<<std::endl;
    }

    std::cout<<" Before... Beam  lengths = ";
    for (unsigned int b=0; b<m_adaptiveinterpolation->getNumBeams(); b++)
        std::cout<<" "<<m_adaptiveinterpolation->getLength(b);
    std::cout<<" "<<std::endl;

#endif


    this->applyNewSampling(newCurvAbs, nodeCurvAbs, x, v);

 #ifdef DEBUG
    std::cout<<"After.... Beam  lengths = ";

    for (unsigned int b=0; b<m_adaptiveinterpolation->getNumBeams(); b++)
        std::cout<<" "<<m_adaptiveinterpolation->getLength(b);
    std::cout<<" "<<std::endl;
#endif

    nodeCurvAbs.clear();
    for( unsigned int i=0; i< newCurvAbs.size(); i++)
    {
        nodeCurvAbs.push_back(newCurvAbs[i]);

    }

    datax->endEdit();
    datav->endEdit();



}


//////*************** PRIVATE FUNCTIONS ****************//

// this function calls ComputeTotalBendingRotationAngle on the beams between xmin and xmax
template <class DataTypes>
void SutureController<DataTypes>::computeBendingAngle(Real& angle, const Real& xmin, const Real& xmax, const Real& dx_comput, const VecCoord& Pos)
{

    //test for verification
    if (xmax < xmin || dx_comput==0.0){
        serr<<"wrong parameters in computeBendingAngle function"<<sendl;
        angle=0.0;
        return;
    }

    if (xmax > m_adaptiveinterpolation->getRestTotalLength()){
        serr<<" in computeBendingAngle : max > getRestTotalLength"<<sendl;
        angle=0.0;
        return;
    }


    unsigned int idBeamMin, idBeamMax;
    Real baryCoordMin, baryCoordMax;
    Transform Tnode0, Tnode1;

    m_adaptiveinterpolation->getBeamAtCurvAbs(xmin, idBeamMin, baryCoordMin);
    m_adaptiveinterpolation->getBeamAtCurvAbs(xmax, idBeamMax, baryCoordMax);
    m_adaptiveinterpolation->computeTransform2(idBeamMin,Tnode0,Tnode1, Pos);


    angle=0.0;

    if (idBeamMin==idBeamMax)
    {
        m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(angle, dx_comput, Tnode0, Tnode1, m_adaptiveinterpolation->getLength(idBeamMin), baryCoordMin, baryCoordMax);
         //std::cout<<" idBeamMin==idBeamMax, angle ="<<angle<<std::endl;
        return;
    }

    if (idBeamMin>idBeamMax)
    {
        serr<<"ERROR in computeBendingAngle"<<sendl;
        return;
    }

    //////////


    // compute the angle for the first beam

    Real angleBeam=0.0;
    m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(angleBeam, dx_comput, Tnode0, Tnode1, m_adaptiveinterpolation->getLength(idBeamMin), baryCoordMin, 1.0);
    angle+=angleBeam;

    unsigned int b=idBeamMin+1;
    while(b<idBeamMax)
    {
        m_adaptiveinterpolation->computeTransform2(b,Tnode0,Tnode1,Pos);
        m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(angleBeam, dx_comput, Tnode0, Tnode1, m_adaptiveinterpolation->getLength(b), 0.0, 1.0);
        angle+=angleBeam;
        b++;
    }

    m_adaptiveinterpolation->computeTransform2(idBeamMax,Tnode0,Tnode1,Pos);
    m_adaptiveinterpolation->ComputeTotalBendingRotationAngle(angleBeam, dx_comput, Tnode0, Tnode1, m_adaptiveinterpolation->getLength(idBeamMax), 0.0, baryCoordMax);

    angle+=angleBeam;

}


// this function computes the tangent value on a series of discrete points (store also the curv_abs of these discrete points)
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
void SutureController<DataTypes>::detectRigidBeams(sofa::helper::vector<Real> &newCurvAbs)
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
void SutureController<DataTypes>::applyNewSampling(sofa::helper::vector<Real> &newCurvAbs, sofa::helper::vector<Real> &oldCurvAbs, VecCoord &x, VecDeriv &v)
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
    _topology->addPoint( xDof[0], xDof[1], xDof[2] );
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
            _topology->addEdge( (int)(numNodes-1), (int)(numNodes-1) );

            L =  newCurvAbs[s+1] - newCurvAbs[s];
            m_adaptiveinterpolation->addBeam(s, L, newCurvAbs[s], newCurvAbs[s+1] ,GravityCenter_H_interpol0, GravityCenter_H_interpol1 );

#ifdef DEBUG
            std::cout<<" add rigid segment on segment "<<s<<" that links the same dof: [ "<<_topology->getEdge(s)[0]<<" "<<_topology->getEdge(s)[1]
                    <<"]   x curv = ["<< newCurvAbs[s]<<" " <<newCurvAbs[s+1]<<"]"<<std::endl;
            std::cout<<"GravityCenter_H_interpol0 ="<<GravityCenter_H_interpol0<<"   - GravityCenter_H_interpol1"<<GravityCenter_H_interpol1<<std::endl;
#endif
        }

        else // the beam is deformable
        {


            L =  newCurvAbs[s+1] - newCurvAbs[s];

            // ADD the beam in the topology and in the interpolation

            _topology->addEdge( (int)(numNodes-1), (int)(numNodes) );
            m_adaptiveinterpolation->addBeam(s, L, newCurvAbs[s], newCurvAbs[s+1] ,0.0 );
#ifdef DEBUG
            std::cout<<" add deformable segment on segment "<<s<<" that links the same dof: [ "<<_topology->getEdge(s)[0]<<" "<<_topology->getEdge(s)[1]
                    <<"]   x curv = ["<< newCurvAbs[s]<<" " <<newCurvAbs[s+1]<<"]"<<std::endl;
#endif
            // ADD a DOF for the second node of the beam
            numNodes++;
            xDof.getCenter()     = vec_global_H_node[s+1].getOrigin();
            xDof.getOrientation()= vec_global_H_node[s+1].getOrientation();
            x.push_back(xDof);

            v.push_back(vec_global_Vel_node[s+1]);

            _topology->addPoint( xDof[0], xDof[1], xDof[2] );


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


    if(!verifyRigidCurveSegmentSort()){
        serr<<" WARNING : rigidCurveSegments are not correctly sorted !"<<sendl;
    }


    Real sutureLength= m_adaptiveinterpolation->getRestTotalLength();
    Real dx=sutureLength/100.0;
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
                this->computeBendingAngle(BendingAngle, L_add, L_add+L_beam_straight, dx, x);
            else
                this->computeBendingAngle(BendingAngle, L_add, xP_noticeable[part+1], dx, x);




            Real L_beam = L_beam_straight;
            if (BendingAngle>maxBendingAngle.getValue())
            {
                L_beam*=maxBendingAngle.getValue()/BendingAngle;
            }
            L_add+=L_beam;

        }
        newCurvAbs_notSecure.push_back(xP_noticeable[part+1]);

    }







    this->rigidCurveSegments.clear();
    /*

    std::pair<Real, Real> rigidSegment;

    rigidSegment.first=8.5*sutureLength/10.0;
    rigidSegment.second = 9.5*sutureLength/10.0;
    this->rigidCurveSegments.push_back(rigidSegment);

    addRigidCurvAbs(newCurvAbs_notSecure, threshold.getValue());

    */


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
