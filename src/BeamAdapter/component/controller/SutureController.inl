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
#pragma once

#include <sofa/core/visual/VisualParams.h>
#include <sofa/simulation/MechanicalVisitor.h>
#include <sofa/simulation/UpdateMappingVisitor.h>

#include <BeamAdapter/component/controller/SutureController.h>
#include <BeamAdapter/component/WireBeamInterpolation.h>


namespace sofa::component::controller::_suturecontroller_
{

using sofa::core::objectmodel::BaseContext ;

template <class DataTypes>
SutureController<DataTypes>::SutureController(WireBeamInterpolation<DataTypes>* _adaptiveinterpolation)
: d_startingPos(initData(&d_startingPos,Coord(),"startingPos","starting pos for inserting the instrument"))
, d_threshold(initData(&d_threshold, (Real)0.000001, "threshold", "threshold for controller precision which is homogeneous to the unit of length"))
, d_maxBendingAngle(initData(&d_maxBendingAngle, (Real)0.1, "maxBendingAngle", "max bending criterion (in rad) for one beam"))
, d_useDummyController(initData(&d_useDummyController, false, "useDummyController"," use a very simple controller of adaptativity (use for debug)" ))
, d_fixRigidTransforms(initData(&d_fixRigidTransforms, false, "fixRigidTransforms", "fix the sampling and transformations of rigid segments"))
, d_rigidCurvAbs(initData(&d_rigidCurvAbs, "rigidCurvAbs", "pairs of curv abs for beams we want to rigidify"))
, l_adaptiveInterpolation(initLink("interpolation", "Path to the Interpolation component on scene"), _adaptiveinterpolation)
, d_nodeCurvAbs(initData(&d_nodeCurvAbs, "nodeCurvAbs", ""))
, d_curvatureList(initData(&d_curvatureList, "curvatureList", "List of the beams curvature (abscissa - curvature)"))
, d_controlPoints(initData(&d_controlPoints, "controlPoints", "List of the spline control points positions"))
, d_updateOnBeginAnimationStep(initData(&d_updateOnBeginAnimationStep, false, "updateOnBeginAnimationStep", "If true update interpolation and subgraph on beginAnimationStep"))
, d_applyOrientationFirstInCreateNeedle(initData(&d_applyOrientationFirstInCreateNeedle, false, "applyOrientationFirstInCreateNeedle", "if true, it sets first the orientation, then the rotation for a init node of the needle"))
, d_reinitilizeWireOnInit(initData(&d_reinitilizeWireOnInit, false, "reinitilizeWireOnInit"," reinitialize the wire everytime init() is called (for planning purposes)" ))
, d_actualStepNoticeablePoints(initData(&d_actualStepNoticeablePoints, "m_actualStepNoticeablePoints", "points (as curv. absc.) that are to be considered when computing a new sampling in the actual time step"))
, m_topology(0)
{
}

template <class DataTypes>
void SutureController<DataTypes>::init()
{
    f_listening.setValue(true);
    Inherit1::init();

    if (!l_adaptiveInterpolation)
    {
        BaseContext *c=getContext();
        l_adaptiveInterpolation.set(c->get<WInterpolation>(BaseContext::Local));
    }

    if(!l_adaptiveInterpolation)
        msg_error() << "No Beam Interpolation found, the component can not work.";

    l_adaptiveInterpolation->setControlled(true);

    m_XAbsCollisionPointsBuffer.clear();

    getContext()->get(m_topology);

    if (d_reinitilizeWireOnInit.getValue() || !wireIsAlreadyInitialized())
        initWireModel();
}


template <class DataTypes>
void SutureController<DataTypes>::reinit()
{
    getMechanicalState()->cleanup();
    init();
    applyController();
    updateControlPointsPositions();
}


template <class DataTypes>
void SutureController<DataTypes>::initWireModel()
{
    m_topology->clear();
    m_topology->cleanup();

    // on initialise le "wire" en prenant la position de départ + la forme au repos + la discretisation proposée...
    Transform global_T_init;
    const Coord startPos = d_startingPos.getValue();

    if (d_applyOrientationFirstInCreateNeedle.getValue()) {
        global_T_init.setOrientation(startPos.getOrientation());
        global_T_init.setOrigin(startPos.getCenter());
    }
    else {
        global_T_init.setOrigin(startPos.getCenter());
        global_T_init.setOrientation(startPos.getOrientation());
    }

    type::vector< Real > xP_noticeable;
    type::vector< int > nbP_density;
    l_adaptiveInterpolation->getSamplingParameters(xP_noticeable, nbP_density);

    // computation of the number of node on the structure:
    unsigned int numNodes = 1;

    for (unsigned int i=0; i<nbP_density.size(); i++)
        numNodes += nbP_density[i]; // numBeams between each noticeable point

    auto nodeCurvAbs = sofa::helper::getWriteOnlyAccessor(d_nodeCurvAbs);

    // Initial position of the nodes:
    nodeCurvAbs.clear();

    Real x_curv = 0.0;

    Data<VecCoord>* datax = getMechanicalState()->write(sofa::core::VecCoordId::position());
    Data<VecDeriv>* datav = getMechanicalState()->write(sofa::core::VecDerivId::velocity());
    auto x = sofa::helper::getWriteOnlyAccessor(*datax);
    auto v = sofa::helper::getWriteOnlyAccessor(*datav);

    getMechanicalState()->resize(numNodes);

    x.resize(numNodes);
    v.resize(numNodes);

    x[0].getCenter()= global_T_init.getOrigin();
    x[0].getOrientation()= global_T_init.getOrientation();

    unsigned int id_node=0;
    m_topology->addPoint(x[0][0], x[0][1], x[0][2]);
    nodeCurvAbs.push_back(0.0);

    l_adaptiveInterpolation->clear();

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
            l_adaptiveInterpolation->getRestTransformOnX(global_H_localP, x_curv);

            Transform global_H_P = global_T_init * global_H_localP;

            x[id_node].getCenter()= global_H_P.getOrigin();
            x[id_node].getOrientation()= global_H_P.getOrientation();

            // modif the topology
            m_topology->addEdge( (int)(id_node-1), (int)(id_node) );
            m_topology->addPoint( x[id_node][0], x[id_node][1], x[id_node][2] );

            nodeCurvAbs.push_back(x_curv);

            // add the beam to the m_adaptiveinterpolation
            l_adaptiveInterpolation->addBeam(id_node-1,dx,x_curv-dx,x_curv,0.0  );
        }

    }
}


template <class DataTypes>
bool SutureController<DataTypes>::wireIsAlreadyInitialized()
{
    unsigned int numDofs = 0;

    if (getMechanicalState() != nullptr)
    {
        numDofs = getMechanicalState()->getSize();
        if (numDofs == 0)
            return false;
    }
    else
    {
        msg_error() << "SutureController should have a MechanicalState in its context";
        return false;
    }

    /// When there are rigid segments, # of dofs is different than # of edges and beams
    unsigned int numRigidPts = 0;
    helper::ReadAccessor< Data< std::set< Real > > > rigidCurvAbs = d_rigidCurvAbs;
    const auto nbRigidAbs = rigidCurvAbs->size();
    if (nbRigidAbs>0 && (nbRigidAbs%2)==0)
    {
        const type::vector<Real>& curvAbs = d_nodeCurvAbs.getValue();
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

    if (m_topology != nullptr)
    {
        if ((unsigned int)m_topology->getNbPoints() != numDofs || (unsigned int)m_topology->getNbEdges() != (numDofs+numRigidPts-1))
            return false;
    }
    else
    {
        msg_error() << "SutureController should have a topology container in its context";
        return false;
    }

    if (l_adaptiveInterpolation != nullptr)
    {
        if (l_adaptiveInterpolation->getNumBeams() != (numDofs + numRigidPts - 1))
            return false;
    }
    else
    {
         msg_error() << "SutureController should have a WireBeamInterpolation in its context";
        return false;
    }

    if (d_nodeCurvAbs.getValue().size() != (numDofs + numRigidPts))
        return false;

    return true;
}


template <class DataTypes>
void SutureController<DataTypes>::recreateTopology()
{
    m_topology->cleanup();
    m_topology->clear();
}



template <class DataTypes>
void SutureController<DataTypes>::addNodesAndEdge(unsigned int num, Real &xend)
{
    unsigned int numNodes = getMechanicalState()->getSize();
    getMechanicalState()->resize( numNodes + num);
    unsigned int numBeams = l_adaptiveInterpolation->getNumBeams();

    for (unsigned int i=0; i<num; i++)
    {
        m_topology->addPoint( 0.0, 0.0, 0.0 );
        m_topology->addEdge( (int)(numNodes-1+i), (int)(numNodes+i) );
        l_adaptiveInterpolation->addBeam(numBeams-1+i,xend,0.0,xend,0.0  ); // the parameters will be coerrected in applyNewSampling
    }
}


template <class DataTypes>
void SutureController<DataTypes>::removeNodesAndEdge(unsigned int num)
{
    unsigned int numNodes = getMechanicalState()->getSize();
    unsigned int numBeams = l_adaptiveInterpolation->getNumBeams();

    m_topology->cleanup();
    m_topology->clear();

    m_topology->addPoint(0.0,0.0,0.0);
    for (int i=0; i<(int) (numBeams-num); i++)
    {
        m_topology->addPoint((double)(i+1),0.0,0.0);
        m_topology->addEdge(i,i+1);
    }

    getMechanicalState()->resize( numNodes - num);
}


template <class DataTypes>
void SutureController<DataTypes>::onBeginAnimationStep(const double dt)
{
    SOFA_UNUSED(dt);
    if (d_updateOnBeginAnimationStep.getValue())
        applyController();

    // Propagate modifications
    MechanicalProjectPositionAndVelocityVisitor(core::MechanicalParams::defaultInstance(), getContext()->getTime(), sofa::core::VecCoordId::position(),sofa::core::VecDerivId::velocity()); // apply projective constraints
    MechanicalPropagateOnlyPositionAndVelocityVisitor(core::MechanicalParams::defaultInstance(), getContext()->getTime(),sofa::core::VecCoordId::position(),sofa::core::VecDerivId::velocity()).execute( getContext() );
    simulation::UpdateMappingVisitor(core::ExecParams::defaultInstance()).execute(getContext());
}


template <class DataTypes>
void SutureController<DataTypes>::onEndAnimationStep(const double /*dt*/)
{
    if (!d_updateOnBeginAnimationStep.getValue())
    {
        applyController();
    }

    updateControlPointsPositions();
}


template <class DataTypes>
void SutureController<DataTypes>::dummyController(type::vector<Real> &newCurvAbs)
{
    static int compteur = 0;
    Real length = d_nodeCurvAbs.getValue()[d_nodeCurvAbs.getValue().size() - 1];
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

    // between compteur 10 and 90 A node is added
    if( compteur > 10 && compteur < 90)
    {
        Real totalLength = newCurvAbs[newCurvAbs.size()-1];
        newCurvAbs[newCurvAbs.size()-1] = totalLength + decalage;
        newCurvAbs.push_back(totalLength );
    }

    compteur++;
    if (compteur==100)
        compteur=0;


    /// Rigidify a part of the beam
    m_rigidCurveSegments.clear();

    std::pair<Real, Real> rigidSegment;

    rigidSegment.first=3.2*length/10.0;
    rigidSegment.second = 5.5*length/10.0;
    m_rigidCurveSegments.push_back(rigidSegment);

    rigidSegment.first = 7.5*length/10.0;
    rigidSegment.second = 10*length/10.0;
    m_rigidCurveSegments.push_back(rigidSegment);

    addRigidCurvAbs(newCurvAbs, 0.0001);

    Real totalLength = newCurvAbs[newCurvAbs.size()-1];
    m_listOfImposedNodesOnXcurv.push_back(0.23*totalLength);
    m_listOfImposedNodesOnXcurv.push_back(0.57*totalLength);
}


template <class DataTypes>
void SutureController<DataTypes>::addRigidCurvAbs(type::vector<Real> &newCurvAbs, const Real &tol)
{
    type::vector<Real> newCurvAbsBuf=newCurvAbs;

    newCurvAbs.clear();

    unsigned int iterator=1;
    newCurvAbs.push_back(newCurvAbsBuf[0]);

    for (unsigned int i=0; i<m_rigidCurveSegments.size(); i++)
    {

        while(m_rigidCurveSegments[i].first > newCurvAbsBuf[iterator])
        {
            newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
            iterator++;
        }

        // newCurvAbsBuf[iterator-1]  <  rigidCurveSegments[i].first <  newCurvAbsBuf[iterator]
        //=>  if rigidCurveSegments[i].first ~ newCurvAbsBuf[iterator-1]
        // or if rigidCurveSegments[i].first ~ newCurvAbsBuf[iterator]  => do not add the point
        if (m_rigidCurveSegments[i].first - newCurvAbsBuf[iterator-1] > tol && newCurvAbsBuf[iterator] - m_rigidCurveSegments[i].first > tol)
        {
            newCurvAbs.push_back(m_rigidCurveSegments[i].first);
        }



        while(m_rigidCurveSegments[i].second > newCurvAbsBuf[iterator])
        {
            newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
            iterator++;
        }

        // newCurvAbsBuf[iterator-1]  <  rigidCurveSegments[i].second <  newCurvAbsBuf[iterator]
        //=>  if rigidCurveSegments[i].second ~ newCurvAbsBuf[iterator-1]
        // or if rigidCurveSegments[i].second ~ newCurvAbsBuf[iterator]  => do not add the point
        if (m_rigidCurveSegments[i].second - newCurvAbsBuf[iterator-1] > tol && newCurvAbsBuf[iterator] - m_rigidCurveSegments[i].second > tol)
        {
            newCurvAbs.push_back(m_rigidCurveSegments[i].second);
        }

    }


    while(iterator< newCurvAbsBuf.size())
    {
        newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
        iterator++;
    }
}



template <class DataTypes>
void SutureController<DataTypes>::addImposedCurvAbs(type::vector<Real> &newCurvAbs, const Real &tol)
{
    m_listOfImposedNodesOnXcurv.sort();
    m_listOfImposedNodesOnXcurv.unique();

    ListRealIterator it_xcurv_imposed;

    type::vector<Real> newCurvAbsBuf=newCurvAbs;
    newCurvAbs.clear();
    unsigned int iterator=1;
    newCurvAbs.push_back(newCurvAbsBuf[0]);

    for (it_xcurv_imposed=m_listOfImposedNodesOnXcurv.begin(); it_xcurv_imposed!=m_listOfImposedNodesOnXcurv.end(); it_xcurv_imposed++)
    {
        while( (*it_xcurv_imposed) > newCurvAbsBuf[iterator]) // newCurvAbsBuf=  [ 0 2 5 7]    // xcurvImposed= [ 2.1  3  4.9 ]
        {
            newCurvAbs.push_back(newCurvAbsBuf[iterator]) ;
            iterator++;
        }

        if ( newCurvAbsBuf[iterator] - (*it_xcurv_imposed)  > tol  &&  (*it_xcurv_imposed) - newCurvAbsBuf[iterator-1] > tol) // cas 2 par exemple =>  xcurvImposed=3  newCurvAbs[iterator]= 5
        {
            newCurvAbs.push_back( (*it_xcurv_imposed) ); // newCurvAbsBuf=  [ 0 2 3 ... ]
        }
        else if( (*it_xcurv_imposed) - newCurvAbsBuf[iterator-1] < tol ) //  cas 1 par exemple => xcurvImposed=2.1  newCurvAbs
        {
            newCurvAbs[iterator-1] = (*it_xcurv_imposed); // newCurvAbsBuf= [0 2.1 ...]
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
}



template <class DataTypes>
void SutureController<DataTypes>::applyController()
{
    Data<VecCoord>* datax = getMechanicalState()->write(sofa::core::VecCoordId::position());
    Data<VecDeriv>* datav = getMechanicalState()->write(sofa::core::VecDerivId::velocity());
    auto x = sofa::helper::getWriteOnlyAccessor(*datax);
    auto v = sofa::helper::getWriteOnlyAccessor(*datav);
    type::vector<Real> newCurvAbs;

    storeRigidSegmentsTransformations();

    if (d_useDummyController.getValue())
        dummyController(newCurvAbs);
    else
        computeSampling(newCurvAbs, x.wref());

    addImposedCurvAbs(newCurvAbs, 0.0001);

    verifyRigidSegmentsSampling(newCurvAbs);

    applyNewSampling(newCurvAbs, d_nodeCurvAbs.getValue(), x.wref(), v.wref());
    verifyRigidSegmentsTransformations();
    m_prevRigidCurvSegments = m_rigidCurveSegments;

    auto nodeCurvAbs = sofa::helper::getWriteOnlyAccessor(d_nodeCurvAbs);
    nodeCurvAbs.wref().assign(newCurvAbs.begin(), newCurvAbs.end());
}


//////*************** PRIVATE FUNCTIONS ****************//

template <class DataTypes>
typename SutureController<DataTypes>::Real SutureController<DataTypes>::computeBendingAngle(const Real& xmin, const Real& xmax, const Real& dx_comput, const VecCoord& Pos)
{
    //test for verification
    if (xmax < xmin || dx_comput==0.0){
        dmsg_error()<<"Wrong parameters in computeBendingAngle function";
        return 0.0;
    }

    if (xmax > l_adaptiveInterpolation->getRestTotalLength()){
        dmsg_error()<<"In computeBendingAngle : max > getRestTotalLength";
        return 0.0;
    }


    unsigned int idBeamMin, idBeamMax;
    Real baryCoordMin, baryCoordMax;
    Transform Tnode0, Tnode1;

    l_adaptiveInterpolation->getBeamAtCurvAbs(xmin, idBeamMin, baryCoordMin);
    l_adaptiveInterpolation->getBeamAtCurvAbs(xmax, idBeamMax, baryCoordMax);
    l_adaptiveInterpolation->computeTransform(idBeamMin,Tnode0,Tnode1, Pos);

    if (idBeamMin==idBeamMax)
        return l_adaptiveInterpolation->ComputeTotalBendingRotationAngle(dx_comput, Tnode0, Tnode1, l_adaptiveInterpolation->getLength(idBeamMin), baryCoordMin, baryCoordMax);

    if (idBeamMin>idBeamMax)
    {
        dmsg_error()<<"In computeBendingAngle";
        return 0.0;
    }

    // compute the angle for the first beam
    Real angle = l_adaptiveInterpolation->ComputeTotalBendingRotationAngle(dx_comput, Tnode0, Tnode1, l_adaptiveInterpolation->getLength(idBeamMin), baryCoordMin, 1.0);

    unsigned int b=idBeamMin+1;
    while(b<idBeamMax)
    {
        l_adaptiveInterpolation->computeTransform(b,Tnode0,Tnode1,Pos);
        angle += l_adaptiveInterpolation->ComputeTotalBendingRotationAngle(dx_comput, Tnode0, Tnode1, l_adaptiveInterpolation->getLength(b), 0.0, 1.0);
        b++;
    }

    l_adaptiveInterpolation->computeTransform(idBeamMax,Tnode0,Tnode1,Pos);
    angle += l_adaptiveInterpolation->ComputeTotalBendingRotationAngle(dx_comput, Tnode0, Tnode1, l_adaptiveInterpolation->getLength(idBeamMax), 0.0, baryCoordMax);

    return angle;
}


template <class DataTypes>
void SutureController<DataTypes>::computeTangentOnDiscretePoints(type::vector<Vec3> TangTable, type::vector<Real> xTable,  unsigned int numDiscretePoints, const VecCoord& Pos)
{

    TangTable.clear();
    Transform Tnode0, Tnode1;
    Real baryCoord;
    unsigned int beam;

    Real dx = l_adaptiveInterpolation->getRestTotalLength()/(numDiscretePoints-1);
    Real x=dx;

    // compute intial tang for the beginning of the wire:
    l_adaptiveInterpolation->computeTransform(0,Tnode0,Tnode1, Pos);
    Vec3 t = Tnode0.getOrientation().rotate(Vec3(1.0,0.0,0.0));
    TangTable.push_back(t);
    xTable.push_back(0.0);

    for (unsigned int p=0; p<(numDiscretePoints-1) ; p++)
    {
        l_adaptiveInterpolation->getBeamAtCurvAbs(x, beam, baryCoord);
        l_adaptiveInterpolation->computeTransform(beam,Tnode0,Tnode1, Pos);
        l_adaptiveInterpolation->getTangent(t, baryCoord, Tnode0,Tnode1,l_adaptiveInterpolation->getLength(beam) );

        TangTable.push_back(t);
        xTable.push_back(x);

        x+=dx;

    }
}

template <class DataTypes>
void SutureController<DataTypes>::detectRigidBeams(const type::vector<Real> &newCurvAbs)
{
    unsigned int seg=0;

    m_rigidBeamList.clear();

    for (unsigned int i=1; i<newCurvAbs.size(); i++)
    {
        /// Rigidification
        //(1) look if there is a rigidCurveSegment...
        if(seg < m_rigidCurveSegments.size() )
        {
            // (2) look if [newCurvAbs[i-1] newCurvAbs[i] are contained in the next rigidified segment....
            if(newCurvAbs[i-1]+d_threshold.getValue() > m_rigidCurveSegments[seg].first && newCurvAbs[i]-d_threshold.getValue() < m_rigidCurveSegments[seg].second )
                m_rigidBeamList.push_back(true);
            else
                m_rigidBeamList.push_back(false);

            if( newCurvAbs[i]+d_threshold.getValue() > m_rigidCurveSegments[seg].second)
                seg++;
        }
        else
            m_rigidBeamList.push_back(false);

    }
}



// When a new sampling is defined in "newCurvAbs", the position and the velocity needs to be "re-interpolated"

template <class DataTypes>
void SutureController<DataTypes>::applyNewSampling(const type::vector<Real> &newCurvAbs, const type::vector<Real> &oldCurvAbs, VecCoord &x, VecDeriv &v)
{
    VecCoord x_buf=x;
    VecDeriv v_buf=v;

    x.clear();
    v.clear();

    detectRigidBeams(newCurvAbs);


    /// interpolation of the position and velocities
    Transform global_H_interpol;
    Deriv v_interpol;
    unsigned int j=0;
    Vec3 null(0,0,0);

    m_vecGlobalHNode.clear();
    m_vecGlobalVelNode.clear();

    l_adaptiveInterpolation->InterpolateTransformAndVelUsingSpline(0,0.0,null,x_buf, v_buf, global_H_interpol, v_interpol);
    m_vecGlobalHNode.push_back(global_H_interpol);
    m_vecGlobalVelNode.push_back(v_interpol);


    for (unsigned int i=1; i<newCurvAbs.size()-1; i++)
    {

        while(newCurvAbs[i]>oldCurvAbs[j])
        {
                j++;
                if (j>=oldCurvAbs.size())
                {
                        dmsg_warning() << " j ="<<j<<">=oldCurvAbs.size()";
                        return;
                }
        }
        Real L = l_adaptiveInterpolation->getLength(j-1);
        Real L0 = newCurvAbs[i] - oldCurvAbs[j-1];
        Real ratio=L0/L;

        l_adaptiveInterpolation->InterpolateTransformAndVelUsingSpline(j-1,ratio,null,x_buf, v_buf, global_H_interpol, v_interpol);
        m_vecGlobalHNode.push_back(global_H_interpol);
        m_vecGlobalVelNode.push_back(v_interpol);

    }
    l_adaptiveInterpolation->InterpolateTransformAndVelUsingSpline(oldCurvAbs.size()-2,1.0,null,x_buf, v_buf, global_H_interpol, v_interpol);
    m_vecGlobalHNode.push_back(global_H_interpol);
    m_vecGlobalVelNode.push_back(v_interpol);

    /// Compute a gravity center for the rigid segments
    m_vecGlobalHGravityCenter.clear();
    type::vector<Deriv> vec_Vel_gravityCenter;


    Transform global_H_gravityC;
    Deriv vel_gravityC;
    bool rigidification;
    rigidification=false;
    unsigned int Rseg=0;
    Real length_of_rigidSegment=0;


    for (unsigned int s=0; s<m_rigidBeamList.size();s++)
    {
        if(m_rigidBeamList[s]) // the beam is rigidified
        {
            if (!rigidification) //begining of the rigidification
            {

                length_of_rigidSegment = m_rigidCurveSegments[Rseg].second  - m_rigidCurveSegments[Rseg].first;
                Rseg++;
                rigidification=true;
                global_H_gravityC.clear();

            }

            Real length_of_beam = newCurvAbs[s+1]-newCurvAbs[s];
            Real weight = (length_of_beam/(length_of_rigidSegment*2));

            global_H_gravityC.setOrigin( m_vecGlobalHNode[s+1].getOrigin()*weight
                                         + m_vecGlobalHNode[s].getOrigin()*weight
                                         + global_H_gravityC.getOrigin() );


            vel_gravityC.getVCenter() += m_vecGlobalVelNode[s+1].getVCenter()*weight+ m_vecGlobalVelNode[s].getVCenter()*weight;
            vel_gravityC.getVOrientation() += m_vecGlobalVelNode[s+1].getVOrientation()*weight+ m_vecGlobalVelNode[s].getVOrientation()*weight;

        }
        else
        {
            if (rigidification)
            {
                Vec3 pos_G = global_H_gravityC.getOrigin();
                global_H_gravityC.setOrientation(m_vecGlobalHNode[s].getOrientation() );
                global_H_gravityC.setOrigin(pos_G);

                SpatialVector VelNode_in_global, VelNode_in_Node, VelGravityC_inGravityC, VelGravity_inGlobal;
                VelNode_in_global.setAngularVelocity( m_vecGlobalVelNode[s].getVOrientation());
                VelNode_in_global.setLinearVelocity(  m_vecGlobalVelNode[s].getVCenter());

                // Projection of the velocity of the node in the frame of the node
                VelNode_in_Node.setLinearVelocity( m_vecGlobalHNode[s].backProjectVector(VelNode_in_global.getLinearVelocity()) );
                VelNode_in_Node.setAngularVelocity( m_vecGlobalHNode[s].backProjectVector(VelNode_in_global.getAngularVelocity() ) );

                // Transport of the velocity of the node to the gravity center (rigid link)
                Transform gravityC_H_Node = global_H_gravityC.inversed()*m_vecGlobalHNode[s];
                VelGravityC_inGravityC = gravityC_H_Node*VelNode_in_Node;

                // Projection of the velocity of the gravity center in the global frame
                VelGravity_inGlobal.setLinearVelocity(global_H_gravityC.projectVector( VelGravityC_inGravityC.getLinearVelocity() ) );
                VelGravity_inGlobal.setAngularVelocity(  global_H_gravityC.projectVector( VelGravityC_inGravityC.getAngularVelocity() ) );

                vel_gravityC.getVCenter()       = VelGravity_inGlobal.getLinearVelocity();
                vel_gravityC.getVOrientation()  = VelGravity_inGlobal.getAngularVelocity();

                m_vecGlobalHGravityCenter.push_back(global_H_gravityC);
                vec_Vel_gravityCenter.push_back(vel_gravityC);
                rigidification=false; // end of the rigidification

                global_H_gravityC.clear();
            }
        }
    }
    if (rigidification) // rigidification at the end tip
    {

        unsigned s=m_rigidBeamList.size()-1;

        Vec3 pos_G = global_H_gravityC.getOrigin();
        global_H_gravityC.setOrientation(m_vecGlobalHNode[s+1].getOrientation() );
        global_H_gravityC.setOrigin(pos_G);


        /// Computation of the velocities (can be put in an other function)
        SpatialVector VelNode_in_global, VelNode_in_Node, VelGravityC_inGravityC, VelGravity_inGlobal;
        VelNode_in_global.setAngularVelocity( m_vecGlobalVelNode[s].getVOrientation());
        VelNode_in_global.setLinearVelocity(  m_vecGlobalVelNode[s].getVCenter());

        /// Projection of the velocity of the node in the frame of the node
        VelNode_in_Node.setLinearVelocity( m_vecGlobalHNode[s].backProjectVector(VelNode_in_global.getLinearVelocity()) );
        VelNode_in_Node.setAngularVelocity( m_vecGlobalHNode[s].backProjectVector(VelNode_in_global.getAngularVelocity() ) );

        /// Transport of the velocity of the node to the gravity center (rigid link)
        Transform gravityC_H_Node = global_H_gravityC.inversed()*m_vecGlobalHNode[s];
        VelGravityC_inGravityC = gravityC_H_Node*VelNode_in_Node;

        /// projection of the velocity of the gravity center in the global frame
        VelGravity_inGlobal.setLinearVelocity(global_H_gravityC.projectVector( VelGravityC_inGravityC.getLinearVelocity() ) );
        VelGravity_inGlobal.setAngularVelocity(  global_H_gravityC.projectVector( VelGravityC_inGravityC.getAngularVelocity() ) );

        vel_gravityC.getVCenter()       = VelGravity_inGlobal.getLinearVelocity();
        vel_gravityC.getVOrientation()  = VelGravity_inGlobal.getAngularVelocity();

        m_vecGlobalHGravityCenter.push_back(global_H_gravityC);
        vec_Vel_gravityCenter.push_back(vel_gravityC);
        rigidification=false; // end of the rigidification
    }


    /// Set the beam and the topology
    recreateTopology();
    unsigned int numNodes=0;
    l_adaptiveInterpolation->clear();

    numNodes++;
    Coord xDof;
    xDof.getCenter()     = m_vecGlobalHNode[0].getOrigin();
    xDof.getOrientation()= m_vecGlobalHNode[0].getOrientation();
    x.push_back(xDof);
    v.push_back(m_vecGlobalVelNode[0]);
    m_topology->addPoint( xDof[0], xDof[1], xDof[2] );
    Real L;

    Rseg=0;
    rigidification=false;
    for (unsigned int s=0; s<m_rigidBeamList.size();s++)
    {
        if(m_rigidBeamList[s]) // the beam is rigidified
        {

            if (!rigidification) //begining of the rigidification
            {
                // the last element of vector x is replaced by the gravity center of the rigid zone
                x.pop_back();
                xDof.getCenter()     = m_vecGlobalHGravityCenter[Rseg].getOrigin();
                xDof.getOrientation()= m_vecGlobalHGravityCenter[Rseg].getOrientation();
                x.push_back(xDof);
                v.pop_back();
                v.push_back( vec_Vel_gravityCenter[Rseg] );
                Rseg++;
                rigidification=true;

                // add a transformation between the node 1 and the gravity center on previous beam
                if(s)
                {
                    Transform GravityCenter_H_node1 = m_vecGlobalHGravityCenter[Rseg-1].inversed()*m_vecGlobalHNode[s];
                    l_adaptiveInterpolation->setTransformBetweenDofAndNode(s-1,GravityCenter_H_node1,1);
                }

            }

            Transform GravityCenter_H_interpol0 = m_vecGlobalHGravityCenter[Rseg-1].inversed()*m_vecGlobalHNode[s];
            Transform GravityCenter_H_interpol1 = m_vecGlobalHGravityCenter[Rseg-1].inversed()*m_vecGlobalHNode[s+1];

            /// Add a beam that is on a segment that links the same dof
            m_topology->addEdge( (int)(numNodes-1), (int)(numNodes-1) );

            L =  newCurvAbs[s+1] - newCurvAbs[s];
            l_adaptiveInterpolation->addBeam(s, L, newCurvAbs[s], newCurvAbs[s+1] ,GravityCenter_H_interpol0, GravityCenter_H_interpol1 );
        }

        else // the beam is deformable
        {
            L =  newCurvAbs[s+1] - newCurvAbs[s];

            // ADD the beam in the topology and in the interpolation
            m_topology->addEdge( (int)(numNodes-1), (int)(numNodes) );
            l_adaptiveInterpolation->addBeam(s, L, newCurvAbs[s], newCurvAbs[s+1] ,0.0 );

            // ADD a DOF for the second node of the beam
            numNodes++;
            xDof.getCenter()     = m_vecGlobalHNode[s+1].getOrigin();
            xDof.getOrientation()= m_vecGlobalHNode[s+1].getOrientation();
            x.push_back(xDof);

            v.push_back(m_vecGlobalVelNode[s+1]);

            m_topology->addPoint( xDof[0], xDof[1], xDof[2] );


            if (rigidification) // end of the rigidification
            {
                rigidification=false;

                // add a transformation between the node 0 and the gravity center on current beam
                Transform GravityCenter_H_node0 = m_vecGlobalHGravityCenter[Rseg-1].inversed()*m_vecGlobalHNode[s];
                l_adaptiveInterpolation->setTransformBetweenDofAndNode(s,GravityCenter_H_node0,0);
            }

        }
    }

    getMechanicalState()->resize(x.size());
}

template <class DataTypes>
bool SutureController<DataTypes>::verifyRigidCurveSegmentSort()
{

    if(m_rigidCurveSegments.size()==0)
        return true;

    for (unsigned int seg=0; seg<m_rigidCurveSegments.size()-1; seg++)
    {
        if(m_rigidCurveSegments[seg].second > m_rigidCurveSegments[seg+1].first)
            return false;
    }
    return true;

}


template <class DataTypes>
void SutureController<DataTypes>::computeSampling(type::vector<Real> &newCurvAbs, VecCoord &x)
{
    type::vector<Real> xP_noticeable;
    type::vector<int> nbP_density;

    l_adaptiveInterpolation->getSamplingParameters(xP_noticeable, nbP_density);

    helper::ReadAccessor< Data< type::vector<Real> > > actualNoticeable = d_actualStepNoticeablePoints;
    if (!actualNoticeable.empty()) {
        xP_noticeable.insert(xP_noticeable.end(), actualNoticeable.begin(), actualNoticeable.end());
        std::sort(xP_noticeable.begin(), xP_noticeable.end());

        int density = nbP_density[0];
        nbP_density.clear();
        Real beamLength = l_adaptiveInterpolation->getRestTotalLength();
        for (size_t i = 1; i < xP_noticeable.size(); i++) {
            Real diff = xP_noticeable[i] - xP_noticeable[i-1];
            nbP_density.push_back(int((diff/beamLength)*Real(density)+1));
        }

        helper::WriteAccessor< Data< type::vector<Real> > > wActualNoticeable = d_actualStepNoticeablePoints;
        wActualNoticeable.clear();
    }

    if(xP_noticeable.size()<2){
        dmsg_error() <<" xP_noticeable_buf.size()= "<<xP_noticeable.size();
        return;
    }

    std::vector<Real> beamsCurvature;
    unsigned int nbBeams = l_adaptiveInterpolation->getNumBeams();
    beamsCurvature.resize(nbBeams);

    helper::WriteAccessor< Data< type::vector<Vec2> > > curvatureList = d_curvatureList;
    curvatureList.clear();
    curvatureList.resize(nbBeams);
    // Computing the curvature of each beam (from the previous timestep)
    for(unsigned int b=0; b<nbBeams; ++b)
    {
        Real beamLength = l_adaptiveInterpolation->getLength(b);
        l_adaptiveInterpolation->getAbsCurvXFromBeam(b, curvatureList[b][0]);
        Transform Tnode0, Tnode1;
        l_adaptiveInterpolation->computeTransform(b, Tnode0, Tnode1, x);

        curvatureList[b][1] = beamsCurvature[b] = l_adaptiveInterpolation->ComputeTotalBendingRotationAngle(beamLength / 5, Tnode0, Tnode1, beamLength, 0.0, 1.0);
    }

    type::vector<Real> newCurvAbs_notSecure;
    newCurvAbs_notSecure.clear();
    Real currentCurvAbs = 0.0, currentAngle = 0.0, maxAngle = d_maxBendingAngle.getValue();
    unsigned int currentBeam = 0;
    for(unsigned int part=0; part<nbP_density.size(); part++)
    {
        Real maxBeamLength = (xP_noticeable[part+1] - xP_noticeable[part]) / nbP_density[part];

        while(currentCurvAbs < xP_noticeable[part+1]-d_threshold.getValue())
        {
            newCurvAbs_notSecure.push_back(currentCurvAbs);
            Real maxCurvAbs = std::min(currentCurvAbs + maxBeamLength, xP_noticeable[part+1]);

            while(currentBeam < nbBeams)
            {
                Real beamStart, beamEnd, beamLength;
                l_adaptiveInterpolation->getAbsCurvXFromBeam(currentBeam, beamStart, beamEnd);
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

    if(!verifyRigidCurveSegmentSort())
        dmsg_warning()<< "RigidCurveSegments are not correctly sorted.";

    m_rigidCurveSegments.clear();
    helper::ReadAccessor< Data< set< Real > > > rigidCurvAbs = d_rigidCurvAbs;
    const auto nb = rigidCurvAbs->size();

    if(nb>0 && (nb%2)==0) // Make sure we have pairs of curv abs
    {
        RealConstIterator it;
        for(it=rigidCurvAbs->begin(); it!=rigidCurvAbs->end();)
        {
            Real start, end;
            start = *it++;
            end = *it++;
            m_rigidCurveSegments.push_back(std::make_pair(start, end));
        }
        addRigidCurvAbs(newCurvAbs_notSecure, 0.0001);
    }

    /// Verify that there are no beams with null length
    /// DEBUG: should not remove the last abscissa (spoils the object)
    newCurvAbs.clear();
    newCurvAbs.push_back(newCurvAbs_notSecure[0]);
    for (unsigned int i=1; i<newCurvAbs_notSecure.size()-1; i++)
    {
        if (newCurvAbs_notSecure[i] > newCurvAbs_notSecure[i-1]+d_threshold.getValue())
            newCurvAbs.push_back(newCurvAbs_notSecure[i]);
    }
    size_t lastSec = newCurvAbs.size()-1;
    size_t lastNSec = newCurvAbs_notSecure.size()-1;
    if (newCurvAbs_notSecure[lastNSec] < newCurvAbs[lastSec]+d_threshold.getValue()) {
        newCurvAbs.pop_back();
    }
    newCurvAbs.push_back(newCurvAbs_notSecure[lastNSec]);
}

template <class DataTypes>
void SutureController<DataTypes>::verifyRigidSegmentsSampling(type::vector<Real> &newCurvAbs)
{
    // Making sure we keep the same sampling in the rigid segments from one timestep to the next
    if(!d_fixRigidTransforms.getValue())
        return;

    const type::vector<Real> &oldCurvAbs = d_nodeCurvAbs.getValue();
    typename type::vector<Real>::iterator newIter, newIter2;
    typename type::vector<Real>::const_iterator oldIter, oldIter2;
    newIter = newCurvAbs.begin();
    oldIter = oldCurvAbs.begin();

    typename type::vector< std::pair<Real, Real> >::const_iterator rigidIter;
    // For each segment
    for(rigidIter = m_rigidCurveSegments.begin(); rigidIter != m_rigidCurveSegments.end(); ++rigidIter)
    {
        if(std::find(m_prevRigidCurvSegments.begin(), m_prevRigidCurvSegments.end(), *rigidIter) == m_prevRigidCurvSegments.end())
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
    if(!d_fixRigidTransforms.getValue())
        return;

    m_prevRigidTransforms.clear();

    typename type::vector< std::pair<Real, Real> >::const_iterator rigidIter;
    // For each rigid segment
    for(rigidIter = m_prevRigidCurvSegments.begin(); rigidIter != m_prevRigidCurvSegments.end(); ++rigidIter)
    {
        double rigidStart = rigidIter->first, rigidEnd = rigidIter->second;

        // For all curv abs in this segment
        unsigned int beamId, lastBeamId;
        Real bary;
        l_adaptiveInterpolation->getBeamAtCurvAbs(rigidStart, beamId, bary);
        l_adaptiveInterpolation->getBeamAtCurvAbs(rigidEnd, lastBeamId, bary);

        while(beamId < lastBeamId)
        {
            // Save the transformation
            //  (we consider that we don't cut rigid beams so transformations are the same for 2 beams in the same curv abs)
            Real curvAbs0, curvAbs1;
            l_adaptiveInterpolation->getAbsCurvXFromBeam(beamId, curvAbs0, curvAbs1);

            Transform t0, t1;
            l_adaptiveInterpolation->getDOFtoLocalTransform(beamId, t0, t1);
            m_prevRigidTransforms[curvAbs0] = t0;
            m_prevRigidTransforms[curvAbs1] = t1;

            ++beamId;
        }
    }
}

template <class DataTypes>
void SutureController<DataTypes>::verifyRigidSegmentsTransformations()
{
    if(!d_fixRigidTransforms.getValue())
        return;

    typename type::vector< std::pair<Real, Real> >::const_iterator rigidIter;
    // For each rigid segment
    for(rigidIter = m_rigidCurveSegments.begin(); rigidIter != m_rigidCurveSegments.end(); ++rigidIter)
    {
        if(std::find(m_prevRigidCurvSegments.begin(), m_prevRigidCurvSegments.end(), *rigidIter) == m_prevRigidCurvSegments.end())
            continue;	// If this is a new segment, don't modify it

        Real rigidStart = rigidIter->first, rigidEnd = rigidIter->second;

        // For all curv abs in this segment
        unsigned int beamId, lastBeamId;
        Real bary;
        l_adaptiveInterpolation->getBeamAtCurvAbs(rigidStart, beamId, bary);
        l_adaptiveInterpolation->getBeamAtCurvAbs(rigidEnd, lastBeamId, bary);

        while(beamId < lastBeamId)
        {
            // Load the transformations and replace what was computed this timestep
            Real curvAbs0, curvAbs1;
            l_adaptiveInterpolation->getAbsCurvXFromBeam(beamId, curvAbs0, curvAbs1);

            typename std::map<Real, Transform>::const_iterator iter;
            iter = m_prevRigidTransforms.find(curvAbs0);
            if(iter != m_prevRigidTransforms.end())
                l_adaptiveInterpolation->setTransformBetweenDofAndNode(beamId, iter->second, false);

            iter = m_prevRigidTransforms.find(curvAbs1);
            if(iter != m_prevRigidTransforms.end())
                l_adaptiveInterpolation->setTransformBetweenDofAndNode(beamId, iter->second, true);

            ++beamId;
        }
    }
}

template <class DataTypes>
void SutureController<DataTypes>::updateControlPointsPositions()
{
    auto ctrlPts = sofa::helper::getWriteOnlyAccessor(d_controlPoints);
    ctrlPts.clear();

    unsigned int numBeams = l_adaptiveInterpolation->getNumBeams();
    Transform global_H0_local, global_H1_local;
    const VecCoord& x = getMechanicalState()->write(sofa::core::VecCoordId::position())->getValue();
    for (unsigned int b = 0; b < numBeams; b++)
    {
        l_adaptiveInterpolation->computeTransform(b, global_H0_local, global_H1_local, x);
        Coord pt;
        pt.getCenter() = global_H0_local.getOrigin();
        pt.getOrientation() = global_H0_local.getOrientation();
        ctrlPts.push_back(pt);
    }
    Coord pt;
    pt.getCenter() = global_H1_local.getOrigin();
    pt.getOrientation() = global_H1_local.getOrientation();
    ctrlPts.push_back(pt);
}

template <class DataTypes>
void SutureController<DataTypes>::insertActualNoticeablePoint(Real _absc)
{
    helper::WriteAccessor< Data< type::vector<Real> > > actualNoticeable = d_actualStepNoticeablePoints;
    actualNoticeable.push_back(_absc);
}

template <class DataTypes>
void SutureController<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowBehaviorModels()) return;

    if (m_rigidCurveSegments.size() != m_vecGlobalHGravityCenter.size())
    {
        dmsg_error() <<"In draw function rigidCurveSegments.size() ="<< m_rigidCurveSegments.size() <<" != vec_global_H_gravityCenter.size() = "<<m_vecGlobalHGravityCenter.size();
    }

    for (unsigned int i=0; i<m_vecGlobalHGravityCenter.size(); i++)
    {
        Real Length = m_rigidCurveSegments[i].second - m_rigidCurveSegments[i].first;
        Vec3 sizeArrows (Length/4, Length/8, Length/8);

        vparams->drawTool()->drawFrame(m_vecGlobalHGravityCenter[i].getOrigin(), m_vecGlobalHGravityCenter[i].getOrientation(), sizeArrows);
    }
}

} // namespace sofa::component::controller::_suturecontroller_

