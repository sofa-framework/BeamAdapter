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
// C++ Implementation : InterventionalRadiology
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGY_INL
#define SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGY_INL

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <SofaBaseTopology/EdgeSetGeometryAlgorithms.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/system/FileRepository.h>

#include "InterventionalRadiology.h"

namespace sofa
{

namespace component
{

namespace controller
{

namespace _interventionalradiology_
{

using helper::vector;
using core::objectmodel::BaseContext;
using helper::WriteAccessor;
using core::objectmodel::KeypressedEvent;
using core::objectmodel::MouseEvent;


template <class DataTypes>
InterventionalRadiology<DataTypes>::InterventionalRadiology()
    : d_instrumentsPath(initData(&d_instrumentsPath,"instruments", "List of paths to WireInterpolation components on the scene"))
    , d_controlledInstrument(initData(&d_controlledInstrument, 0, "controlledInstrument", "provide the id of the interventional radiology instrument which is under control: press contr + number to change it"))
    , d_xTip(initData(&d_xTip,"xtip", "curvilinear abscissa of the tip of each interventional radiology instrument, this has to be smaller than the lenght of one beam"))
    , d_rotationInstrument(initData(&d_rotationInstrument,"rotationInstrument", "angle of rotation for each interventional radiology instrument"))
    , d_step(initData(&d_step,(Real)0.1,"step","base step when changing beam length"))
    , d_angularStep(initData(&d_angularStep,(Real)(3.1416/20.0),"angularStep","base step when changing beam angle"))
    , d_speed(initData(&d_speed,(Real)0.0,"speed","continuous beam length increase/decrease"))
    , d_startingPos(initData(&d_startingPos,Coord(),"startingPos","starting pos for inserting the instrument"))
    , d_threshold(initData(&d_threshold, (Real)0.01, "threshold", "threshold for controller precision which is homogeneous to the unit of length used in the simulation"))
    //    , d_rigidCurvAbs(initData(&d_rigidCurvAbs, "rigidCurvAbs", "pairs of curv abs for beams we want to rigidify"))
    , d_indexFirstNode(initData(&d_indexFirstNode, (unsigned int) 0, "indexFirstNode", "first node (should be fixed with restshape)"))
    , d_curvAbs(initData(&d_curvAbs,"CurvAbs", "curvi-linear abscissa of the DOFs" ))
    , d_totalLengths(initData(&d_totalLengths,"totalLengths", "the total length of wireRestShaps of each instrument" ))
{
    m_fixedConstraint = NULL;

    // the controller must listen to the event (in particular BeginAnimationStep event)
    if (!f_listening.isSet()) f_listening.setValue(true);
}


template <class DataTypes>
void InterventionalRadiology<DataTypes>::init()
{
    BaseContext* context = getContext();

    //get the pointers of the WireBeamInterpolations
    const vector<std::string>& instrumentPathList = d_instrumentsPath.getValue();
    if (!instrumentPathList.empty())
    {
        for (unsigned int i=0;i<instrumentPathList.size();++i)
        {
            WBeamInterpolation * wbinterpol = NULL;
            context->get(wbinterpol,instrumentPathList[i]);
            if( wbinterpol)
                m_instrumentsList.push_back(wbinterpol)  ;
            else
                msg_error() << "Interpolation of instrument "<<instrumentPathList[i]<< "  not found.";
        }
    }
    if(m_instrumentsList.empty() || m_instrumentsList.size() == 0)
    {        
        msg_error() << "The list of interpolationBeam is not found.";
        //WBeamInterpolation * wbinterpol= context->get< WBeamInterpolation >(BaseContext::Local);
        //m_instrumentsList.push_back(wbinterpol);
        msg_error()<<"No instrument found ( no WireBeamInterpolation) or ";
        msg_error()<<"No Beam Interpolation found !!! the component can not work.";
        return;
    }

    m_activatedPointsBuf.clear();

    if((d_xTip.getValue().size()!= m_instrumentsList.size()) || (d_rotationInstrument.getValue().size()!= m_instrumentsList.size())){
        msg_error()<<"Please, check the dimention of rotationInstrument, x_ti and m_instrumentsList ";
        return;
    }

    for(unsigned int i=0; i<m_instrumentsList.size(); i++)
        m_instrumentsList[i]->setControlled(true);

    context->get(m_fixedConstraint);
    if(m_fixedConstraint==NULL)
        msg_error()<<"No fixedConstraint found.";

    m_nodeCurvAbs.clear();
    m_nodeCurvAbs.push_back(0.0);

    m_idInstrumentCurvAbsTable.clear();
    vector<int> listInit;

    for(unsigned int i=0; i<m_instrumentsList.size(); i++)
        listInit.push_back(int(i));

    m_idInstrumentCurvAbsTable.push_back(listInit);

//    vector<Real>& totalLengths =  *d_totalLengths.beginEdit();
    auto totalLengths = getWriteOnlyAccessor(d_totalLengths);
    totalLengths.clear();

    for (unsigned int i = 0; i < m_instrumentsList.size();i++) {
        Real lengh = m_instrumentsList[i]->getRestTotalLength();
        totalLengths.push_back(lengh);
    }
    //    d_totalLengths.endEdit();
    //std::cout<< "IIIIIIIIIIIIIIIIII=========================> "<< d_totalLengths.getValue()<< std::endl;

    Inherit::init();
}


template <class DataTypes>
void InterventionalRadiology<DataTypes>::bwdInit()
{
    // assign the starting pos to each point of the Mechanical State
    Coord stPos =d_startingPos.getValue();
    stPos.getOrientation().normalize();
    d_startingPos.setValue(stPos);

    WriteAccessor<Data<VecCoord> > x = *getMechanicalState()->write(core::VecCoordId::position());
    msg_info() << " =========> x.size :"<< x.size();
    for(unsigned int i=0; i<x.size(); i++)
        x[i] = d_startingPos.getValue();
    m_numControlledNodes = x.size();


    //verify  if the totalLength is 0, move the first instrument
    Real totalLengthCombined = 0.0;
    auto xTip = getWriteAccessor(d_xTip);
    for (unsigned int m = 0; m < xTip.size(); ++m) {
        if (xTip[m] > totalLengthCombined)
            totalLengthCombined = xTip[m];
    }
    if(totalLengthCombined<0.0001)
        xTip[0]=0.0001;
    applyInterventionalRadiology();


}

template <class DataTypes>
void InterventionalRadiology<DataTypes>::reinit()
{
}


template <class DataTypes>
void InterventionalRadiology<DataTypes>::onBeginAnimationStep(const double dt)
{
    //Check if update is necessary
    if (d_old_xTip == d_xTip.getValue() && d_old_rotationInstrument == d_rotationInstrument.getValue()) return;

    //Perform Update
    applyInterventionalRadiology();

    //Store new values
    d_old_xTip = d_xTip.getValue();
    d_old_rotationInstrument = d_rotationInstrument.getValue();
}


template <class DataTypes>
void InterventionalRadiology<DataTypes>::interventionalRadiologyComputeSampling(vector<Real> &newCurvAbs,
                                                                                vector< vector<int> > &idInstrumentTable,
                                                                                const vector<Real> &xBegin)
{
    // Step 1 = put the noticeable Nodes (keyPoints from wireRestShape)
    double maxAbsLength=0.0;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        vector<Real> xP_noticeable_I;
        vector< int > density_I;
        m_instrumentsList[i]->getSamplingParameters(xP_noticeable_I, density_I);

        for (unsigned int j=0; j<xP_noticeable_I.size(); j++)
        {
            //compute the corresponding abs curv of this "noticeable point" on the combined intrument
            Real curvAbs_xP = xBegin[i] + xP_noticeable_I[j];
            if (curvAbs_xP>0.0)   // all the noticiable point that have a negative curv abs are not simulated => considered as outside of the patient...
            {
                newCurvAbs.push_back(curvAbs_xP);

                if (curvAbs_xP > maxAbsLength)
                    maxAbsLength=curvAbs_xP;
            }
        }
    }


    // Step 2 => add the beams given the sampling parameters
    Real xSampling = 0.0;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        vector<Real> xPNoticeableI;
        vector< int > density_I;
        m_instrumentsList[i]->getSamplingParameters(xPNoticeableI, density_I);

        for (unsigned int j=0; j<density_I.size(); j++){

            //compute the corresponding abs curv of this "noticeable point" on the combined intrument
            Real curvAbsxP = xBegin[i] + xPNoticeableI[j+1];

            // use density parameter (size = xP_noticeable_I -1 )
            if (curvAbsxP > xSampling && density_I[j]>0)
            {
                Real ratio = (Real)density_I[j] / (xPNoticeableI[j+1]  - xPNoticeableI[j]) ;
                int  numNewNodes = (int)floor( (curvAbsxP- xSampling)  *ratio) ;

                for (int k=0; k<numNewNodes; k++)
                {
                    newCurvAbs.push_back( xPNoticeableI[j+1] + xBegin[i]  - (k+1) * (1/ratio) );
                }
                xSampling = curvAbsxP;
            }
        }
    }

    sortCurvAbs(newCurvAbs, idInstrumentTable);
}


template <class DataTypes>
void InterventionalRadiology<DataTypes>::activateBeamListForCollision( vector< vector<int> > &idInstrumentTable)
{
    // 0. useful for rigidification
    //    const vector<Real>  *rigidCurvAbs = &d_rigidCurvAbs.getValue();

    // 1. clear the information related to the collision for each instrument
    for (unsigned int instr=0; instr<m_instrumentsList.size(); instr++)
        m_instrumentsList[instr]->clearCollisionOnBeam();

    // 2.   for each node along the structure, detect the "visible" instrument (the one with the largest section is supposed to be the first in the list)
    //      if visible, the beam is assigned for collision on this instrument
    for (unsigned int p=0; p<idInstrumentTable.size()-1; p++)
    {
        unsigned int instr0;
        if(idInstrumentTable.size()==1 )
        {
            instr0=idInstrumentTable[0][0];
            m_instrumentsList[ instr0 ]->addCollisionOnBeam(p);
        }
    }
}


template <class DataTypes>
void InterventionalRadiology<DataTypes>::applyInterventionalRadiology()
{
    /// Create vectors with the CurvAbs of the noticiable points and the id of the corresponding instrument


    /// STEP 1
    /// Find the total length of the COMBINED INSTRUMENTS and the one for which xtip > 0 (so the one which are simulated)
    ///

    helper::AdvancedTimer::stepBegin("step1");
    vector<Real> newCurvAbs;
    vector<Real> xbegin;
    Real totalLengthCombined=0.0;
    getTotalLengthCombined(newCurvAbs, xbegin, totalLengthCombined);
    helper::AdvancedTimer::stepEnd("step1");

    /// STEP 2:
    /// get the noticeable points that need to be simulated
    // Fill=> newCurvAbs which provides a vector with curvilinear abscissa of each simulated node
    //     => id_instrument_table which provides for each simulated node, the id of all instruments which belong this node
    //     => xbegin (theoritical curv abs of the beginning point of the instrument (could be negative) xbegin= xtip - intrumentLength)
    helper::AdvancedTimer::stepBegin("step2");
    vector<vector<int>> idInstrumentTable;
    interventionalRadiologyComputeSampling(newCurvAbs,idInstrumentTable, xbegin);
    helper::AdvancedTimer::stepEnd("step2");


    /// STEP 3
    /// Re-interpolate the positions and the velocities
    helper::AdvancedTimer::stepBegin("step3");
    unsigned int nnode=newCurvAbs.size(); // number of simulated nodes
    unsigned int nbeam=nnode - 1; // number of simulated beams

    unsigned int nnode_old= m_nodeCurvAbs.size(); // previous number of simulated nodes;
    auto datax = this->getMechanicalState()->write(core::VecCoordId::position());
    //    VecCoord& x = *datax->beginEdit();
    auto x = sofa::helper::write(*datax);    /// WriteAccessor

    VecCoord xbuf ;
    xbuf = x.ref();    ///

    vector<Real> modifiedCurvAbs;
    totalLengthIsChanging(newCurvAbs, modifiedCurvAbs, idInstrumentTable);

    Real xmax_prev = m_nodeCurvAbs[m_nodeCurvAbs.size()-1];

    for (unsigned int p=0; p <= nbeam; p++)
    {
        int idP = m_numControlledNodes - nnode + p;
        Real xabs = modifiedCurvAbs[p];

        // 2 cases:  TODO : remove first case
        //1. the abs curv is further than the previous state of the instrument
        //2. this is not the case and the node position can be interpolated using previous step positions
        if(xabs > xmax_prev + d_threshold.getValue())
        {
            std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<< std::endl;
        }
        else
        {
            // case 2 (the node position can be interpolated straightfully using previous step positions)
            unsigned int p0=0;
            while(p0<m_nodeCurvAbs.size())
            {
                if((m_nodeCurvAbs[p0]+d_threshold.getValue())>xabs)
                    break;
                p0++;
            }

            unsigned int idP0 =  m_numControlledNodes - nnode_old + p0 ;

            msg_info() << "idP : "<< idP << " idP0 :"<< idP0;

            if(fabs(m_nodeCurvAbs[p0]-xabs)<d_threshold.getValue()){
                msg_info() << " ==++++> I am there! ";
                x[idP] = xbuf[idP0];
            }
            else
            {
                // the node must be interpolated using beam interpolation
                //find the instrument
                int id = m_idInstrumentCurvAbsTable[p0][0];
                //find the good beam (TODO: do not work if xbegin of one instrument >0)
                int b = static_cast<int>(p0) - 1;
                // test to avoid wrong indices
                if (b<0)
                    x[p]=d_startingPos.getValue();
                else
                {
                    Transform global_H_interpol;
                    Real ratio = (xabs - m_nodeCurvAbs[b])/ (m_nodeCurvAbs[p0]-m_nodeCurvAbs[b]);
                    Transform Global_H_local0(xbuf[idP0-1].getCenter(),xbuf[idP0-1].getOrientation() ),
                            Global_H_local1(xbuf[idP0].getCenter(),xbuf[idP0].getOrientation() );

                    Real L = m_nodeCurvAbs[p0] - m_nodeCurvAbs[b];

                    m_instrumentsList[id]->InterpolateTransformUsingSpline(global_H_interpol, ratio, Global_H_local0, Global_H_local1 ,L);

                    x[idP].getCenter() = global_H_interpol.getOrigin();
                    x[idP].getOrientation() = global_H_interpol.getOrientation();

                }
            }
        }
    }
    helper::AdvancedTimer::stepEnd("step3");


    /// STEP 4
    /// Assign the beams
    helper::AdvancedTimer::stepBegin("step4");
    unsigned int numEdges= m_numControlledNodes-1;

    // verify that there is a sufficient number of Edges in the topology
    if(numEdges<nbeam)
    {
        if (f_printLog.getValue())
            msg_error()<<"Not enough edges in the topology.";
        nbeam=numEdges;
    }


    for (unsigned int b=0; b<nbeam; b++)
    {
        Real x0 = newCurvAbs[b];
        Real x1 = newCurvAbs[b+1];
        for (unsigned int i=0; i<m_instrumentsList.size(); i++)
        {
            Real xmax = d_xTip.getValue()[i];
            Real xmin = xbegin[i];

            Real eps= d_threshold.getValue();

            if (x0>(xmin-eps) && x0<(xmax+eps) && x1>(xmin-eps) && x1<(xmax+eps))
            {
                BaseMeshTopology::EdgeID eID = (BaseMeshTopology::EdgeID)(numEdges-nbeam + b );

                Real length = x1 - x0;
                Real x0_local = x0-xmin;
                Real x1_local = x1-xmin;

                Real theta = d_rotationInstrument.getValue()[i];

                m_instrumentsList[i]->addBeam(eID, length, x0_local, x1_local, theta);

            }
        }
    }
    helper::AdvancedTimer::stepEnd("step4");

    /// STEP 5
    /// Fix the not simulated nodes
    helper::AdvancedTimer::stepBegin("step5");
    unsigned int firstSimulatedNode = m_numControlledNodes - nbeam;

    //1. Fix the nodes (beginning of the instruments) that are not "out"
    fixFirstNodesWithUntil(firstSimulatedNode);

    /// STEP 6
    /// Activate Beam List for collision of each instrument
    activateBeamListForCollision(idInstrumentTable);

    m_nodeCurvAbs = newCurvAbs;
    m_idInstrumentCurvAbsTable = idInstrumentTable;
    helper::AdvancedTimer::stepEnd("step5");
}

template <class DataTypes>
void InterventionalRadiology<DataTypes>::totalLengthIsChanging(const vector<Real>& newNodeCurvAbs,
                                                               vector<Real>& modifiedNodeCurvAbs,
                                                               const vector< vector<int> >& newTable)
{

    // If total length is changing, it means that we need to simulate the fact that some beam will "move" during the time step.
    // thus during the interpolation process, where the position of the nodes is based on the previous position of the wire,
    // we initialize some points at a x_curv ref pos without the motion (computed by DLength)
    // due to the elasticity of the beam, the point will then naturally go the position that reespects the newNodeCurvAbs

    Real dLength = newNodeCurvAbs[ newNodeCurvAbs.size()-1] - m_nodeCurvAbs[m_nodeCurvAbs.size() - 1];
    modifiedNodeCurvAbs = newNodeCurvAbs;

    // we look for the last value in the CurvAbs
    if(fabs(dLength) > d_threshold.getValue())
    {
        unsigned int i=newTable.size()-1;
        while (i>0 && newTable[i].size()==1)
        {
            modifiedNodeCurvAbs[i]-=dLength;

            // force modifiedNode to be "locally" sorted
            if(modifiedNodeCurvAbs[i]<modifiedNodeCurvAbs[i-1])
            {
                modifiedNodeCurvAbs[i] = modifiedNodeCurvAbs[i-1]+ d_threshold.getValue();
            }

            i--;
        }
    }
}

template <class DataTypes>
void InterventionalRadiology<DataTypes>::sortCurvAbs(vector<Real> &curvAbs,
                                                     vector<vector<int> >& idInstrumentTable)
{
    //    vector<Real> &newCurvAbs =(*d_curvAbs.beginEdit());
    auto newCurvAbs  = getWriteOnlyAccessor(d_curvAbs);
    Real eps = d_threshold.getValue();

    // here we sort CurvAbs
    // a buffer verctor: newCurvAbs will be filled by iteratively removing the minimal values found in CurvAbs and push_back them in newCurvAbs
    // a threshod is used to remove the values that are "too" similars...
    // TODO: could be improve by function "sort" ???

    newCurvAbs.clear();
    while(curvAbs.size()>0)
    {
        Real xMin = 1.0e30;
        unsigned int index_min=0;
        for(unsigned int j=0; j<curvAbs.size(); j++)
        {
            if(xMin > curvAbs[j] )
            {
                xMin = curvAbs[j];
                index_min=j;
            }
        }

        helper::removeIndex( curvAbs, index_min);

        // verify using a eps that the same node does not already exist
        Real xLast;
        if (newCurvAbs.size()>0 )
            xLast= newCurvAbs[newCurvAbs.size()-1];
        else
            xLast=-1.0e30;

        if( fabs(xMin - xLast) > eps)
        {
            newCurvAbs.push_back(xMin);
        }
    }

    curvAbs = newCurvAbs.ref();

    // here we build a table that provides the list of each instrument for each dof in the list of newCurvAbs
    // dofs can be shared by several instruments
    idInstrumentTable.clear();

    for (unsigned int i=0; i<newCurvAbs.size(); i++)
    {
        vector<int> listNewNode;
        for (unsigned int id=0; id<m_instrumentsList.size(); id++)
        {
            Real xEnd= d_xTip.getValue()[id];
            Real xBegin = xEnd - m_instrumentsList[id]->getRestTotalLength();

            if ( xBegin<(newCurvAbs[i]+eps) && xEnd >(newCurvAbs[i]-eps) )
                listNewNode.push_back(id);
        }

        if(listNewNode.size() ==0)
            msg_error()<<"No instrument find for curvAbs"<<newCurvAbs[i];

        idInstrumentTable.push_back(listNewNode);

    }
}


template <class DataTypes>
void InterventionalRadiology<DataTypes>::fixFirstNodesWithUntil(unsigned int firstSimulatedNode)
{
    WriteAccessor<Data<VecCoord> > xMstate = *getMechanicalState()->write(core::VecCoordId::position());
    WriteAccessor<Data<VecDeriv> > vMstate = *getMechanicalState()->write(core::VecDerivId::velocity());

    // set the position to startingPos for all the nodes that are not simulated
    // and add a fixedConstraint
    m_fixedConstraint->clearConstraints();
    for(unsigned int i=0; i<firstSimulatedNode-1 ; i++)
    {
        xMstate[i]=d_startingPos.getValue();
        vMstate[i].clear();
        m_fixedConstraint->addConstraint(i);
    }
    d_indexFirstNode = firstSimulatedNode-1 ;
}


template <class DataTypes>
std::string InterventionalRadiology<DataTypes>::getTemplateName() const
{
    return templateName(this);
}

template <class DataTypes>
std::string InterventionalRadiology<DataTypes>::templateName(const InterventionalRadiology<DataTypes>* thisClass)
{
    SOFA_UNUSED(thisClass);
    return DataTypes::Name();
}

//template <class DataTypes>
//bool InterventionalRadiology<DataTypes>::activePoint(int index, core::CollisionModel * cm)
//{
//    SOFA_UNUSED(cm);
//    return m_activatedPointsBuf[index];
//}

//template <class DataTypes>
//bool InterventionalRadiology<DataTypes>::activeLine(int index, core::CollisionModel * cm)
//{
//    SOFA_UNUSED(cm);
//    return m_activatedPointsBuf[index+1];
//}

template <class DataTypes>
bool InterventionalRadiology<DataTypes>::modifyTopology(void)
{
    return false;
}

template <class DataTypes>
void InterventionalRadiology<DataTypes>::getInstrumentList(vector<fem::WireBeamInterpolation<DataTypes>*>& list)
{
    list = m_instrumentsList;
}

template <class DataTypes>
const vector< vector<int> >& InterventionalRadiology<DataTypes>::get_id_instrument_curvAbs_table() const
{
    return m_idInstrumentCurvAbsTable;
}

template <class DataTypes>
int InterventionalRadiology<DataTypes>::getTotalNbEdges() const
{
    return getContext()->getMeshTopology()->getNbEdges();
}


} // namespace _interventionalradiology_

} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGY_INL */
