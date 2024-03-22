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
// C++ Implementation : InterventionalRadiologyController
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

#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/component/topology/container/dynamic/EdgeSetGeometryAlgorithms.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/AdvancedTimer.h>
#include <sofa/helper/system/FileRepository.h>
#include <sofa/type/vector_algorithm.h>


#include <BeamAdapter/component/controller/InterventionalRadiologyController.h>


namespace sofa::component::controller::_interventionalradiologycontroller_
{

using type::vector;
using core::objectmodel::BaseContext;
using helper::WriteAccessor;
using core::objectmodel::KeypressedEvent;
using core::objectmodel::MouseEvent;
using namespace sofa::beamadapter;


template <class DataTypes>
InterventionalRadiologyController<DataTypes>::InterventionalRadiologyController()
: d_instrumentsPath(initData(&d_instrumentsPath,"instruments", "List of paths to WireInterpolation components on the scene"))
, d_controlledInstrument(initData(&d_controlledInstrument, 0, "controlledInstrument", "provide the id of the interventional radiology instrument which is under control: press contr + number to change it"))
, d_xTip(initData(&d_xTip,"xtip", "curvilinear abscissa of the tip of each interventional radiology instrument"))
, d_rotationInstrument(initData(&d_rotationInstrument,"rotationInstrument", "angle of rotation for each interventional radiology instrument"))
, d_step(initData(&d_step,(Real)0.1,"step","base step when changing beam length"))
, d_angularStep(initData(&d_angularStep,(Real)(3.1416/20.0),"angularStep","base step when changing beam angle"))
, d_speed(initData(&d_speed,(Real)0.0,"speed","continuous beam length increase/decrease"))
, d_startingPos(initData(&d_startingPos,Coord(),"startingPos","starting pos for inserting the instrument"))
, d_threshold(initData(&d_threshold, (Real)0.01, "threshold", "threshold for controller precision which is homogeneous to the unit of length used in the simulation"))
, d_rigidCurvAbs(initData(&d_rigidCurvAbs, "rigidCurvAbs", "pairs of curv abs for beams we want to rigidify"))
, d_motionFilename(initData(&d_motionFilename, "motionFilename", "text file that includes tracked motion from optical sensor"))
, d_indexFirstNode(initData(&d_indexFirstNode, (unsigned int) 0, "indexFirstNode", "first node (should be fixed with restshape)"))
{
    m_fixedConstraint = nullptr;
    m_sensored =false;
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::init()
{
    BaseContext* context = getContext();
    this->mState = nullptr;

    //get the pointers of the WireBeamInterpolations
    const type::vector<std::string>& instrumentPathList = d_instrumentsPath.getValue();
    if (instrumentPathList.empty())
    {
        getContext()->template get<WBeamInterpolation>(&m_instrumentsList, sofa::core::objectmodel::BaseContext::Local);
    }
    else
    {
        for (unsigned int i=0;i<instrumentPathList.size();++i)
        {
            WBeamInterpolation * wbinterpol = nullptr;
            context->get(wbinterpol,instrumentPathList[i]);
            if( wbinterpol)
                m_instrumentsList.push_back(wbinterpol)  ;
            else
                msg_error() << "Interpolation of instrument "<<instrumentPathList[i]<< "  not found.";
        }
    }

    if (m_instrumentsList.empty()) {
        msg_error() << "No instrument found (no WireBeamInterpolation)! the component can not work and will be set to Invalid.";
        sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }
    else
    {
        msg_info() << m_instrumentsList.size() << " instrument(s) found (WireBeamInterpolation)";
    }

     m_activatedPointsBuf.clear();

     if(d_speed.getValue()>0)
     {
         m_FF=true;
         m_RW=false;
         m_sensored = false;
     }
    if (!d_motionFilename.getValue().empty())
    {
        m_FF = true; m_sensored = true; m_currentSensorData = 0;
        loadMotionData(d_motionFilename.getValue());
    }

    auto x_instr_tip = sofa::helper::getWriteOnlyAccessor(d_xTip);
    x_instr_tip.resize(m_instrumentsList.size());

    auto angle_Instrument = sofa::helper::getWriteOnlyAccessor(d_rotationInstrument);
    angle_Instrument.resize(m_instrumentsList.size());

    for(unsigned int i=0; i<m_instrumentsList.size(); i++)
        m_instrumentsList[i]->setControlled(true);

    context->get(m_fixedConstraint);
    if(m_fixedConstraint==nullptr)
        msg_error()<<"No fixedConstraint found.";

    // the controller must listen to the event (in particular BeginAnimationStep event)
    if (!f_listening.isSet())
    {
        f_listening.setValue(true);
    }

    m_nodeCurvAbs.clear();
    m_idInstrumentCurvAbsTable.clear();
    m_nodeCurvAbs.push_back(0.0);
    type::vector<int> listInit;

    for(unsigned int i=0; i<m_instrumentsList.size(); i++)
        listInit.push_back(int(i));

    m_idInstrumentCurvAbsTable.push_back(listInit);

    Inherit::init();

    sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}

template<class DataTypes>
void InterventionalRadiologyController<DataTypes>::loadMotionData(std::string filename)
{
    if (!helper::system::DataRepository.findFile(filename))
    {
        msg_error() << "File " << filename << " not found.";
        return;
    }
    std::ifstream file(filename.c_str());

    std::string line;
    Vec3 result;
    while( std::getline(file,line) )
    {
        if (line.empty())
            continue;
        std::istringstream values(line);
        values >> result[0] >> result[1] >> result[2];
        result[0] /= 1000;
        m_sensorMotionData.push_back(result);
    }

    file.close();
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::bwdInit()
{
    // assign the starting pos to each point of the Mechanical State
    Coord stPos =d_startingPos.getValue();
    stPos.getOrientation().normalize();
    d_startingPos.setValue(stPos);

    if (!this->mState) {
        msg_error() << "No MechanicalState found. The component can not work and will be set to Invalid.";
        sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        return;
    }

    WriteAccessor<Data<VecCoord> > x = *this->mState->write(core::VecCoordId::position());
    for(unsigned int i=0; i<x.size(); i++)
        x[i] = d_startingPos.getValue();
    m_numControlledNodes = x.size();

    sofa::Size nbrBeam = 0;
    for (unsigned int i = 0; i < m_instrumentsList.size(); i++)
    {
        type::vector<Real> xP_noticeable_I;
        type::vector< int > density_I;
        m_instrumentsList[i]->getSamplingParameters(xP_noticeable_I, density_I);

        for (auto nb : density_I)
            nbrBeam += nb;
    }

    if (nbrBeam > m_numControlledNodes)
    {
        msg_warning() << "Parameter missmatch: According to the list of controlled instrument controlled. The number of potential beams: "
            << nbrBeam << " exceed the number of degree of freedom in the MechanicalObject: " << m_numControlledNodes << ". This could lead to unespected behavior.";
    }
        
    applyInterventionalRadiologyController();

    sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
}


/*!
 * \todo fix the mouse event with better controls
 */
template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::onMouseEvent(MouseEvent * mev)
{
    SOFA_UNUSED(mev);
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::onKeyPressedEvent(KeypressedEvent *kev)
{
    if (m_useBeamActions)
        return;

    /// Control keys for interventonal Radiology simulations:
    switch(kev->getKey())
    {
        case 'D':
            applyAction(BeamAdapterAction::DROP_TOOL);
            break;
        case '2':
            applyAction(BeamAdapterAction::USE_TOOL_2);
            break;
        case '1':
            applyAction(BeamAdapterAction::USE_TOOL_1);
            break;
        case '0':
            applyAction(BeamAdapterAction::USE_TOOL_0);
            break;
        case 20: // droite = 20
            applyAction(BeamAdapterAction::SPIN_RIGHT);
            break;
        case 18: // gauche = 18
            applyAction(BeamAdapterAction::SPIN_LEFT);
            break;
        case 19: // fleche haut = 19
            applyAction(BeamAdapterAction::MOVE_FORWARD);
            break;
        case 21: // bas = 21
            applyAction(BeamAdapterAction::MOVE_BACKWARD);
            break;
        case '*':
            {
                if(m_RW)
                {
                    m_RW=false;

                }
                else
                {
                    m_FF = true;
                }
            }
            break;
        case '/':
            {
                if(m_FF)
                {
                    m_FF=false;
                }
                else
                {
                    m_RW = true;
                }
            }
            break;
    }
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::onBeginAnimationStep(const double dt)
{
    SOFA_UNUSED(dt);

    if (m_useBeamActions)
        return;

    BaseContext* context = getContext();
    auto xInstrTip = sofa::helper::getWriteOnlyAccessor(d_xTip);
    if(m_FF || m_RW)
    {
        int id = d_controlledInstrument.getValue();
        if (id >= (int)xInstrTip.size())
        {
            msg_warning()<<"Controlled Instument num "<<id<<" does not exist (size ="<< xInstrTip.size() <<") use instrument 0 instead";
            id=0;
        }
        if (m_FF)
        {
            if (!m_sensored)
                xInstrTip[id] += d_speed.getValue() * context->getDt();
        else
            {
                unsigned int newSensorData = m_currentSensorData + 1;

                while( m_sensorMotionData[newSensorData][0] < context->getTime() )
                {
                    m_currentSensorData = newSensorData;
                    newSensorData++;
                }
                if(newSensorData >= m_sensorMotionData.size())
                {
                    xInstrTip[id] = 0;
                }
                else
                {
                    xInstrTip[id] += m_sensorMotionData[m_currentSensorData][1];
                }
            }
        }
        if (m_RW)
        {
            xInstrTip[id] -= d_speed.getValue()* context->getDt();
            // verif min x :
            if ( xInstrTip[id] < 0.0)
            {
                xInstrTip[id] = 0.0;
                m_RW = false;
            }
        }
    }

    /// The tip of the instrument can not be further than its total length
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
        if (xInstrTip[i] > m_instrumentsList[i]->getRestTotalLength() )
            xInstrTip[i] = m_instrumentsList[i]->getRestTotalLength();

    applyInterventionalRadiologyController();
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::applyAction(sofa::beamadapter::BeamAdapterAction action)
{
    int id = d_controlledInstrument.getValue();
    if (id >= int(m_instrumentsList.size()))
    {
        msg_warning() << "Controlled Instrument num " << id << " does not exist (size =" << m_instrumentsList.size() << ").";
        return;
    }

    switch (action)
    {
    case BeamAdapterAction::NO_ACTION:
        break;
    case BeamAdapterAction::MOVE_FORWARD:
    {
        auto xInstrTip = sofa::helper::getWriteOnlyAccessor(d_xTip);
        xInstrTip[id] += d_step.getValue();
        break;
    }
    case BeamAdapterAction::MOVE_BACKWARD:
    {
        auto xInstrTip = sofa::helper::getWriteOnlyAccessor(d_xTip);
        xInstrTip[id] -= d_step.getValue();
        break;
    }
    case BeamAdapterAction::SPIN_RIGHT:
    {
        auto rotInstrument = sofa::helper::getWriteOnlyAccessor(d_rotationInstrument);
        rotInstrument[id] += d_angularStep.getValue();
        break;
    }
    case BeamAdapterAction::SPIN_LEFT:
    {
        auto rotInstrument = sofa::helper::getWriteOnlyAccessor(d_rotationInstrument);
        rotInstrument[id] -= d_angularStep.getValue();
        break;
    }
    case BeamAdapterAction::SWITCH_NEXT_TOOL:
    {
        if (id + 1 >= int(m_instrumentsList.size()))
            msg_warning() << "Switching to next tool is not possible, no more instrument in list.";
        else
            d_controlledInstrument.setValue(id + 1);
        break;
    }
    case BeamAdapterAction::SWITCH_PREVIOUS_TOOL:
    {
        if (id == 0)
            msg_warning() << "Switching to previous tool is not possible, already controlling first instrument.";
        else
            d_controlledInstrument.setValue(id - 1);
        break;
    }
    case BeamAdapterAction::USE_TOOL_0:
    {
        d_controlledInstrument.setValue(0);
        break;
    }
    case BeamAdapterAction::USE_TOOL_1:
    {
        if (1 >= m_instrumentsList.size())
            msg_warning() << "Controlled Instument num 1 do not exist (size =" << m_instrumentsList.size() << ") do not change the instrument id";
        else
            d_controlledInstrument.setValue(1);
        break;
    }
    case BeamAdapterAction::USE_TOOL_2:
    {
        if (2 >= m_instrumentsList.size())
            msg_warning() << "Controlled Instument num 2 do not exist (size =" << m_instrumentsList.size() << ") do not change the instrument id";
        else
            d_controlledInstrument.setValue(2);
        break;
    }
    case BeamAdapterAction::DROP_TOOL:
    {
        msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
    }
    }
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::computeInstrumentsCurvAbs(type::vector<Real> &newCurvAbs,
    const type::vector<Real>& tools_xBegin,
    const Real& totalLength)
{
    // Step 1 => put the noticeable Nodes
    // Step 2 => add the beams given the sampling parameters
    double maxAbsLength=0.0;
    Real xSampling = 0.0;
    for (sofa::Size i=0; i<m_instrumentsList.size(); i++)
    {
        type::vector<Real> xP_noticeable_I;
        type::vector< int > density_I;
        m_instrumentsList[i]->getSamplingParameters(xP_noticeable_I, density_I); // sampling of the different section of this instrument

        // check each interval of noticeable point to see if they go out (>0) and use corresponding density to sample the interval.
        for (int j=0; j<(int)(xP_noticeable_I.size()-1); j++)
        {
            const Real xP = xP_noticeable_I[j];
            const Real nxP = xP_noticeable_I[j + 1];

            //compute the corresponding abs curv of this "noticeable point" on the combined intrument deployed. 
            const Real curvAbs_xP = tools_xBegin[i] + xP; // xBegin = xend - instrument Total Length
            const Real curvAbs_nxP = tools_xBegin[i] + nxP;

            // In any case, the key points are added as soon as they are deployed
            if (curvAbs_xP > 0) {
                newCurvAbs.push_back(curvAbs_xP);
            }
            
            // compute interval between next point and previous one (0 for the first iter)
            const Real curvAbs_interval = (curvAbs_nxP - xSampling);

            if (curvAbs_interval > 0)
            {
                // compute the number of point of the emerged interval (if all the interval is emerged, i.e >0 , numNewNodes == density[j])
                Real ratio = Real(density_I[j]) / (nxP - xP);
                int numNewNodes = int(floor(curvAbs_interval * ratio)); // if density == 0, no sampling (numNewNodes == 0) 

                // Add the new points using reverse order iterator as they are computed deduce to next noticeable point
                for (int k = numNewNodes; k>0; k--)
                {
                    auto value = curvAbs_nxP - (k / ratio);
                    newCurvAbs.push_back(value);
                }

                xSampling = curvAbs_nxP;
            }
        }

        // After the end of the for loop above, we just have to process the
        // instrument last key point
        const Real lastxP = xP_noticeable_I[xP_noticeable_I.size()-1];
        const Real curvAbs_lastxP = tools_xBegin[i] + lastxP;
        if (curvAbs_lastxP > 0) {
            newCurvAbs.push_back(curvAbs_lastxP);
        }
    }


    // Step 1(bis) = add Nodes the curv_abs of the rigid parts border
    // When there are rigid segments, # of dofs is different than # of edges and beams
    const type::vector< Real > *rigidCurvAbs = &d_rigidCurvAbs.getValue();
    const auto nb = rigidCurvAbs->size();

    bool begin=true;
    if(nb>0 && (nb%2)==0)	// Make sure we have pairs of curv abs
    {
        RealConstIterator it;
        for(it=rigidCurvAbs->begin(); it!=rigidCurvAbs->end();)
        {
            Real abs;
            abs = *it++;
            if (abs < totalLength) // verify that the rigidified part is not outside the model
            {
                if(begin)
                {
                    newCurvAbs.push_back(abs-d_step.getValue());
                }
                else
                {
                    if (abs+d_step.getValue()<maxAbsLength)
                        newCurvAbs.push_back(abs+d_step.getValue());
                }
                begin = !begin;
                newCurvAbs.push_back(abs);
            }
        }
    }

    sortCurvAbs(newCurvAbs);
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::interventionalRadiologyCollisionControls(type::vector<Real> &xPointList,
                                                                                            type::vector<int> &idInstrumentList,
                                                                                            type::vector<int> &removeEdge)
{
    if(idInstrumentList.size() != xPointList.size())
    {
        msg_error()<<"The list do not have the same size";
        return;
    }

    // we enter the point from the tip to the end of the combined instrument
    // x_abs_curv provides the value of the curv abs along the combined instrument
    unsigned int node= m_nodeCurvAbs.size()-1;
    Real xAbsCurv = m_nodeCurvAbs[node];
    int firstInstruOnx = m_idInstrumentCurvAbsTable[node][0];

    type::vector<unsigned int> segRemove;

    for (unsigned int it=0; it<m_instrumentsList.size(); it++)
        segRemove.push_back(0);

    for (int i = static_cast<int>(xPointList.size()) - 1; i>=0; i--)
    {
        //1.  we determin if the poin ument
        int instrumentId = idInstrumentList[i];

        // x_max for the instrument that is controlled (not dropped part)
        Real xMaxControlled = m_instrumentsList[instrumentId]->getRestTotalLength();

        if (xPointList[i]>xMaxControlled)
        {
            unsigned int idInstr = idInstrumentList[i];
            segRemove[idInstr] = i;
            continue;
        }

        // 2. we assign the value of the curv abs for the point and the corresponding instrument
        Real xTipFirstInstruOnx = d_xTip.getValue()[firstInstruOnx];
        Real xBegin = xTipFirstInstruOnx - m_instrumentsList[firstInstruOnx]->getRestTotalLength();
        xPointList[i] = xAbsCurv - xBegin; // provides the "local" curv absc of the point (on the instrument reference)
        idInstrumentList[i] = firstInstruOnx;

        // 3. we look for the collision sampling of the current instrument in order to "place" the following point
        Real xIncr;
        m_instrumentsList[firstInstruOnx]->getCollisionSampling(xIncr, xPointList[i]);
        xAbsCurv -= xIncr;

        // the following point could not have x_abs_curv<0;
        if (xAbsCurv<0.0)
        {
            xAbsCurv=0.0;
            continue;
        }

        // the following point can be place on an other instrument
        while (node > 0 && xAbsCurv < m_nodeCurvAbs[node-1])
        {
            node--; // we change the beam support...
            if( m_idInstrumentCurvAbsTable[node][0] != firstInstruOnx)
            {
                // instrument has changed !!
                firstInstruOnx = m_idInstrumentCurvAbsTable[node][0];
                xAbsCurv = m_nodeCurvAbs[node];
                break;
            }

        }

    }

    for (unsigned int it=0; it<m_instrumentsList.size(); it++)
    {
        if(segRemove[it]!=0)
            removeEdge.push_back(segRemove[it]);
    }

    // A  way to detect if a collision point is "activated" or not=> look at its curv_abs  and if > 0, it is active
    // first, we need to compute abs_curv_point
    type::vector<Real> absCurvPoint;
    absCurvPoint.clear();

    for (unsigned int i=0; i<xPointList.size(); i++)
    {
        int instrumentId = idInstrumentList[i];

        // x_max for the instrument that is controlled (not dropped part)
        Real xMaxInstrument = m_instrumentsList[instrumentId]->getRestTotalLength();
        Real xTipInstrument = d_xTip.getValue()[instrumentId];
        Real xPoint= xPointList[i] - xMaxInstrument + xTipInstrument;

        absCurvPoint.push_back( xPoint );
    }

    // x point < epsilon... it is not activated`
    m_activatedPointsBuf.clear();
    m_activatedPointsBuf.push_back(false);
    for (unsigned int i=1; i<absCurvPoint.size(); i++)
    {
        Real xMaxInstrument = m_instrumentsList[idInstrumentList[i]]->getRestTotalLength();

        if (absCurvPoint[i] < std::numeric_limits<float>::epsilon() * xMaxInstrument
            || fabs(absCurvPoint[i] - absCurvPoint[i - 1]) < std::numeric_limits<float>::epsilon() * xMaxInstrument)
            m_activatedPointsBuf.push_back(false);
        else
            m_activatedPointsBuf.push_back(true);
    }
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::activateBeamListForCollision( type::vector<Real> &curv_abs,
                                                                                 type::vector< type::vector<int> > &idInstrumentTable)
{
    // 0. useful for rigidification
    const type::vector<Real>  *rigidCurvAbs = &d_rigidCurvAbs.getValue();

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
        else
        {
            instr0=idInstrumentTable[p+1][0];// get the first instrument

            ////// beam p should be assigned to instrument num instr0
            //3 .  Before assignement, verification that the beam is not on a rigidified part !
            bool rigid=false;
            RealConstIterator it = rigidCurvAbs->begin();
            while (it!=rigidCurvAbs->end())
            {
                if(curv_abs[p+1] <= (*it) )
                {
                    break;
                }
                else
                {
                    rigid = !rigid;
                    it++;
                }
            }
            if(!rigid)
            {
                m_instrumentsList[ instr0 ]->addCollisionOnBeam(p);
            }
        }
    }
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::applyInterventionalRadiologyController()
{
    const Real& threshold = d_threshold.getValue();

    
    // Create vectors with the CurvAbs of the noticiable points and the id of the corresponding instrument
    type::vector<Real> newCurvAbs;
    type::vector<type::vector<int>> idInstrumentTable;

    // ## STEP 1: Find the total length of the COMBINED INSTRUMENTS and the one for which xtip > 0 (so the one which are simulated)
    helper::AdvancedTimer::stepBegin("step1");
    Real totalLengthCombined=0.0;
    type::vector<Real> tools_xBegin, tools_xEnd;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        const Real& xend= d_xTip.getValue()[i];
        tools_xEnd.push_back(xend);
        tools_xBegin.push_back(xend - m_instrumentsList[i]->getRestTotalLength());
        
        if (xend> totalLengthCombined)
        {
            totalLengthCombined = xend;
        }

        // clear the present interpolation of the beams
        m_instrumentsList[i]->clear();
    }

    // create the first node (on x=0)
    if (totalLengthCombined > 0.0)
    {        
        newCurvAbs.push_back(0.0);
    }

    /// Some verif of step 1
    // if the totalLength is 0, move the first instrument
    if (totalLengthCombined < std::numeric_limits<float>::epsilon())
    {
        auto x = sofa::helper::getWriteOnlyAccessor(d_xTip);
        x[0] = std::numeric_limits<float>::epsilon();
        applyInterventionalRadiologyController();
        return;
    }
    helper::AdvancedTimer::stepEnd("step1");


    // ## STEP 2: Get the noticeable points that need to be simulated
    //     => Fill newCurvAbs which provides a vector with curvilinear abscissa of each simulated node    
    //     => xbegin (theoritical curv abs of the beginning point of the instrument (could be negative) xbegin= xtip - intrumentLength)
    helper::AdvancedTimer::stepBegin("step2");    
    computeInstrumentsCurvAbs(newCurvAbs, tools_xBegin, totalLengthCombined);

    //     => id_instrument_table which provides for each simulated node, the id of all instruments which belong this node
    fillInstrumentCurvAbsTable(newCurvAbs, tools_xBegin, tools_xEnd, idInstrumentTable);
    helper::AdvancedTimer::stepEnd("step2");

    // ## STEP 3: Re-interpolate the positions and the velocities
    helper::AdvancedTimer::stepBegin("step3");
    //    => Change curv if totalLength has changed: modifiedCurvAbs = newCurvAbs - current motion (Length between new and old tip curvAbs)
    type::vector<Real> modifiedCurvAbs; // This buffer will contain all deployed curvAbs minus current motion to mimic previous curvAbs (with 2 points with nearly the same abs at start) 
    totalLengthIsChanging(newCurvAbs, modifiedCurvAbs, idInstrumentTable); 

    //    => Get write access to current nodes/dofs
    Data<VecCoord>* datax = this->getMechanicalState()->write(core::VecCoordId::position());
    auto x = sofa::helper::getWriteOnlyAccessor(*datax);
    VecCoord xbuf = x.ref();

    sofa::Size nbrCurvAbs = newCurvAbs.size(); // number of simulated nodes
    if (nbrCurvAbs > x.size())
    {
        msg_warning() << "Parameters missmatch. There are more curv abscisses '" << nbrCurvAbs << "' than the number of dof: " << x.size();
        nbrCurvAbs = x.size();
    }

    const sofa::Size prev_nbrCurvAbs = m_nodeCurvAbs.size(); // previous number of simulated nodes;

    const sofa::Size nbrUnactiveNode = (m_numControlledNodes > nbrCurvAbs) ? m_numControlledNodes - nbrCurvAbs : 0; // m_numControlledNodes == nbr Dof | nbr of CurvAbs > 0
    const sofa::Size prev_nbrUnactiveNode = (m_numControlledNodes > prev_nbrCurvAbs) ? m_numControlledNodes - prev_nbrCurvAbs : 0;

    for (sofa::Index xId = 0; xId < nbrCurvAbs; xId++)
    {
        const sofa::Index globalNodeId = nbrUnactiveNode + xId; // position of the curvAbs in the dof buffer filled by the end
        const Real xCurvAbs = modifiedCurvAbs[xId];

        if ((xCurvAbs - std::numeric_limits<float>::epsilon()) > m_nodeCurvAbs.back() + threshold)
        {
            msg_warning() << "Case 1 should never happen while using totalLengthIsChanging. xCurvAbs = " << xCurvAbs 
                << " > m_nodeCurvAbs.back() = " << m_nodeCurvAbs.back() << " + threshold: " << threshold << "\n"
                << "\n | newCurvAbs: " << newCurvAbs                
                << "\n | modifiedCurvAbs: " << modifiedCurvAbs
                << "\n | previous nodeCurvAbs: " << m_nodeCurvAbs;
        }

        // The node position is not further than previous state, it can be interpolated straightfully using previous step positions
        sofa::Index prev_xId = 0;
        for (prev_xId = 0; prev_xId < m_nodeCurvAbs.size(); prev_xId++)
        {
            // if old_curvAbs[id] + threshold > current xabs, use this id to interpolate new curvAbs
            if ((m_nodeCurvAbs[prev_xId] + threshold) > xCurvAbs)
                break;
        }

        sofa::Index prev_globalNodeId = prev_nbrUnactiveNode + prev_xId;
        const Real prev_xCurvAbs = m_nodeCurvAbs[prev_xId];

        if (fabs(prev_xCurvAbs - xCurvAbs) < threshold)
        {
            x[globalNodeId] = xbuf[prev_globalNodeId]; // xBuf all initialised at start with d_startingPos
        }
        else
        {
            // the node must be interpolated using beam interpolation
            //find the instrument
            int id = m_idInstrumentCurvAbsTable[prev_xId][0];
            //find the good beam (TODO: do not work if xbegin of one instrument >0)
            int b = prev_xId - 1;
               
            // test to avoid wrong indices
            if (b < 0)
                x[globalNodeId] = d_startingPos.getValue();
            else
            {
                Transform global_H_interpol;
                const Real L = prev_xCurvAbs - m_nodeCurvAbs[b];
                Real baryCoef = 1.0;
                if (L < std::numeric_limits<float>::epsilon()) {
                    msg_error() << "Two consecutives curvAbs with the same position. Length is null. Using barycenter coefficient: baryCoef = 1";
                }
                else {
                    baryCoef = (xCurvAbs - m_nodeCurvAbs[b]) / L;
                }

                Transform Global_H_local0(xbuf[prev_globalNodeId - 1].getCenter(), xbuf[prev_globalNodeId - 1].getOrientation());
                Transform Global_H_local1(xbuf[prev_globalNodeId].getCenter(), xbuf[prev_globalNodeId].getOrientation());

                m_instrumentsList[id]->InterpolateTransformUsingSpline(global_H_interpol, baryCoef, Global_H_local0, Global_H_local1, L);

                x[globalNodeId].getCenter() = global_H_interpol.getOrigin();
                x[globalNodeId].getOrientation() = global_H_interpol.getOrientation();
            }
        }
    }
    helper::AdvancedTimer::stepEnd("step3");


    // ## STEP 4: Assign the beams
    helper::AdvancedTimer::stepBegin("step4");
    sofa::Size nbrBeam = newCurvAbs.size() - 1; // number of simulated beams
    const sofa::Size numEdges = m_numControlledNodes - 1;
    
    if (numEdges < nbrBeam) // verify that there is a sufficient number of Edge in the topology : TODO if not, modify topo !
    {
        msg_error() << "Not enough edges in the topology. Only: " << numEdges << " while nbrBeam = " << nbrBeam << ". Will simulate only " << numEdges << " beams.";
        nbrBeam = numEdges;
    }


    const type::vector<Real>& rotInstruments = d_rotationInstrument.getValue();
    for (unsigned int b=0; b< nbrBeam; b++)
    {
        const Real& x0 = newCurvAbs[b];
        const Real& x1 = newCurvAbs[b+1];

        for (unsigned int i=0; i<m_instrumentsList.size(); i++)
        {
            const Real& xmax = tools_xEnd[i];
            const Real& xmin = tools_xBegin[i];

            if (x0>(xmin- threshold) && x0<(xmax+ threshold) && x1>(xmin- threshold) && x1<(xmax+ threshold))
            {
                BaseMeshTopology::EdgeID eID = (BaseMeshTopology::EdgeID)(numEdges - nbrBeam + b);

                Real length = x1 - x0;
                Real x0_local = x0-xmin;
                Real x1_local = x1-xmin;

                Real theta = rotInstruments[i];

                m_instrumentsList[i]->addBeam(eID, length, x0_local, x1_local,theta );
            }
        }
    }
    helper::AdvancedTimer::stepEnd("step4");


    // ## STEP 5: Fix the not simulated nodes
    helper::AdvancedTimer::stepBegin("step5");
    unsigned int firstSimulatedNode = m_numControlledNodes - nbrBeam;

    //    => 1. Fix the nodes (beginning of the instruments) that are not "out"
    fixFirstNodesWithUntil(firstSimulatedNode);

    //    => 2. Fix the node that are "fixed". When there are rigid segments, # of dofs is different than # of edges and beams
    const std::vector< Real > *rigidCurvAbs = &d_rigidCurvAbs.getValue();

    bool rigid=false;
    unsigned int nbRigidAbs = rigidCurvAbs->size();
    if (nbRigidAbs>0 && (nbRigidAbs%2)==0)
    {
        RealConstIterator it = rigidCurvAbs->begin();

        for (unsigned int i=0; i<newCurvAbs.size(); i++)
        {
            if (newCurvAbs[i] < ((*it)+ std::numeric_limits<float>::epsilon()) && newCurvAbs[i] > ((*it)- std::numeric_limits<float>::epsilon())) // node= border of the rigid segment
            {
                if (!rigid)
                {
                    rigid=true;
                    m_fixedConstraint->addConstraint(firstSimulatedNode+i);
                }
                else
                {
                    m_fixedConstraint->addConstraint(firstSimulatedNode+i);
                    rigid=false;
                }
                it++;

                if(it==rigidCurvAbs->end())
                    break;
            }
            else
            {
                if(rigid)
                    m_fixedConstraint->addConstraint(firstSimulatedNode+i);
            }
        }
    }
    helper::AdvancedTimer::stepEnd("step5");

    // ## STEP 6: Activate Beam List for collision of each instrument
    activateBeamListForCollision(newCurvAbs,idInstrumentTable);

    // ## STEP 7: Save new computed Data
    m_nodeCurvAbs = newCurvAbs;
    m_idInstrumentCurvAbsTable = idInstrumentTable;
}



template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::totalLengthIsChanging(const type::vector<Real>& newNodeCurvAbs,
                                                                         type::vector<Real>& modifiedNodeCurvAbs,
                                                                         const type::vector< type::vector<int> >& newTable)
{

    // If total length is changing, it means that we need to simulate the fact that some beam will "move" during the time step.
    // thus during the interpolation process, where the position of the nodes is based on the previous position of the wire,
    // we initialize some points at a x_curv ref pos without the motion (computed by DLength)
    // due to the elasticity of the beam, the point will then naturally go the position that reespects the newNodeCurvAbs

    Real dLength = newNodeCurvAbs.back() - m_nodeCurvAbs.back();
    modifiedNodeCurvAbs = newNodeCurvAbs;

    // we look for the last value in the CurvAbs
    if (fabs(dLength) <= d_threshold.getValue())
        return;

    for (unsigned int i = 0; i < newTable.size(); i++)
    {
        if (newTable[i].size() == 1)
        {
            modifiedNodeCurvAbs[i] -= dLength;
            if (i > 1 && modifiedNodeCurvAbs[i] < modifiedNodeCurvAbs[i - 1])
            {
                modifiedNodeCurvAbs[i] = modifiedNodeCurvAbs[i - 1];
            }
        }
    }
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::sortCurvAbs(type::vector<Real> &curvAbs)
{
    // here we sort CurvAbs   
    std::sort(curvAbs.begin(), curvAbs.end());

    // a threshold is used to remove the values that are "too" close...
    const auto threshold = d_threshold.getValue();
    auto it = std::unique(curvAbs.begin(), curvAbs.end(), [threshold](const Real v1, const Real v2) {
        return fabs(v1 - v2) < threshold;
    });
    curvAbs.erase(it, curvAbs.end());
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::fillInstrumentCurvAbsTable(const type::vector<Real>& curvAbs, 
    const type::vector<Real>& tools_xBegin,
    const type::vector<Real>& tools_xEnd,
    type::vector< type::vector<int> >& idInstrumentTable)
{
    // here we build a table that provides the list of each instrument for each dof in the list of curvAbs
    // dofs can be shared by several instruments
    idInstrumentTable.clear();
    idInstrumentTable.resize(curvAbs.size());

    const auto threshold = d_threshold.getValue();
    for (unsigned int id = 0; id < m_instrumentsList.size(); id++)
    {
        // Get instrument absciss range
        Real xEnd = tools_xEnd[id];
        Real xBegin = tools_xBegin[id];

        // enlarge range to ensure to considere borders in absisses comparisons
        xBegin -= threshold;
        xEnd += threshold;

        // check curvAbs sorted value, if value is inside [xBegin, xBegin] of the tool add it to instrumentList. 
        for (unsigned int i = 0; i < curvAbs.size(); i++)
        {
            if (curvAbs[i] < xBegin) // still not inside range
                continue;

            if (curvAbs[i] > xEnd) // exit range
                break;

            idInstrumentTable[i].push_back(id);
        }
    }
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::fixFirstNodesWithUntil(unsigned int firstSimulatedNode)
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
bool InterventionalRadiologyController<DataTypes>::modifyTopology(void)
{
    return false;
}

template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::getInstrumentList(type::vector<fem::WireBeamInterpolation<DataTypes>*>& list)
{
    list = m_instrumentsList;
}

template <class DataTypes>
const type::vector< type::vector<int> >& InterventionalRadiologyController<DataTypes>::get_id_instrument_curvAbs_table() const
{
    return m_idInstrumentCurvAbsTable;
}

template <class DataTypes>
int InterventionalRadiologyController<DataTypes>::getTotalNbEdges() const
{
    return getContext()->getMeshTopology()->getNbEdges();
}


} // namespace sofa::component::controller::_interventionalradiologycontroller_


