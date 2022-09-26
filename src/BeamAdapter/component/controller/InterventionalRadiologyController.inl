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

namespace sofa::component::controller
{

namespace _interventionalradiologycontroller_
{

using type::vector;
using core::objectmodel::BaseContext;
using helper::WriteAccessor;
using core::objectmodel::KeypressedEvent;
using core::objectmodel::MouseEvent;


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
, d_curvAbs(initData(&d_curvAbs,"CurvAbs", "curvi-linear abscissa of the DOFs" ))
{
    m_fixedConstraint = nullptr;
    m_dropCall = false;
    m_sensored =false;
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::init()
{
    BaseContext* context = getContext();

    //get the pointers of the WireBeamInterpolations
    const type::vector<std::string>& instrumentPathList = d_instrumentsPath.getValue();
    if (instrumentPathList.empty())
    {
        WBeamInterpolation * wbinterpol= context->get< WBeamInterpolation >(BaseContext::Local);
        m_instrumentsList.push_back(wbinterpol);
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

    if(m_instrumentsList.empty())
        msg_error()<<"No Beam Interpolation found !!! the component can not work.";

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

    if (m_instrumentsList.size() == 0)
    {
        msg_error()<<"No instrument found ( no WireBeamInterpolation).";
        return;
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

    /// List of the instrument for which a "DROPPED" was proceeed TODO
    m_droppedInstruments.clear();


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

    reinit();
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
        std::cout << line << std::endl;
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

    WriteAccessor<Data<VecCoord> > x = *getMechanicalState()->write(core::VecCoordId::position());
    for(unsigned int i=0; i<x.size(); i++)
        x[i] = d_startingPos.getValue();
    m_numControlledNodes = x.size();

    applyInterventionalRadiologyController();
    reinit();
}

template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::reinit()
{
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::computeVertexT()
{
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
    /// Control keys for interventonal Radiology simulations:
    switch(kev->getKey())
    {
        case 'D':
            m_dropCall = true;
            break;

        case '2':
            {
                if (2 >= (int)m_instrumentsList.size() && f_printLog.getValue() )
                    msg_warning()<<"Controlled Instument num 2 do not exist (size ="<< m_instrumentsList.size() <<") do not change the instrument id";
                else
                    d_controlledInstrument.setValue(2);
            }
            break;

        case '1':
            {
                if (1 >= (int)m_instrumentsList.size() && f_printLog.getValue() )
                    msg_warning()<<"Controlled Instument num 1 do not exist (size ="<< m_instrumentsList.size() <<") do not change the instrument id";
                else
                    d_controlledInstrument.setValue(1);
            }
            break;

        case '0':
            d_controlledInstrument.setValue(0);
            break;

        case 20: // droite = 20
            {    
                auto rotInstrument = sofa::helper::getWriteOnlyAccessor(d_rotationInstrument);
                int id = d_controlledInstrument.getValue();
                rotInstrument[id] += d_angularStep.getValue();
            }
            break;
        case 18: // gauche = 18
            {
                int id = d_controlledInstrument.getValue();
                auto rotInstrument = sofa::helper::getWriteOnlyAccessor(d_rotationInstrument);
                rotInstrument[id] -= d_angularStep.getValue();
            }
            break;

        case 19: // fleche haut = 19
            {
                int id = d_controlledInstrument.getValue();
                auto xInstrTip = sofa::helper::getWriteOnlyAccessor(d_xTip);
                if (id >= (int)xInstrTip.size())
                {
                    msg_warning()<<"Controlled Instument num "<<id<<" does not exist (size ="<< xInstrTip.size() <<") use instrument 0 instead";
                    id=0;
                }
                xInstrTip[id] += d_step.getValue();
            }
            break;

        case 21: // bas = 21
            {
                int id = d_controlledInstrument.getValue();
                auto xInstrTip = sofa::helper::getWriteOnlyAccessor(d_xTip);
                if (id >= (int)xInstrTip.size())
                {
                    msg_warning()<<"Controlled Instument num "<<id<<" does not exist (size ="<< xInstrTip.size() <<") use instrument 0 instead.";
                    id=0;
                }
                xInstrTip[id] -= d_step.getValue();
            }
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
void InterventionalRadiologyController<DataTypes>::processDrop(unsigned int &previousNumControlledNodes,
                                                               unsigned int &segRemove)
{
    int ci = int(d_controlledInstrument.getValue());
    Real xMinOutLocal= 0.0;

    Real xBegin=0.0;

    // Quelque soit le resultat du process, le drop call est traite ici
    m_dropCall = false;

    // Step1 : quel est l'abscisse curviligne ou l'instrument controllé est seul ?
    for (unsigned int i=0; i<m_nodeCurvAbs.size(); i++)
    {
        // on parcourt toutes abscisse curv jusqu'à trouver un endroit où il n'y a qu'un seul instrument
        if (m_idInstrumentCurvAbsTable[i].size()==1)
        {
            // on vérifie qu'il s'agit de celui qui est controle
            if( ci ==m_idInstrumentCurvAbsTable[i][0])
            {
                 xBegin = d_xTip.getValue()[ci] - m_instrumentsList[ci]->getRestTotalLength();
                 xMinOutLocal = m_nodeCurvAbs[i] - xBegin;
                 break;
            }
            else
            {
                msg_error()<<" The control instrument is not out, drop is impossible.";
                return;
            }
        }
    }

    if(xMinOutLocal<=0.0)
    {
         msg_error()<<" x_min_out_local <= 0.0 The control instrument is not out, drop is impossible.";
         return;
    }

    // Step2 : on verifie que cette abscisse curviligne est compatible avec celle de l'instrument
    // (on ne peut pas casser un instrument s'il est à l'intérieur d'un autre instrument)
    int numBeamsNotUnderControlled = 0;
    Real xBreak;
    if( m_instrumentsList[ci]->breaksInTwo(xMinOutLocal, xBreak, numBeamsNotUnderControlled) )
    {
        msg_error()<<"Breaks in two process activated.";

        // for now, we simply suppress one more beam !
        m_numControlledNodes -= (numBeamsNotUnderControlled + 1);

        auto xEnds = sofa::helper::getWriteOnlyAccessor(d_xTip);
        xEnds[ci] =  xBegin +  xBreak;
    }
    else
        return;

    // Step3 : on remet à jour les abscisse curviligne des noeuds en virant toutes celles qui correspondent à la partie
    // cassée
    Real eps=d_threshold.getValue();
    for (unsigned int i=0; i<m_nodeCurvAbs.size(); i++)
    {
        if( m_nodeCurvAbs[i] > (xBegin +  xBreak + eps) )
        {
            type::removeIndex(m_nodeCurvAbs,i);
            type::removeIndex(m_idInstrumentCurvAbsTable, i);
            i--;
        }
    }
    segRemove = 1;
    previousNumControlledNodes =m_numControlledNodes;
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::interventionalRadiologyComputeSampling(type::vector<Real> &newCurvAbs,
                                                                                          type::vector< type::vector<int> > &idInstrumentTable,
                                                                                          const type::vector<Real> &xBegin,
                                                                                          const Real& xend)

{
    // Step 1 = put the noticeable Nodes
    double maxAbsLength=0.0;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        type::vector<Real> xP_noticeable_I;
        type::vector< int > density_I;
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
            if (abs<xend) // verify that the rigidified part is not outside the model
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

    // Step 2 => add the beams given the sampling parameters
    Real xSampling = 0.0;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        type::vector<Real> xPNoticeableI;
        type::vector< int > density_I;
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

        if (absCurvPoint[i] < 0.00000001*xMaxInstrument || fabs(absCurvPoint[i] - absCurvPoint[i-1])<0.00000001*xMaxInstrument)
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
    /// Create vectors with the CurvAbs of the noticiable points and the id of the corresponding instrument
    type::vector<Real> newCurvAbs;

    /// In case of drop:
    unsigned int previousNumControlledNodes = m_numControlledNodes;
    unsigned int seg_remove = 0;

    if (m_dropCall)
        processDrop(previousNumControlledNodes, seg_remove);

    /// STEP 1
    /// Find the total length of the COMBINED INSTRUMENTS and the one for which xtip > 0 (so the one which are simulated)
    Real xend;
    helper::AdvancedTimer::stepBegin("step1");
    Real totalLengthCombined=0.0;
    type::vector<Real> xbegin;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
         xend= d_xTip.getValue()[i];
        Real xb = xend - m_instrumentsList[i]->getRestTotalLength();
        xbegin.push_back(xb);

        if (xend> totalLengthCombined)
        {
            totalLengthCombined=xend;
        }

        // clear the present interpolation of the beams
        m_instrumentsList[i]->clear();

        if( xend > 0.0)
        {
            // create the first node (on x=0)
            newCurvAbs.push_back(0.0);
        }
    }

    /// Some verif of step 1
    // if the totalLength is 0, move the first instrument
    if(totalLengthCombined<0.0001)
    {
        auto x = sofa::helper::getWriteOnlyAccessor(d_xTip);
        x[0]=0.0001;
        applyInterventionalRadiologyController();
        return;
    }
    helper::AdvancedTimer::stepEnd("step1");


    /// STEP 2:
    /// get the noticeable points that need to be simulated
    // Fill=> newCurvAbs which provides a vector with curvilinear abscissa of each simulated node
    //     => id_instrument_table which provides for each simulated node, the id of all instruments which belong this node
    //     => xbegin (theoritical curv abs of the beginning point of the instrument (could be negative) xbegin= xtip - intrumentLength)
    helper::AdvancedTimer::stepBegin("step2");
    type::vector<type::vector<int>> idInstrumentTable;
    interventionalRadiologyComputeSampling(newCurvAbs,idInstrumentTable, xbegin, totalLengthCombined);
    helper::AdvancedTimer::stepEnd("step2");


    /// STEP 3
    /// Re-interpolate the positions and the velocities
    helper::AdvancedTimer::stepBegin("step3");
    unsigned int nbeam=newCurvAbs.size()-1; // number of simulated beams
    unsigned int nnode=newCurvAbs.size(); // number of simulated nodes

    unsigned int nnode_old= m_nodeCurvAbs.size(); // previous number of simulated nodes;
    Data<VecCoord>* datax = this->getMechanicalState()->write(core::VecCoordId::position());                
    auto x = sofa::helper::getWriteOnlyAccessor(*datax);

    VecCoord xbuf =x;

    type::vector<Real> modifiedCurvAbs;

    totalLengthIsChanging(newCurvAbs, modifiedCurvAbs, idInstrumentTable);

    Real xmax_prev = m_nodeCurvAbs[m_nodeCurvAbs.size()-1];

    for (unsigned int p=0; p<nbeam+1; p++)
    {
        int idP = m_numControlledNodes-nnode + p;
        Real xabs = modifiedCurvAbs[p];

        // 2 cases:  TODO : remove first case
            //1. the abs curv is further than the previous state of the instrument
            //2. this is not the case and the node position can be interpolated using previous step positions
        if(xabs > xmax_prev + d_threshold.getValue())
        {
            msg_error_when(f_printLog.getValue())
                << "case 1 should never happen ==> avoid using totalLengthIsChanging ! xabs = " << xabs << " - xmax_prev = " << xmax_prev
                << "newCurvAbs  = " << newCurvAbs << "  previous nodeCurvAbs" << m_nodeCurvAbs
                << "modifiedCurvAbs =" << modifiedCurvAbs;
            // case 1 (the abs curv is further than the previous state of the instrument)
            // verifier qu'il s'agit bien d'un instrument qu'on est en train de controller
            // interpoler toutes les positions "sorties" de l'instrument en supprimant l'ajout de dx qu'on vient de faire
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

            int idP0 =  previousNumControlledNodes + seg_remove - nnode_old + p0 ;

            if(fabs(m_nodeCurvAbs[p0]-xabs)<d_threshold.getValue())
                x[idP] = xbuf[idP0];
            else
            {
                // the node must be interpolated using beam interpolation
                    //find the instrument
                int id = m_idInstrumentCurvAbsTable[p0][0];
                //find the good beam (TODO: do not work if xbegin of one instrument >0)
                int b = p0-1;
                // test to avoid wrong indices
                if (b<0)
                    x[p]=d_startingPos.getValue();
                else
                {
                    Transform global_H_interpol;
                    Real ratio = (xabs - m_nodeCurvAbs[b])/ (m_nodeCurvAbs[p0]-m_nodeCurvAbs[b]);
                    Transform Global_H_local0(xbuf[idP0-1].getCenter(),xbuf[idP0-1].getOrientation() ), Global_H_local1(xbuf[idP0].getCenter(),xbuf[idP0].getOrientation() );

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

    // verify that there is a sufficient number of Edge in the topology : TODO if not, modify topo !
    if(numEdges<nbeam)
    {
        if (f_printLog.getValue())
        {
            msg_error()<<"Not enough edges in the topology.";
        }
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

                m_instrumentsList[i]->addBeam(eID, length, x0_local, x1_local,theta );

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

    //2. Fix the node that are "fixed"
    // When there are rigid segments, # of dofs is different than # of edges and beams
    const std::vector< Real > *rigidCurvAbs = &d_rigidCurvAbs.getValue();

    bool rigid=false;
    unsigned int nbRigidAbs = rigidCurvAbs->size();
    if (nbRigidAbs>0 && (nbRigidAbs%2)==0)
    {
        RealConstIterator it = rigidCurvAbs->begin();

        for (unsigned int i=0; i<newCurvAbs.size(); i++)
        {
            if (newCurvAbs[i] < ((*it)+0.001) && newCurvAbs[i] > ((*it)-0.001)) // node= border of the rigid segment
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

    /// STEP 6
    /// Activate Beam List for collision of each instrument
    activateBeamListForCollision(newCurvAbs,idInstrumentTable);

    m_nodeCurvAbs = newCurvAbs;
    m_idInstrumentCurvAbsTable = idInstrumentTable;
    helper::AdvancedTimer::stepEnd("step5");
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
void InterventionalRadiologyController<DataTypes>::sortCurvAbs(type::vector<Real> &curvAbs,
                                                               type::vector<type::vector<int> >& idInstrumentTable)
{
    // here we sort CurvAbs   
    std::sort(curvAbs.begin(), curvAbs.end());

    // copy threshold in variable member for std::unique
    m_threshold = d_threshold.getValue();

    // a threshod is used to remove the values that are "too" similars...
    auto it = std::unique(curvAbs.begin(), curvAbs.end(), [this](const Real v1, const Real v2) {
        return fabs(v1 - v2) < this->m_threshold;
    });
    curvAbs.erase(it, curvAbs.end());

    // here we build a table that provides the list of each instrument for each dof in the list of curvAbs
    // dofs can be shared by several instruments
    idInstrumentTable.clear();
    idInstrumentTable.resize(curvAbs.size());

    for (unsigned int id = 0; id < m_instrumentsList.size(); id++)
    {
        // Get instrument absciss range
        Real xEnd = d_xTip.getValue()[id];
        Real xBegin = xEnd - m_instrumentsList[id]->getRestTotalLength();

        // enlarge range to ensure to considere borders in absisses comparisons
        xBegin -= m_threshold;
        xEnd += m_threshold;

        // check curvAbs sorted value, if value is inside [xBegin, xBegin] of the tool add it to instrumentList. 
        for (unsigned int i = 0; i < curvAbs.size(); i++)
        {
            auto xCurv = curvAbs[i];
            if (xCurv < xBegin) // still not inside range
                continue;

            if (xCurv > xEnd) // exit range
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


} // namespace _interventionalradiologycontroller_

} // namespace sofa::component::controller
