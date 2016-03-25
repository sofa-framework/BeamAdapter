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

#ifndef SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGYCONTROLLER_INL
#define SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGYCONTROLLER_INL

#include "InterventionalRadiologyController.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <SofaBaseTopology/EdgeSetGeometryAlgorithms.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/helper/AdvancedTimer.h>



namespace sofa
{

namespace component
{



namespace controller
{

template <class DataTypes>
InterventionalRadiologyController<DataTypes>::InterventionalRadiologyController()
: m_instrumentsPath(initData(&m_instrumentsPath,"instruments", "List of paths to WireInterpolation components on the scene"))
, controlledInstrument(initData(&controlledInstrument, 0, "controlledInstrument", "provide the id of the interventional radiology instrument which is under control: press contr + number to change it"))
, xtip(initData(&xtip,"xtip", "curvilinear abscissa of the tip of each interventional radiology instrument"))
, rotationInstrument(initData(&rotationInstrument,"rotationInstrument", "angle of rotation for each interventional radiology instrument"))
, step(initData(&step,(Real)0.1,"step","base step when changing beam length"))
, angularStep(initData(&angularStep,(Real)(3.1416/20.0),"angularStep","base step when changing beam angle"))
, speed(initData(&speed,(Real)0.0,"speed","continuous beam length increase/decrease"))
, startingPos(initData(&startingPos,Coord(),"startingPos","starting pos for inserting the instrument"))
, threshold(initData(&threshold, (Real)0.01, "threshold", "threshold for controller precision which is homogeneous to the unit of length used in the simulation"))
, m_rigidCurvAbs(initData(&m_rigidCurvAbs, "rigidCurvAbs", "pairs of curv abs for beams we want to rigidify"))
, motionFilename(initData(&motionFilename, "motionFilename", "text file that includes tracked motion from optical sensor"))
, indexFirstNode(initData(&indexFirstNode, (unsigned int) 0, "indexFirstNode", "first node (should be fixed with restshape)"))
{
    //edgeSetInNode=true;
    _fixedConstraint = NULL;
    this->dropCall = false;
    sensored =false;

    /*
    if (!_instrumentlist.empty())
    {
        //m_instrumentsList.resize( _instrumentlist.size() );
        for(unsigned int id=0;id<_instrumentlist.size();id++)
        {
            m_instrumentsList.push_back(_instrumentlist[id]);
                         sout<<"====================================== "<<id<<". "<<m_instrumentsList[id]->getName()<<"  radius ="<<m_instrumentsList[id]->radius.getValue()<<sendl;
        }
    }
    */


}


/*
    const helper::vector<std::string>& precondNames = f_preconditioners.getValue();
    if (precondNames.empty() || !use_precond.getValue()) {
        c->get<sofa::core::behavior::LinearSolver>(&solvers,BaseContext::SearchDown);
    } else {
        for (unsigned int i=0;i<precondNames.size();++i) {
            sofa::core::behavior::LinearSolver* s = NULL;
            c->get(s, precondNames[i]);
            if (s) solvers.push_back(s);
            else serr << "Solver \"" << precondNames[i] << "\" not found." << sendl;
        }
    }

    for (unsigned int i=0;i<solvers.size();++i) {
            if (solvers[i] && solvers[i] != this) {
                    this->preconditioners.push_back(solvers[i]);
            }
    }
*/



template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::init()
{
    sofa::core::objectmodel::BaseContext* context = this->getContext();

    //// get the pointers of the WireBeamInterpolations ///////

    const helper::vector<std::string>& instrumentPathList = m_instrumentsPath.getValue();
    if (instrumentPathList.empty()) {
        WBeamInterpolation * wbinterpol= context->get< WBeamInterpolation >(core::objectmodel::BaseContext::Local);
        m_instrumentsList.push_back(wbinterpol);
        // TODO: edit m_instrumentsPath with the name of wbinterpol
    } else {

        for (unsigned int i=0;i<instrumentPathList.size();++i)
        {
            WBeamInterpolation * wbinterpol = NULL;
            context->get(wbinterpol,instrumentPathList[i]);
            if( wbinterpol)
                m_instrumentsList.push_back(wbinterpol)  ;
            else
                serr << "Interpolation of instrument "<<instrumentPathList[i]<< "  not found"<<sendl;

        }
    }

    if(m_instrumentsList.empty())
        serr<<" no Beam Interpolation found !!! the component can not work"<<sendl;

    /////////////////


     activated_Points_buf.clear();

     if(speed.getValue()>0)
     {
         FF=true;
         RW=false;
         sensored = false;
     }
    if (!motionFilename.getValue().empty())
    {
        FF = true; sensored = true; currentSensorData = 0;
        loadMotionData(motionFilename.getValue());
    }


    // get the topology & tools for topology changes
    _topology = context->getMeshTopology();


    //context->get< sofa::component::fem::WireBeamInterpolation<DataTypes> > ( &m_instrumentsList, core::objectmodel::BaseContext::SearchDown );
    if (m_instrumentsList.size() == 0)
    {
        if (this->f_printLog.getValue())
        {
            serr<<" No instrument found ( no WireBeamInterpolation) !"<<sendl;
        }
        return;
    }
    sofa::helper::vector<Real> &x_instr_tip = (*this->xtip.beginEdit());
    x_instr_tip.resize(m_instrumentsList.size());
    if (this->f_printLog.getValue())
    {
        sout << "Size of xtip (curvilinear abscissa of the tip of each interventional radiology instrument) is now " << x_instr_tip.size() << sendl;
        for (unsigned int i =0; i< x_instr_tip.size(); i++)
            sout << "xtip[" << i << "] = " << x_instr_tip[i] << sendl;

    }
    this->xtip.endEdit();

    sofa::helper::vector<Real> &angle_Instrument = (*this->rotationInstrument.beginEdit());
    angle_Instrument.resize(m_instrumentsList.size());
    this->rotationInstrument.endEdit();

    if (this->f_printLog.getValue()) sout<<"  +++++++++++++++++ \n Init Interventional Radiology  : instruments found (should be sorted by decreasing radius)"<<sendl;
    for(unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        if (this->f_printLog.getValue())  sout<<" "<<i<<". "<<m_instrumentsList[i]->getName()<<"  radius ="<<m_instrumentsList[i]->radius.getValue()<<sendl;
        m_instrumentsList[i]->setControlled(true);
        /// \todo : automatically sort by radius
    }
    //sout<<"  +++++++++++++++++ "<<sendl;

    context->get(_fixedConstraint);
    if(_fixedConstraint==NULL)
    {
        if (this->f_printLog.getValue())
        {
            serr<<" No _fixedConstraint found !"<<sendl;
        }
        /// \todo: dynamically add a fixedConstraint in the node :
    }


    ////////// List of the instrument for which a "DROPPED" was proceeed TODO
    this->droppedInstruments.clear();



    // the controller must listen to the event (in particular BeginAnimationStep event)
    this->f_listening.setValue(true);


    nodeCurvAbs.clear();
    id_instrument_curvAbs_table.clear();
    nodeCurvAbs.push_back(0.0);
    sofa::helper::vector<int> list_init;

    for(unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        list_init.push_back(i);
    }
    id_instrument_curvAbs_table.push_back(list_init);

    Inherit::init();


    /// \todo : VERIFIER QUE LA TOPOLOGIE EST CELLE D'UN WIRE


    reinit();

}

template<class DataTypes>
void InterventionalRadiologyController<DataTypes>::loadMotionData(std::string filename)
{
    if (!sofa::helper::system::DataRepository.findFile(filename))
    {
        std::cerr << "File " << filename << " not found " << std::endl;
        return;
    }
    std::ifstream file(filename.c_str());

    std::string line;
    Vec3 result;
    while( std::getline(file,line) )
    {
        if (line.empty()) continue;
        std::cout << line << std::endl;
        std::istringstream values(line);
        values >> result[0] >> result[1] >> result[2];
        result[0] /= 1000;
        sensorMotionData.push_back(result);
    }

    file.close();

    std::cout << "Motion sensor: " << sensorMotionData.size() << " events recorded" << std::endl;
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::bwdInit()
{
    // assign the starting pos to each point of the Mechanical State //
    Coord stPos =startingPos.getValue();
    stPos.getOrientation().normalize();
    startingPos.setValue(stPos);

    //VecCoord& x = *(this->getMechanicalState()->read(sofa::core::ConstVecCoordId::position())->getValue());
    helper::WriteAccessor<Data<VecCoord> > x = *this->getMechanicalState()->write(sofa::core::VecCoordId::position());
    for(unsigned int i=0; i<x.size(); i++)
    {
        x[i] = startingPos.getValue();
    }
    numControlledNodes = x.size();

    applyInterventionalRadiologyController();
    reinit();


}

template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::reinit()
{
    applyController();
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::computeVertexT()
{
    // voir definition dans edgeSetController

}


/*!
 * \todo fix the mouse event with better controls
 */
template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::onMouseEvent(core::objectmodel::MouseEvent * /*mev*/)
{
//    int PosX = mev->getPosX();
//    int PosY = mev->getPosY();


//    //Translation input
//    Real PosYcorr = 0.0;
//    int idy = controlledInstrument.getValue();
//    sofa::helper::vector<Real> &x_instr_tip = (*this->xtip.beginEdit());
//    if (idy >= (int)x_instr_tip.size()){
//		if (this->f_printLog.getValue())
//			serr<<"WARNING controlled Instument num "<<idy<<" do not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead"<<sendl;
//        idy=0;
//    }
////    PosYcorr = -PosY*0.2;
////    x_instr_tip[idy] += PosYcorr;
//
//    if (this->f_printLog.getValue())
//    {
//    	sout << "[onMouseEvent] Size of xtip is now " << x_instr_tip.size() << sendl;
//    	for (unsigned int i =0; i< x_instr_tip.size(); i++)
//    		sout << "xtip[" << i << "] = " << x_instr_tip[i] << sendl;
//
//    }
//
//    this->xtip.endEdit();



    //Rotation input
//    Real PosXcorr = 0.0;
//    int idx = controlledInstrument.getValue();
//    sofa::helper::vector<Real> &rot_instrument = (*this->rotationInstrument.beginEdit());
//    PosXcorr = PosX*0.015;
//    rot_instrument[idx] += PosXcorr;

}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::onKeyPressedEvent(core::objectmodel::KeypressedEvent *kev)
{

    ///////////////////////////////// Control keys for interventonal Radiology simulations:
    switch(kev->getKey())
    {
        case 'D':
            this->dropCall = true;
            break;


        case '3':
            {
                if (3 >= (int)m_instrumentsList.size() && this->f_printLog.getValue() ){
                    serr<<"WARNING controlled Instument num 3 do not exist (size ="<< m_instrumentsList.size() <<") do not change the instrument id"<<sendl;
                }
                else
                    this->controlledInstrument.setValue(3);
            }
            break;


        case '2':
            {
                if (2 >= (int)m_instrumentsList.size() && this->f_printLog.getValue() ){
                    serr<<"WARNING controlled Instument num 2 do not exist (size ="<< m_instrumentsList.size() <<") do not change the instrument id"<<sendl;
                }
                else
                    this->controlledInstrument.setValue(2);
            }
            break;

        case '1':
            {
                if (1 >= (int)m_instrumentsList.size() && this->f_printLog.getValue() ){
                    serr<<"WARNING controlled Instument num 1 do not exist (size ="<< m_instrumentsList.size() <<") do not change the instrument id"<<sendl;
                }
                else
                    this->controlledInstrument.setValue(1);
            }
            break;

        case '0':
            this->controlledInstrument.setValue(0);
            break;

        case 'A':
            {

                int id = controlledInstrument.getValue();
                sofa::helper::vector<Real> &rot_instrument = (*this->rotationInstrument.beginEdit());
                rot_instrument[id] += angularStep.getValue();
                this->rotationInstrument.endEdit();

            }
            break;
        case 'E':
            {

                int id = controlledInstrument.getValue();
                sofa::helper::vector<Real> &rot_instrument = (*this->rotationInstrument.beginEdit());
                rot_instrument[id] -= angularStep.getValue();
                this->rotationInstrument.endEdit();

            }
            break;

        case '+':
            {
                int id = controlledInstrument.getValue();
                sofa::helper::vector<Real> &x_instr_tip = (*this->xtip.beginEdit());
                if (id >= (int)x_instr_tip.size()){
                    serr<<"WARNING controlled Instument num "<<id<<" does not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead"<<sendl;
                    id=0;
                }
                x_instr_tip[id] += step.getValue();
                /*for (unsigned int i =0; i<= id; i++)
                {
                    x_instr_tip[i] += step.getValue();
                }

                for (unsigned int i =0; i< 10; i++)
                {
                    serr<<"The id is "<<id<<" ##############################################################"<<sendl;
                }*/


                this->xtip.endEdit();
            }
            break;

        case '-':
            {
                int id = controlledInstrument.getValue();
                sofa::helper::vector<Real> &x_instr_tip = (*this->xtip.beginEdit());
                if (id >= (int)x_instr_tip.size()){
                    serr<<"WARNING controlled Instument num "<<id<<" does not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead"<<sendl;
                    id=0;
                }
                x_instr_tip[id] -= step.getValue();

                this->xtip.endEdit();
            }
            break;

        case '*':
            {
                if(RW)
                {
                    RW=false;

                }
                else
                {
                    FF = true;
                }
            }
            break;

        case '/':
            {
                if(FF)
                {
                    FF=false;
                }
                else
                {
                    RW = true;
                }
            }
            break;


    }
    if (this->f_printLog.getValue())
    {
        std::cout<<" key :"<<kev->getKey()<<" controlledInstrument is "<< this->controlledInstrument.getValue()<<std::endl;
    }

}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::onBeginAnimationStep(const double /*dt*/)
{
    sofa::core::objectmodel::BaseContext* context = this->getContext();
    sofa::helper::vector<Real> &x_instr_tip = (*this->xtip.beginEdit());
    if(FF || RW)
    {
        int id = controlledInstrument.getValue();
        if (id >= (int)x_instr_tip.size()){
            serr<<"WARNING controlled Instument num "<<id<<" does not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead"<<sendl;
            id=0;
        }
        if (FF)
        {
            if (!sensored)
                x_instr_tip[id] += speed.getValue() * context->getDt();
            /*{for (unsigned int i=0;i<=id;i++)
                {
                    x_instr_tip[i] += speed.getValue() * context->getDt();
                }
                for (unsigned int i=0;i<=id;i++)
                {
                    serr<<"Je passe ici ####################################################################"<<sendl;
                }}*/
        else
            {
                // x_instr_tip[id] += sensorMotionData[1000 * context->getTime()][1];
                unsigned int newSensorData = currentSensorData + 1;

                // std::cout << "Current time is " << context->getTime() << std::endl;
                while( sensorMotionData[newSensorData][0] < context->getTime() )
                {
                    currentSensorData = newSensorData;
                    newSensorData++;
                }
                if(newSensorData >= sensorMotionData.size())
                {
                    x_instr_tip[id] = 0;
                }
                else
                {
                    std::cout << "Current data is " << currentSensorData << " // " << sensorMotionData[currentSensorData][1] << std::endl;
                    x_instr_tip[id] += sensorMotionData[currentSensorData][1];
                    std::cout << "motion is --- " << x_instr_tip[id] << std::endl;
                }
            }
        }
        if (RW)
        {
            x_instr_tip[id] -= speed.getValue()* context->getDt();
            // verif min x :
            if ( x_instr_tip[id] < 0.0)
            {
                x_instr_tip[id] = 0.0;
                RW = false;
            }
        }
    }


    ///////// The tip of the instrument can not be further than its total length
    ////// \todo => dropp !!
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
        if (x_instr_tip[i] > this->m_instrumentsList[i]->getRestTotalLength() )
            x_instr_tip[i] = this->m_instrumentsList[i]->getRestTotalLength();
    if (this->f_printLog.getValue())
    {
        sout << "[onBeginAnimationStep] Size of xtip is now " << x_instr_tip.size() << sendl;
        for (unsigned int i =0; i< x_instr_tip.size(); i++)
            sout << "xtip[" << i << "] = " << x_instr_tip[i] << sendl;

    }

    this->xtip.endEdit();

    applyInterventionalRadiologyController();
    //VecElementID
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::processDrop(unsigned int &previousNumControlledNodes, unsigned int &seg_remove)
{
    int ci = (int) controlledInstrument.getValue();
    if ( this->f_printLog.getValue())
        std::cout<<" \n +++++ Drop called on instrument "<<ci<<" \n +++++"<<std::endl;
    Real x_min_out_local= 0.0;

    Real xbegin=0.0;

    // Quelque soit le resultat du process, le drop call est traite ici
    this->dropCall = false;

    // Step1 : quel est l'abscisse curviligne ou l'instrument controllé est seul ?
    for (unsigned int i=0; i<nodeCurvAbs.size(); i++)
    {
        // on parcourt toutes abscisse curv jusqu'à trouver un endroit où il n'y a qu'un seul instrument
        if (id_instrument_curvAbs_table[i].size()==1)
        {
            // on vérifie qu'il s'agit de celui qui est controle
            if( ci ==id_instrument_curvAbs_table[i][0])
            {
                 xbegin = this->xtip.getValue()[ci] - this->m_instrumentsList[ci]->getRestTotalLength();
                 x_min_out_local = nodeCurvAbs[i] - xbegin;
                 break;
            }
            else
            {
                sout<<" The control instrument is not out, drop is impossible"<<sendl;
                return;
            }
        }
    }

    if(x_min_out_local<=0.0)
    {
         sout<<" x_min_out_local <=0.0 The control instrument is not out, drop is impossible"<<sendl;
         return;
    }

    // Step2 : on verifie que cette abscisse curviligne est compatible avec celle de l'instrument
    // (on ne peut pas casser un instrument s'il est à l'intérieur d'un autre instrument)
    int numBeamsNotUnderControlled=0.0;
    Real x_break;
    if( m_instrumentsList[ci]->breaksInTwo(x_min_out_local, x_break, numBeamsNotUnderControlled) )
    {
        sout<<" BREAKS IN TWO PROCESS ACTIVATED !!!"<<sendl;

        // TODO => ADD POINT  ( coord : DUPLICATE POINT [numNodes - 1 - numBeamsNotUnderControlled] )
        //                    ( renumbering so that it follows the point used for duplication )
        //                    ( modif seg[numSeg - numBeamsNotUnderControlled] and to replace old node by the new one)


        // for now, we simply suppress one more beam !
        numControlledNodes -= (numBeamsNotUnderControlled + 1);

        sofa::helper::vector<Real> &xends = (*this->xtip.beginEdit());
        xends[ci] =  xbegin +  x_break;
        this->xtip.endEdit();
    }
    else
        return;





    // Step3 : on remet à jour les abscisse curviligne des noeuds en virant toutes celles qui correspondent à la partie
    // cassée
    Real eps=threshold.getValue();
    for (unsigned int i=0; i<nodeCurvAbs.size(); i++)
    {
        if( nodeCurvAbs[i] > (xbegin +  x_break + eps) )
        {
            sofa::helper::removeIndex(nodeCurvAbs,i);
            sofa::helper::removeIndex(id_instrument_curvAbs_table, i);
            i--;
        }
    }
    seg_remove = 1;
    previousNumControlledNodes =numControlledNodes;
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::interventionalRadiologyComputeSampling(sofa::helper::vector<Real> &newCurvAbs, sofa::helper::vector< sofa::helper::vector<int> > &id_instrument_table,
                                                                               const sofa::helper::vector<Real> &xbegin, const Real& xend)

{

    // Step 1 = put the noticeable Nodes
    double maxAbsLength=0.0;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        helper::vector<Real> xP_noticeable_I;
         helper::vector< int > density_I;
        this->m_instrumentsList[i]->getSamplingParameters(xP_noticeable_I, density_I);

        for (unsigned int j=0; j<xP_noticeable_I.size(); j++){

            //compute the corresponding abs curv of this "noticeable point" on the combined intrument
            Real curvAbs_xP = xbegin[i] + xP_noticeable_I[j];
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
    helper::ReadAccessor< Data< helper::set< Real > > > rigidCurvAbs = m_rigidCurvAbs;
    int nb = rigidCurvAbs->size();

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
                    newCurvAbs.push_back(abs-this->step.getValue());
                }
                else
                {
                    if (abs+this->step.getValue()<maxAbsLength)
                        newCurvAbs.push_back(abs+this->step.getValue());
                }
                begin = !begin;
                newCurvAbs.push_back(abs);


            }
        }
    }




    // Step 2 => add the beams given the sampling parameters

    Real x_sampling = 0.0;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
        helper::vector<Real> xP_noticeable_I;
         helper::vector< int > density_I;
        this->m_instrumentsList[i]->getSamplingParameters(xP_noticeable_I, density_I);


        for (unsigned int j=0; j<density_I.size(); j++){

            //compute the corresponding abs curv of this "noticeable point" on the combined intrument
            Real curvAbs_xP = xbegin[i] + xP_noticeable_I[j+1];

            // use density parameter (size = xP_noticeable_I -1 )
            if (curvAbs_xP > x_sampling && density_I[j]>0)
            {

               // std::cout<<" ratio computation : xP_noticeable_I[j+1] ="<<xP_noticeable_I[j+1]<<" - xP_noticeable_I[j] ="<<xP_noticeable_I[j]<<std::endl;
                Real ratio =(Real) density_I[j] / (xP_noticeable_I[j+1]  - xP_noticeable_I[j]) ;
                int  numNewNodes = (int) floor( (curvAbs_xP- x_sampling)  *ratio) ;

                for (int k=0; k<numNewNodes; k++)
                {
                    newCurvAbs.push_back( xP_noticeable_I[j+1] + xbegin[i]  - (k+1) * (1/ratio) );

                }


                x_sampling = curvAbs_xP;

            }
        }
    }


    sortCurvAbs(newCurvAbs, id_instrument_table);

    if (this->f_printLog.getValue())
            std::cout<<" noticeable points :["<<newCurvAbs<<"] -   id_instrument_table."
            <<id_instrument_table.size() <<"  : "<<id_instrument_table<<std::endl;
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::interventionalRadiologyCollisionControls(sofa::helper::vector<Real> &x_point_list,
                                        sofa::helper::vector<int> &id_instrument_list, sofa::helper::vector<int> &removeEdge)
{


    if(id_instrument_list.size() != x_point_list.size())
    {
                if ( this->f_printLog.getValue() )
                        serr<<" pb: interventionalRadiologyCollisionControls : the list do not have the same size !!"<<sendl;
        return;
    }



      // we enter the point from the tip to the end of the combined instrument

    // x_abs_curv provides the value of the curv abs along the combined instrument
    unsigned int node= nodeCurvAbs.size()-1;
    Real x_abs_curv = nodeCurvAbs[node];
    int first_instru_on_x = id_instrument_curvAbs_table[node][0];

    sofa::helper::vector<unsigned int> SegRemove;

    for (unsigned int it=0; it<m_instrumentsList.size(); it++)
        SegRemove.push_back(0);


    for (int i=x_point_list.size()-1; i>=0; i--)
    {

        //1.  we determin if the poin ument
        int instrument_id = id_instrument_list[i];

        // x_max for the instrument that is controlled (not dropped part)
        Real x_max_controlled = this->m_instrumentsList[instrument_id]->getRestTotalLength();

        if (x_point_list[i]>x_max_controlled)
        {
            int id_instr=id_instrument_list[i];
            SegRemove[id_instr] = i;
            continue;
        }



        // 2. we assign the value of the curv abs for the point and the corresponding instrument
        Real xtip_first_instru_on_x = this->xtip.getValue()[first_instru_on_x];
        Real x_begin = xtip_first_instru_on_x - m_instrumentsList[first_instru_on_x]->getRestTotalLength();
        x_point_list[i] = x_abs_curv - x_begin; // provides the "local" curv absc of the point (on the instrument reference)
        id_instrument_list[i] = first_instru_on_x;

        // 3. we look for the collision sampling of the current instrument in order to "place" the following point
        Real x_incr;
        this->m_instrumentsList[first_instru_on_x]->getCollisionSampling(x_incr, x_point_list[i]);
        x_abs_curv -= x_incr;

        // the following point could not have x_abs_curv<0;
        if (x_abs_curv<0.0)
        {
            x_abs_curv=0.0;
            continue;
        }

        // the following point can be place on an other instrument
        while (node > 0 && x_abs_curv < nodeCurvAbs[node-1])
        {
            node--; // we change the beam support...
            if( id_instrument_curvAbs_table[node][0] != first_instru_on_x)
            {
                // instrument has changed !!
                first_instru_on_x = id_instrument_curvAbs_table[node][0];
                x_abs_curv = nodeCurvAbs[node];
                break;
            }

        }

    }


    for (unsigned int it=0; it<m_instrumentsList.size(); it++)
    {
        if(SegRemove[it]!=0)
        {
            if (this->f_printLog.getValue())
                std::cout<<" ADD point " <<SegRemove[it]<<" in removeEdge list"<<std::endl;
            removeEdge.push_back(SegRemove[it]);
        }
    }




          // A  way to detect if a collision point is "activated" or not=> look at its curv_abs  and if > 0, it is active
          // first, we need to compute abs_curv_point
    sofa::helper::vector<Real> abs_curv_point;
    abs_curv_point.clear();

    for (unsigned int i=0; i<x_point_list.size(); i++)
    {

        int instrument_id = id_instrument_list[i];

        // x_max for the instrument that is controlled (not dropped part)
        Real x_max_instrument = this->m_instrumentsList[instrument_id]->getRestTotalLength();

        Real x_tip_instrument = xtip.getValue()[instrument_id];


        Real x_point= x_point_list[i] - x_max_instrument + x_tip_instrument;

        abs_curv_point.push_back( x_point );

    }


    // x point < epsilon... it is not activated`

    activated_Points_buf.clear();
    activated_Points_buf.push_back(false);
    for (unsigned int i=1; i<abs_curv_point.size(); i++)
    {

        Real x_max_instrument = this->m_instrumentsList[id_instrument_list[i]]->getRestTotalLength();

        if (abs_curv_point[i] < 0.00000001*x_max_instrument || fabs(abs_curv_point[i] - abs_curv_point[i-1])<0.00000001*x_max_instrument)
            activated_Points_buf.push_back(false);
        else
            activated_Points_buf.push_back(true);

        //std::cout<<" point "<<i<<" x_point = "<< x_point<< " instrument_id = "<<instrument_id<<" x_max_instrument ="<<x_max_instrument<<"  x_tip_instrument = "<<x_tip_instrument<<std::endl;

    }


}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::activateBeamListForCollision( sofa::helper::vector<Real> &curv_abs, sofa::helper::vector< sofa::helper::vector<int> > &id_instrument_table)
{

#ifdef DEBUG
    std::cout<<"  +++++++++++ \n id_instrument_table :"<<std::endl;
#endif
// 0. useful for rigidification
    helper::ReadAccessor< Data< helper::set< Real > > > rigidCurvAbs = m_rigidCurvAbs;

    unsigned int nbRigidAbs = rigidCurvAbs->size();

    if(curv_abs.size() != id_instrument_table.size())
        serr<<" id_instrument_table.size() is not equal to curv_abs.size()"<<sendl;


// 1. clear the information related to the collision for each instrument
   for (unsigned int instr=0; instr<m_instrumentsList.size(); instr++)
       m_instrumentsList[instr]->clearCollisionOnBeam();

// 2.   for each node along the structure, detect the "visible" instrument (the one with the largest section is supposed to be the first in the list)
//      if visible, the beam is assigned for collision on this instrument

    for (unsigned int p=0; p<id_instrument_table.size()-1; p++)
    {
        unsigned int instr0;
        if(id_instrument_table.size()==1 )
        {
            instr0=id_instrument_table[0][0];
            m_instrumentsList[ instr0 ]->addCollisionOnBeam(p);
        }
        else
        {
            instr0=id_instrument_table[p+1][0];// get the first instrument


        ////// beam p should be assigned to instrument num instr0

 //3 .  Before assignement, verification that the beam is not on a rigidified part !
            bool rigid=false;
            RealConstIterator it = rigidCurvAbs->begin();
#ifdef DEBUG
            std::cout<<"( *it) begin "<<(*it)<<" (*it) end"<< (* rigidCurvAbs->end())<<std::endl;
#endif
            while (it!=rigidCurvAbs->end())
            {
#ifdef DEBUG
                std::cout<<" curv_abs[p+1] =  "<<curv_abs[p+1]<<" (*it)"<<(*it)<<std::endl;
 #endif
                if(curv_abs[p+1] <= (*it) )
                {
                    break;
                }
                else
                {
#ifdef DEBUG
                    std::cout<<"++++++++++\n Rigidification detected "<<std::endl;
#endif
                    rigid = !rigid;
                    it++;
                }
            }
#ifdef DEBUG
            if (rigid)
            {
                std::cout<<" beam "<<p<<" is detected to be rigidified"<<std::endl;
            }
#endif
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



    /////////////// Create vectors with the CurvAbs of the noticiable points and the id of the corresponding instrument
    sofa::helper::vector<Real> newCurvAbs;

    //// in case of drop:
    unsigned int previousNumControlledNodes = numControlledNodes;
    unsigned int seg_remove = 0;


    if (this->dropCall)
    {
        processDrop(previousNumControlledNodes, seg_remove);
    }




    //// STEP 1
    ////find the total length of the COMBINED INSTRUMENTS and the one for which xtip > 0 (so the one which are simulated)
    Real xend;
    sofa::helper::AdvancedTimer::stepBegin("step1");
    Real totalLengthCombined=0.0;
//    int last_instrument=0; //commented to remove compilation warning
    sofa::helper::vector<Real> xbegin;
    for (unsigned int i=0; i<m_instrumentsList.size(); i++)
    {
         xend= this->xtip.getValue()[i];
        Real xb = xend - this->m_instrumentsList[i]->getRestTotalLength();
        xbegin.push_back(xb);

        if (xend> totalLengthCombined)
        {
            totalLengthCombined=xend;
//            last_instrument=i; //commented to remove compilation warning
        }

        // clear the present interpolation of the beams
        this->m_instrumentsList[i]->clear();

        if( xend > 0.0)
        {
            // create the first node (on x=0)
            newCurvAbs.push_back(0.0);

        }
    }

    //// SOME VERIF OF STEP 1
    // if the totalLength is 0, move the first instrument
    if(totalLengthCombined<0.0001)
    {
        sofa::helper::vector<Real> &x = (*this->xtip.beginEdit());
        x[0]=0.0001;
        this->xtip.endEdit();
        applyInterventionalRadiologyController();
        return;
    }
    sofa::helper::AdvancedTimer::stepEnd("step1");


    //// STEP 2:
    //// get the noticeable points that need to be simulated
    /////////////////////////
    // Fill=> newCurvAbs which provides a vector with curvilinear abscissa of each simulated node
    //     => id_instrument_table which provides for each simulated node, the id of all instruments which belong this node
    //     => xbegin (theoritical curv abs of the beginning point of the instrument (could be negative) xbegin= xtip - intrumentLength)
    sofa::helper::AdvancedTimer::stepBegin("step2");
    sofa::helper::vector< sofa::helper::vector<int> > id_instrument_table;
    this->interventionalRadiologyComputeSampling(newCurvAbs,id_instrument_table, xbegin, totalLengthCombined);
    sofa::helper::AdvancedTimer::stepEnd("step2");




    ////STEP 3
       // Re-interpolate the positions and the velocities
    /////////////////////////
    sofa::helper::AdvancedTimer::stepBegin("step3");
    unsigned int nbeam=newCurvAbs.size()-1; // number of simulated beams
    unsigned int nnode=newCurvAbs.size(); // number of simulated nodes



    unsigned int nnode_old= nodeCurvAbs.size(); // previous number of simulated nodes;

    //VecCoord& x = (*this->getMechanicalState()->read(sofa::core::ConstVecCoordId::position())->getValue());
    Data<VecCoord>* datax = this->getMechanicalState()->write(sofa::core::VecCoordId::position());
    VecCoord& x = *datax->beginEdit();

    VecCoord xbuf =x;



    sofa::helper::vector<Real> modifiedCurvAbs;

    this->totalLengthIsChanging(newCurvAbs, modifiedCurvAbs, id_instrument_table);
    Real xmax_prev = nodeCurvAbs[nodeCurvAbs.size()-1];

    for (unsigned int p=0; p<nbeam+1; p++)
    {
        int idP = numControlledNodes-nnode + p;
        Real xabs = modifiedCurvAbs[p];
        if (this->f_printLog.getValue())
            sout <<" interpolate point at abs "<<xabs<< sendl;

        // 2 cases:  TODO : remove first case
            //1. the abs curv is further than the previous state of the instrument
            //2. this is not the case and the node position can be interpolated using previous step positions
        if(xabs > xmax_prev + threshold.getValue())
        {
            if (this->f_printLog.getValue()){
                serr<<"case 1 should never happen ==> avoid using totalLengthIsChanging ! xabs = "<<xabs<<" - xmax_prev = "<<xmax_prev<<sendl;
                serr<<"newCurvAbs  = "<<newCurvAbs<<"  previous nodeCurvAbs"<<nodeCurvAbs<<sendl;
                serr<<"modifiedCurvAbs ="<<modifiedCurvAbs<<sendl;
            }
            // case 1 (the abs curv is further than the previous state of the instrument)
            ////// verifier qu'il s'agit bien d'un instrument qu'on est en train de controller
            ////// interpoler toutes les positions "sorties" de l'instrument en supprimant l'ajout de dx qu'on vient de faire
        }
        else
        {
            // case 2 (the node position can be interpolated straightfully using previous step positions)
            unsigned int p0=0;
            while(p0<this->nodeCurvAbs.size())
            {
                if((nodeCurvAbs[p0]+threshold.getValue())>xabs)
                    break;
                p0++;
            }

            int idP0 =  previousNumControlledNodes + seg_remove - nnode_old + p0 ;

            if(fabs(nodeCurvAbs[p0]-xabs)<threshold.getValue())
            {
                x[idP] = xbuf[idP0];
            }
            else
            {

                // the node must be interpolated using beam interpolation
                    //find the instrument
                int id = id_instrument_curvAbs_table[p0][0];
                //find the good beam (TODO: do not work if xbegin of one instrument >0)
                int b = p0-1;
                    // test to avoid wrong indices
                if (b<0)
                {
                    x[p]=this->startingPos.getValue();
                }
                else
                {
                    // the position is interpolated
                    Transform global_H_interpol;
                    Real ratio = (xabs - nodeCurvAbs[b])/ (nodeCurvAbs[p0]-nodeCurvAbs[b]);
                    unsigned int numBeamNotUnderControl = this->m_instrumentsList[id]->getNumBeamsNotUnderControl();
                    if (this->f_printLog.getValue())
                        std::cout<<" b= "<<b<<" / numBeamNotUnderControl= "<<numBeamNotUnderControl<<std::endl;

                    //unsigned int beamID = b + numBeamNotUnderControl;

                    Transform Global_H_local0(xbuf[idP0-1].getCenter(),xbuf[idP0-1].getOrientation() ), Global_H_local1(xbuf[idP0].getCenter(),xbuf[idP0].getOrientation() );

                    Real L = nodeCurvAbs[p0] - nodeCurvAbs[b];


                    this->m_instrumentsList[id]->InterpolateTransformUsingSpline(global_H_interpol, ratio, Global_H_local0, Global_H_local1 ,L);

                    x[idP].getCenter() = global_H_interpol.getOrigin();
                    x[idP].getOrientation() = global_H_interpol.getOrientation();

                    // TODO: interpolate velocity !!

                    if (this->f_printLog.getValue())
                        std::cout<<" le noeud ["<<idP<<"] est interpolé"<<"  x["<<idP<<"]="<<x[idP]<<std::endl;

                }
            }

        }

    }
    sofa::helper::AdvancedTimer::stepEnd("step3");


 ////STEP 4
 // Assign the beams
 /////////////////////////
    //const sofa::helper::vector<sofa::core::topology::BaseMeshTopology::Edge>* edgeTopo = this->_topology->getEdges();
    // TODO: verifier que les edges sont bien "ordonnés"
    sofa::helper::AdvancedTimer::stepBegin("step4");
    unsigned int numEdges= numControlledNodes-1;

    // verify that there is a sufficient number of Edge in the topology : TODO if not, modify topo !
    if(numEdges<nbeam)
    {
        if (this->f_printLog.getValue())
        {
            serr<<" not enough edges in the topology !!"<<sendl;
        }
        nbeam=numEdges;
    }


    for (unsigned int b=0; b<nbeam; b++)
    {
        Real x0 = newCurvAbs[b];
        Real x1 = newCurvAbs[b+1];
        for (unsigned int i=0; i<m_instrumentsList.size(); i++)
        {
            Real xmax = this->xtip.getValue()[i];
            Real xmin = xbegin[i];

            Real eps= threshold.getValue();

            if (x0>(xmin-eps) && x0<(xmax+eps) && x1>(xmin-eps) && x1<(xmax+eps))
            {
                sofa::core::topology::BaseMeshTopology::EdgeID eID = (sofa::core::topology::BaseMeshTopology::EdgeID)(numEdges-nbeam + b );

                Real length = x1 - x0;
                Real x0_local = x0-xmin;
                Real x1_local = x1-xmin;

                Real theta = this->rotationInstrument.getValue()[i];

                this->m_instrumentsList[i]->addBeam(eID, length, x0_local, x1_local,theta );

            }
        }
    }
    sofa::helper::AdvancedTimer::stepEnd("step4");

    ////STEP 5
    // Fix the not simulated nodes
    /////////////////////////
    sofa::helper::AdvancedTimer::stepBegin("step5");
    unsigned int first_simulated_node = numControlledNodes - nbeam;

    //1. Fix the nodes (beginning of the instruments) that are not "out"
    fixFirstNodesWithUntil(first_simulated_node);

    //2. Fix the node that are "fixed"
    // When there are rigid segments, # of dofs is different than # of edges and beams
    helper::ReadAccessor< Data< helper::set< Real > > > rigidCurvAbs = m_rigidCurvAbs;

    bool rigid=false;
    unsigned int nbRigidAbs = rigidCurvAbs->size();
    if (nbRigidAbs>0 && (nbRigidAbs%2)==0)
    {
        RealConstIterator it = rigidCurvAbs->begin();

        for (unsigned int i=0; i<newCurvAbs.size(); i++)
        {
            std::cout<<"it "<<(*it)<<std::endl;
            if (newCurvAbs[i] < ((*it)+0.001) && newCurvAbs[i] > ((*it)-0.001)) // node= border of the rigid segment
            {
                std::cout<<" border dtected"<<std::endl;

                std::cout<<"  -  abs_curv ="<<  newCurvAbs[i]<<" i="<<i<<" rigid="<<rigid<<std::endl;

                if (!rigid)
                {
                    rigid=true;
                    _fixedConstraint->addConstraint(first_simulated_node+i);
                }
                else
                {
                    _fixedConstraint->addConstraint(first_simulated_node+i);
                    rigid=false;
                }
                it++;

                std::cout<<" rigidCurvAbs->end() ="<<(*rigidCurvAbs->end())<<" it="<<(*it)<<std::endl;
                if(it==rigidCurvAbs->end())
                    break;


            }
            else
            {
                std::cout<<"abs_curv ="<<  newCurvAbs[i]<<" i="<<i<<" rigid="<<rigid<<std::endl;
                if(rigid)
                    _fixedConstraint->addConstraint(first_simulated_node+i);
            }

        }

    }


    ////STEP 6
    // Activate Beam List for collision of each instrument
    /////////////////////////
    activateBeamListForCollision(newCurvAbs,id_instrument_table);







    ///////////// BUF id_instrument_curvAbs_table and nodeCurvAbs

    nodeCurvAbs = newCurvAbs;
    id_instrument_curvAbs_table = id_instrument_table;
    sofa::helper::AdvancedTimer::stepEnd("step5");

    if (this->f_printLog.getValue()){
        std::cout<<"id_instrument_curvAbs_table = \n";
        for (unsigned int i=0; i<id_instrument_curvAbs_table.size();i++)
            std::cout<<id_instrument_curvAbs_table[i]<<" end control"<<std::endl;
    }








    datax->endEdit();
}

template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::totalLengthIsChanging(const sofa::helper::vector<Real>& newNodeCurvAbs, sofa::helper::vector<Real>& modifiedNodeCurvAbs, const sofa::helper::vector< sofa::helper::vector<int> >& newTable)
{

    // If total length is changing, it means that we need to simulate the fact that some beam will "move" during the time step.
    // thus during the interpolation process, where the position of the nodes is based on the previous position of the wire,
    // we initialize some points at a x_curv ref pos without the motion (computed by DLength)
    // due to the elasticity of the beam, the point will then naturally go the position that reespects the newNodeCurvAbs

    int id = controlledInstrument.getValue();
    Real DLength = newNodeCurvAbs[ newNodeCurvAbs.size()-1] - nodeCurvAbs[nodeCurvAbs.size() - 1];
    modifiedNodeCurvAbs = newNodeCurvAbs;



    // we look for the last value in the CurvAbs
    if(fabs(DLength) > threshold.getValue())
    {

        if(newTable[newNodeCurvAbs.size()-1][0] != id)
        {
            //this->f_printLog.setValue(true);

        }

        unsigned int i=newTable.size()-1;
        while (i>0 && newTable[i].size()==1)
        {
            modifiedNodeCurvAbs[i]-=DLength;
            i--;

        }

    }

    /// \todo : be careful !! modifiedNodeCurvAbs is no more necessarily well sorted !!!



}

// sort the curv Abs in the ascending order and avoid doubloon

template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::sortCurvAbs(sofa::helper::vector<Real> &CurvAbs, sofa::helper::vector< sofa::helper::vector<int> >& id_instrument_table)
{
    sofa::helper::vector<Real> newCurvAbs;

    Real eps = threshold.getValue();


    newCurvAbs.clear();

    while(CurvAbs.size()>0)
    {
        Real xmin = 1.0e30;
        unsigned int index_min=0;
        for(unsigned int j=0; j<CurvAbs.size(); j++)
        {
            if(xmin > CurvAbs[j] )
            {
                xmin = CurvAbs[j];
               index_min=j;
            }
        }

        sofa::helper::removeIndex( CurvAbs, index_min);

        // verify using a eps that the same node does not already exist
        Real x_last;
        if (newCurvAbs.size()>0 )
            x_last= newCurvAbs[newCurvAbs.size()-1];
        else
            x_last=-1.0e30;

        if( fabs(xmin - x_last) > eps)
        {
            newCurvAbs.push_back(xmin);
        }


    }

    CurvAbs = newCurvAbs;


    id_instrument_table.clear();

    for (unsigned int i=0; i<newCurvAbs.size(); i++)
    {
        sofa::helper::vector<int> listNewNode;
        for (unsigned int id=0; id<m_instrumentsList.size(); id++)
        {
            Real xend= this->xtip.getValue()[id];
            Real xbegin = xend - this->m_instrumentsList[id]->getRestTotalLength();

            if ( xbegin<(newCurvAbs[i]+eps) && xend >(newCurvAbs[i]-eps) )
            {
                listNewNode.push_back(id);

            }
        }

        if(listNewNode.size() ==0)
        {
            std::cerr<<" \n \n ERROR: no instrument find for curvAbs"<<newCurvAbs[i]<<std::endl;
            Real xend= this->xtip.getValue()[0];
            Real xbegin = xend - this->m_instrumentsList[0]->getRestTotalLength();

            std::cerr<<"Test for instrument [0] xend="<<xend<<"   xbegin"<<xbegin<<std::endl;
        }




        id_instrument_table.push_back(listNewNode);

    }
}


template <class DataTypes>
void InterventionalRadiologyController<DataTypes>::fixFirstNodesWithUntil(unsigned int first_simulated_Node)
{
    //VecCoord& xMstate = (*this->getMechanicalState()->read(sofa::core::ConstVecCoordId::position())->getValue());
    //VecDeriv& vMstate = (*this->getMechanicalState()->getV());

    helper::WriteAccessor<Data<VecCoord> > xMstate = *this->getMechanicalState()->write(sofa::core::VecCoordId::position());
    helper::WriteAccessor<Data<VecDeriv> > vMstate = *this->getMechanicalState()->write(sofa::core::VecDerivId::velocity());


    // set the position to startingPos for all the nodes that are not simulated
    // and add a fixedConstraint
    _fixedConstraint->clearConstraints();
    for(unsigned int i=0; i<first_simulated_Node-1 ; i++)
    {
        xMstate[i]=startingPos.getValue();
        vMstate[i].clear();
        _fixedConstraint->addConstraint(i);
    }
    indexFirstNode = first_simulated_Node-1 ;

}






} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGYCONTROLLER_INL */
