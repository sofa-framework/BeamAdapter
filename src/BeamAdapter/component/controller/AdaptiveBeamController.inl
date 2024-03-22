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
// C++ Implementation : AdaptiveBeamController
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

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
#include <sofa/component/topology/container/dynamic/EdgeSetGeometryAlgorithms.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <BeamAdapter/component/controller/AdaptiveBeamController.h>


namespace sofa::component::controller::_adaptivebeamcontroller_
{

using sofa::core::objectmodel::BaseContext ;
using std::string;
using sofa::core::VecCoordId;
using sofa::core::VecDerivId;

//TODO(dmarchal 2017-05-17) to christian & euallie.
//  This component seems specific to radiology instrument... but its name suggest a very generic behavior
// (remove in 1 year if not answered).
template <class DataTypes>
AdaptiveBeamController<DataTypes>::AdaptiveBeamController()
    //TODO(dmarchal 2017-05-17) a Data<vector of string> is a weird beast. Have to write test to validate behavior (remove in 1 year if not done).
    : d_interpolationPath(initData(&d_interpolationPath,"interpolation", "Path to the Interpolation component on scene"))
    , d_controlledInstrument(initData(&d_controlledInstrument, 0, "controlledInstrument", "provide the id of the interventional radiology instrument which is under control: press contr + number to change it"))
    , d_xtip(initData(&d_xtip,"xtip", "curvilinear abscissa of the tip of each interventional radiology instrument"))
    , d_rotationInstrument(initData(&d_rotationInstrument,"rotationInstrument", "angle of rotation for each interventional radiology instrument"))
    , d_step(initData(&d_step,(Real)0.1,"step","base step when changing beam length"))
    , d_angularStep(initData(&d_angularStep,(Real)(3.1416/20.0),"angularStep","base step when changing beam angle"))
    , d_speed(initData(&d_speed,(Real)0.0,"speed","continuous beam length increase/decrease"))
{
}

template <class DataTypes>
void AdaptiveBeamController<DataTypes>::init()
{
    //////////////////// Initialize internal data structure. ///////////////////////////////////////
    m_xAbsCollisionPointsBuffer.clear();

    //////////////////// Handle the component's attributes /////////////////////////////////////////
    BaseContext* c = this->getContext();
    this->f_listening.setValue(true);

    const type::vector<string>& interpolName = d_interpolationPath.getValue();
    if (interpolName.empty()) {
        m_adaptiveinterpolation = c->get<BInterpolation>(BaseContext::Local);
    } else {
        m_adaptiveinterpolation = c->get<BInterpolation>(d_interpolationPath.getValue()[0]);
    }

    if(m_adaptiveinterpolation==nullptr){
        msg_error() <<"No Beam Interpolation component found. This component is thus de-actiated. ";
    }
    else {
        msg_info() <<"This component operates on '"<<m_adaptiveinterpolation->getName()<<"'." ;
    }

    if(d_speed.getValue()>0)
    {
        FF=true;
        RW=false;
    }


    Inherit::init();

    //TODO(dmarchal 2017-05-17) Please specify who will do that and when this will be done.
    // TODO : VERIFIER QUE LA TOPOLOGIE EST CELLE D'UN WIRE
}


template <class DataTypes>
void AdaptiveBeamController<DataTypes>::reinit()
{
    AdaptiveBeamController<DataTypes>::applyController();
}


template <class DataTypes>
void AdaptiveBeamController<DataTypes>::onMouseEvent(MouseEvent *mev)
{
    /// Retrieve the mouse position
    int PosX = mev->getPosX();
    int PosY = mev->getPosY();

    /// Translation input
    Real PosYcorr = 0.0;
    auto id = d_controlledInstrument.getValue();
    auto x_instr_tip = sofa::helper::getWriteOnlyAccessor(d_xtip);

    if (id >= (int)(x_instr_tip.size())){
        msg_warning() << "The instrument number " << id << " does not exist (size =" << x_instr_tip.size() << ") exiting mouse event.";
        return;
    }

    PosYcorr = -PosY * 0.2;
    x_instr_tip[id] += PosYcorr;

    /// Rotation input
    auto rot_instrument = sofa::helper::getWriteOnlyAccessor(d_rotationInstrument);
    Real PosXcorr = PosX * 0.015;
    rot_instrument[id] += PosXcorr;
}

//TODO(dmarchal 2017-05-17) Christian & Eulalie, how can user know which key-mouse behavior he
// can expect ? A kind a of getHelp function may be usefull for controllers ?
// (remove in 1 year if not answered)
template <class DataTypes>
void AdaptiveBeamController<DataTypes>::onKeyPressedEvent(KeypressedEvent *kev)
{
    ///////////////////////////////// Control keys for interventonal Radiology simulations:
    switch(kev->getKey())
    {
    case '0':
        this->d_controlledInstrument.setValue(0);
        break;

    case 'A':
    {
        int id = d_controlledInstrument.getValue();
        auto rot_instrument = sofa::helper::getWriteOnlyAccessor(d_rotationInstrument);
        rot_instrument[id] += d_angularStep.getValue();
    }
        break;
    case 'E':
    {

        int id = d_controlledInstrument.getValue();
        auto rot_instrument = sofa::helper::getWriteOnlyAccessor(d_rotationInstrument);
        rot_instrument[id] -= d_angularStep.getValue();

    }
        break;
    case '+':
    {
        int id = d_controlledInstrument.getValue();
        auto x_instr_tip = sofa::helper::getWriteOnlyAccessor(d_xtip);
        if (id >= (int)x_instr_tip.size()){
            msg_warning() << "Controlled Instument num "<<id<<" do not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead" ;
            id=0;
        }
        x_instr_tip[id] += d_step.getValue();
    }
        break;

    case '-':
    {
        int id = d_controlledInstrument.getValue();
        auto x_instr_tip = sofa::helper::getWriteOnlyAccessor(d_xtip);
        if (id >= (int)x_instr_tip.size()){
            msg_warning() << "Controlled Instument num "<<id<<" do not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead" ;
            id=0;
        }
        x_instr_tip[id] -= d_step.getValue();
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
}

template <class DataTypes>
void AdaptiveBeamController<DataTypes>::onBeginAnimationStep(const double /*dt*/)
{
    applyController();
}


template <class DataTypes>
void AdaptiveBeamController<DataTypes>::applyController()
{
    Data<VecCoord>* datax = this->getMechanicalState()->write(VecCoordId::position());
    Data<VecDeriv>* datav = this->getMechanicalState()->write(VecDerivId::velocity());
    auto x = sofa::helper::getWriteOnlyAccessor(*datax);
    auto v = sofa::helper::getWriteOnlyAccessor(*datav);

    /////// analyse de la beam actuelle :: TODO => use nodeCurvAbs which store this info
    Real totalLength=0.0;
    type::vector<Real> oldCurvAbs, newCurvAbs;
    oldCurvAbs.push_back(0.0);
    unsigned int numBeams = m_adaptiveinterpolation->getNumBeams();
    for (unsigned int b=0; b<numBeams; b++)
    {
        totalLength += m_adaptiveinterpolation->getLength(b);
        oldCurvAbs.push_back(totalLength);
    }

    Real totalSplineLength=0.0;
    type::vector<Real> splineAbs;
    Vec3 P0,P1,P2,P3;
    splineAbs.push_back(0.0);
    for (unsigned int b=0; b<numBeams; b++)
    {
        m_adaptiveinterpolation->getSplinePoints(b,x,P0,P1,P2,P3);
        Vec3 P0P1,P1P2,P2P3;
        P0P1=P1-P0; P1P2=P2-P1; P2P3=P3-P2;
        Real l =  P0P1.norm() + P1P2.norm() + P2P3.norm();
        l+= 30.0*fabs(l-m_adaptiveinterpolation->getLength(b));

        totalSplineLength += l;
        splineAbs.push_back(totalSplineLength);
    }

    Real samplingSpline = totalSplineLength/numBeams;
    Real x_spline=0;
    unsigned int j=0;
    newCurvAbs.push_back(0.0);
    for (unsigned int b=0; b<numBeams-1; b++)
    {
        x_spline+=samplingSpline;
        while(x_spline>splineAbs[j])
            j++;

        // x_spline est entre splineAbs[j-1] et splineAbs[j]
        Real ratio = (x_spline - splineAbs[j-1]) / (splineAbs[j] - splineAbs[j-1]);

        if (ratio<0 || ratio >1 )
            msg_warning() <<"Ratio = "<<ratio<<" while it should be between 0 and 1 " ;

        Real x_curv= oldCurvAbs[j-1] + ratio * (oldCurvAbs[j]  - oldCurvAbs[j-1]);
        newCurvAbs.push_back( x_curv );
    }
    // the last CurvAbs is known:
    newCurvAbs.push_back(  totalLength );

    /////////// Meme modifications a partir de  oldCurvAbs et newCurvAbs
    ///(rmq : pour l'instant, on ne s'occupe pas des 2 extremités du fil)
    j=0;
    for (unsigned int i=1; i<newCurvAbs.size()-1; i++)
    {
        while(newCurvAbs[i]>oldCurvAbs[j])
        {
            j++;
            if (j>=oldCurvAbs.size()) // DEBUG //
            {
                dmsg_warning() <<"j ="<<j<<">=oldCurvAbs.size()" ;
                return;
            }
        }

        Real L = m_adaptiveinterpolation->getLength(j-1);
        Real L0 = newCurvAbs[i] - oldCurvAbs[j-1];
        Real ratio=L0/L;

        Transform global_H_interpol;
        Vec3 null(0,0,0);

        // ATTENTION !!!!! =>  si j<i on utilise une position et une vitesse que l'on vient de modifier !
        // Il faudrait stocker avant de modifier !
        m_adaptiveinterpolation->InterpolateTransformUsingSpline(j-1,ratio,null,x,global_H_interpol);
        x[i].getCenter() = global_H_interpol.getOrigin();
        x[i].getOrientation() = global_H_interpol.getOrientation();
        v[i] = v[j-1] * (1-ratio)  + v[j] * ratio;

        // on definie la longeur des beam en fonction des nouvelles abscisses curvilignes:
        L =  newCurvAbs[i] - newCurvAbs[i-1];
        m_adaptiveinterpolation->setLength(i-1,L);
    }
    //on definie la longueur de la dernière beam:
    Real L = newCurvAbs[newCurvAbs.size()-1] - newCurvAbs[newCurvAbs.size()-2];
    m_adaptiveinterpolation->setLength(newCurvAbs.size()-2,L);
}

} // namespace sofa::component::controller::_adaptivebeamcontroller_


