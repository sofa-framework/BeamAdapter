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

#ifndef SOFA_COMPONENT_CONTROLLER_ADAPTIVEBEAMCONTROLLER_INL
#define SOFA_COMPONENT_CONTROLLER_ADAPTIVEBEAMCONTROLLER_INL

#include "AdaptiveBeamController.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/objectmodel/MouseEvent.h>
#include <sofa/core/objectmodel/KeypressedEvent.h>
//#include <BaseTopology/EdgeSetGeometryAlgorithms.h>
#include <sofa/component/topology/EdgeSetGeometryAlgorithms.h>
#include <sofa/core/behavior/MechanicalState.h>



namespace sofa
{

namespace component
{

namespace controller
{

template <class DataTypes>
AdaptiveBeamController<DataTypes>::AdaptiveBeamController()
: m_adaptiveinterpolation(NULL)
, m_interpolationPath(initData(&m_interpolationPath,"interpolation", "Path to the Interpolation component on scene"))
, controlledInstrument(initData(&controlledInstrument, 0, "controlledInstrument", "provide the id of the interventional radiology instrument which is under control: press contr + number to change it"))
, xtip(initData(&xtip,"xtip", "curvilinear abscissa of the tip of each interventional radiology instrument"))
, rotationInstrument(initData(&rotationInstrument,"rotationInstrument", "angle of rotation for each interventional radiology instrument"))
, step(initData(&step,(Real)0.1,"step","base step when changing beam length"))
, angularStep(initData(&angularStep,(Real)(3.1416/20.0),"angularStep","base step when changing beam angle"))
, speed(initData(&speed,(Real)0.0,"speed","continuous beam length increase/decrease"))
, startingPos(initData(&startingPos,Coord(),"startingPos","starting pos for inserting the instrument"))
, threshold(initData(&threshold, (Real)0.000001, "threshold", "threshold for controller precision which is homogeneous to the unit of length"))
{
	_fixedConstraint = NULL;
	this->dropCall = false;
}



template <class DataTypes>
void AdaptiveBeamController<DataTypes>::init()
{


    ///////// get the Adaptive Interpolation component ///////
    //std::vector<sofa::core::behavior::LinearSolver*> solvers;
    core::objectmodel::BaseContext * c = this->getContext();

    const helper::vector<std::string>& interpolName = m_interpolationPath.getValue();
    if (interpolName.empty()) {
        m_adaptiveinterpolation= c->get<BInterpolation>(core::objectmodel::BaseContext::Local);
    } else {

        m_adaptiveinterpolation = c->get<BInterpolation>(m_interpolationPath.getValue()[0]);
    }

    if(m_adaptiveinterpolation==NULL)
        serr<<" no Beam Interpolation found !!! the component can not work"<<sendl;
    else
        sout<<" interpolation named"<<m_adaptiveinterpolation->getName()<<" found (for "<<this->getName()<<")"<<sendl;


    //////// INIT:

	xAbs_collisionPoints_buf.clear();

	if(speed.getValue()>0)
	{
		FF=true;
		RW=false;
	}

	_topology = this->getContext()->getMeshTopology();
	this->f_listening.setValue(true);

	Inherit::init();


	// TODO : VERIFIER QUE LA TOPOLOGIE EST CELLE D'UN WIRE

}


template <class DataTypes>
void AdaptiveBeamController<DataTypes>::reinit()
{
	applyController();

}


template <class DataTypes>
void AdaptiveBeamController<DataTypes>::onMouseEvent(core::objectmodel::MouseEvent *mev)
{
	int PosX = mev->getPosX();
	int PosY = mev->getPosY();


	//Translation input
	Real PosYcorr = 0.0;
	int idy = controlledInstrument.getValue();
	sofa::helper::vector<Real> &x_instr_tip = (*this->xtip.beginEdit());
	if (idy >= (int)x_instr_tip.size()){
		std::cerr<<"WARNING controlled Instument num "<<idy<<" do not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead"<<std::endl;
		idy=0;
	}
	PosYcorr = -PosY*0.2;
	x_instr_tip[idy] += PosYcorr;
	this->xtip.endEdit();



	//Rotation input
	Real PosXcorr = 0.0;
	int idx = controlledInstrument.getValue();
	sofa::helper::vector<Real> &rot_instrument = (*this->rotationInstrument.beginEdit());
	PosXcorr = PosX*0.015;
	rot_instrument[idx] += PosXcorr;

}


template <class DataTypes>
void AdaptiveBeamController<DataTypes>::onKeyPressedEvent(core::objectmodel::KeypressedEvent *kev)
{

	///////////////////////////////// Control keys for interventonal Radiology simulations:
	switch(kev->getKey())
	{
	case 'D':
		this->dropCall = true;
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
			std::cerr<<"WARNING controlled Instument num "<<id<<" do not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead"<<std::endl;
			id=0;
		}
		x_instr_tip[id] += step.getValue();
		this->xtip.endEdit();
	}
	break;

	case '-':
	{
		int id = controlledInstrument.getValue();
		sofa::helper::vector<Real> &x_instr_tip = (*this->xtip.beginEdit());
		if (id >= (int)x_instr_tip.size()){
			std::cerr<<"WARNING controlled Instument num "<<id<<" do not exist (size ="<< x_instr_tip.size() <<") use instrument 0 instead"<<std::endl;
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

}


template <class DataTypes>
void AdaptiveBeamController<DataTypes>::onBeginAnimationStep(const double /*dt*/)
{
	applyController();
}


template <class DataTypes>
void AdaptiveBeamController<DataTypes>::applyController()
{

	Data<VecCoord>* datax = this->getMechanicalState()->write(sofa::core::VecCoordId::position());
	Data<VecDeriv>* datav = this->getMechanicalState()->write(sofa::core::VecDerivId::velocity());
	VecCoord& x = *datax->beginEdit();
	VecDeriv& v = *datav->beginEdit();

	/////// analyse de la beam actuelle :: TODO => use nodeCurvAbs which store this info
	Real totalLength=0.0;
	sofa::helper::vector<Real> oldCurvAbs, newCurvAbs;
	oldCurvAbs.push_back(0.0);
	unsigned int numBeams = m_adaptiveinterpolation->getNumBeams();
	for (unsigned int b=0; b<numBeams; b++)
	{
		totalLength += m_adaptiveinterpolation->getLength(b);
		oldCurvAbs.push_back(totalLength);
	}

	Real totalSplineLength=0.0;
	sofa::helper::vector<Real> splineAbs;
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
		                                                	  std::cerr<<"WARNING ratio = "<<ratio<<" while it should be between 0 and 1 "<<std::endl;

		                                                  Real x_curv= oldCurvAbs[j-1] + ratio * (oldCurvAbs[j]  - oldCurvAbs[j-1]);
		                                                  newCurvAbs.push_back( x_curv );
	}
	// the last CurvAbs is known:
	newCurvAbs.push_back(  totalLength );


	//newCurvAbs = oldCurvAbs;
	//unsigned int j;

	/////////// Meme modifications a partir de  oldCurvAbs et newCurvAbs
	//(rmq : pour l'instant, on ne s'occupe pas des 2 extremités du fil)
	j=0;
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

	datax->endEdit();
	datav->endEdit();

}





} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_ADAPTIVEBEAMCONTROLLER_INL */
