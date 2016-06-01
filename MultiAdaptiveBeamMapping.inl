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
// C++ Implementation : UnifiedMultiMultiAdaptiveBeamMapping
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_MAPPING_MULTIADAPTIVEBEAMMAPPING_INL
#define SOFA_COMPONENT_MAPPING_MULTIADAPTIVEBEAMMAPPING_INL

#include "MultiAdaptiveBeamMapping.h"
#include <sofa/core/behavior/MechanicalState.h>
#include <string>
#include <sofa/core/Mapping.inl>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/simulation/AnimateBeginEvent.h>


namespace sofa
{

namespace component
{

namespace mapping
{

using namespace sofa::defaulttype;



template <class TIn, class TOut>
MultiAdaptiveBeamMapping< TIn, TOut>::MultiAdaptiveBeamMapping(core::State< In >* from, core::State< Out >* to,InterventionalRadiologyController<TIn>* _ircontroller)
: Inherit(from, to)
, useCurvAbs(initData(&useCurvAbs,true,"useCurvAbs","true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same"))
, m_controlerPath(initData(&m_controlerPath,"ircontroller", "Path to the ircontroller component on scene"))
, m_ircontroller(_ircontroller)
, isBarycentricMapping(false)
{
	this->addAlias(&m_controlerPath, "controller");
}


template <class TIn, class TOut>
MultiAdaptiveBeamMapping< TIn, TOut>::MultiAdaptiveBeamMapping()
: Inherit()
, useCurvAbs(initData(&useCurvAbs,true,"useCurvAbs","true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same"))
, m_controlerPath(initData(&m_controlerPath,"ircontroller", "Path to the ircontroller component on scene"))
, m_ircontroller(NULL)
, isBarycentricMapping(false)
{
	this->addAlias(&m_controlerPath, "controller");
}


template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::apply(const core::MechanicalParams* mparams /* PARAMS FIRST */, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
	VecCoord& out = *dOut.beginEdit();
	if(!isBarycentricMapping)
		out.resize(_xPointList.size());
	else if((int )out.size() != this->getMechTo()[0]->getSize())
		out.resize(this->getMechTo()[0]->getSize());
	dOut.endEdit();
	for (unsigned int subMap=0; subMap<m_subMappingList.size(); subMap++)
	{
		if (this->f_printLog.getValue()){
			std::cout<<"apply ";
			m_subMappingList[subMap]->printIstrumentInfo();//ctn_DEV
		}
		m_subMappingList[subMap]->apply(mparams /* PARAMS FIRST */, dOut, dIn);
	}

}

template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::applyJ(const core::MechanicalParams* mparams /* PARAMS FIRST */, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
	for (unsigned int subMap=0; subMap<m_subMappingList.size(); subMap++)
	{
		if (this->f_printLog.getValue()){
			std::cout<<"applyJ ";
			m_subMappingList[subMap]->printIstrumentInfo();//ctn_DEV
		}
		m_subMappingList[subMap]->applyJ(mparams /* PARAMS FIRST */, dOut, dIn);
	}
}

template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::applyJT(const core::MechanicalParams* mparams /* PARAMS FIRST */, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{
	for (unsigned int subMap=0; subMap<m_subMappingList.size(); subMap++)
	{
		if (this->f_printLog.getValue()){
			std::cout<<"applyJT ";
			m_subMappingList[subMap]->printIstrumentInfo();//ctn_DEV
		}
		m_subMappingList[subMap]->applyJT(mparams /* PARAMS FIRST */, dOut, dIn);
	}
}


// MultiAdaptiveBeamMapping::applyJT(InMatrixDeriv& out, const OutMatrixDeriv& in) //
// this function propagate the constraint through the Adaptive Beam mapping :
// if one constraint along (vector n) with a value (v) is applied on the childModel (like collision model)
// then this constraint is transformed by (Jt.n) with value (v) for the rigid model
// note : the value v is not propagated through the mapping
template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::applyJT(const core::ConstraintParams* cparams /* PARAMS FIRST */, Data<InMatrixDeriv>& dOut, const Data<OutMatrixDeriv>& dIn)
{
	for (unsigned int subMap=0; subMap<m_subMappingList.size(); subMap++)
	{
		if (this->f_printLog.getValue()){
			std::cout<<"applyJTT ";
			m_subMappingList[subMap]->printIstrumentInfo();//ctn_DEV
		}
		m_subMappingList[subMap]->applyJT(cparams /* PARAMS FIRST */, dOut, dIn);
	}
}




template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::handleEvent(sofa::core::objectmodel::Event * event)
{

	if (dynamic_cast<sofa::simulation::AnimateBeginEvent *>(event))
	{
		assignSubMappingFromControllerInfo();
	}
}



template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::assignSubMappingFromControllerInfo()
{
	sofa::helper::vector<int> removeEdgeAtPoint;

	// 1. get the new controls
 //   std::cout<<" _xPointList before "<<_xPointList<<std::endl;
	m_ircontroller->interventionalRadiologyCollisionControls(_xPointList, _idm_instrumentList, removeEdgeAtPoint);
 //   std::cout<<" _xPointList after "<<_xPointList<<std::endl;

	if(!isBarycentricMapping)
	{
		//Case if this is not a barycentric mapping
        // 2. assign each value of xPointList to the corresponding "sub Mapping"
        // i.e = each point from xPointList of the collision model is controlled by a given instrument wich id is provided in _id_m_instrumentList
		sofa::helper::vector< sofa::helper::vector < Vec3 >* > pointsList;
		pointsList.clear();
		pointsList.resize(m_subMappingList.size());

		for (unsigned int i=0; i<m_subMappingList.size(); i++)
		{
			pointsList[i] = m_subMappingList[i]->points.beginEdit();
			pointsList[i]->clear();
			m_subMappingList[i]->clearidPointSubMap();
		}

		for (unsigned int i=0; i< _xPointList.size(); i++)
		{
			unsigned int  id = _idm_instrumentList[i];
			pointsList[ id ]->push_back( Vec3( _xPointList[i], 0, 0) );
			m_subMappingList[id]->addidPointSubMap(i);
		}
		for (unsigned int i=0; i<m_subMappingList.size(); i++)
		{
			m_subMappingList[i]->points.endEdit();
		}

		// handle the possible topological change

		sofa::helper::vector<unsigned int> edgeToRemove;
		edgeToRemove.clear();

		for (unsigned int i=0; i<removeEdgeAtPoint.size(); i++)
		{
			// look if the edge is already suppressed (only one edge around the targetted point)
			sofa::helper::vector< unsigned int > baseEdge(0);
			baseEdge = _topology->getEdgesAroundVertex((sofa::core::topology::BaseMeshTopology::PointID) removeEdgeAtPoint[i]);

			// need to be removed
			if (baseEdge.size() == 2)
			{
				// chose the edge to remove
				for (unsigned i=0; i<baseEdge.size(); i++)
				{
					const sofa::core::topology::BaseMeshTopology::Edge e = _topology->getEdge(baseEdge[i]);
					if( ((int) e[1]== removeEdgeAtPoint[i] && e[1] > e[0]) || ((int) e[0]==removeEdgeAtPoint[i] && e[0] > e[1]) )
					{
						edgeToRemove.push_back(baseEdge[i]);
					}
				}
			}
			else if (baseEdge.size() == 1)
				std::cout<<" ok, the edge is already suppressed"<<std::endl;
			else
				serr<<" WARNING !! baseEdge = "<<baseEdge<<sendl;

			if (edgeToRemove.size()>0)
			{
				_edgeMod->removeEdges(edgeToRemove,false);

			}
		}
	}

	const core::MechanicalParams* _mparams = core::MechanicalParams::defaultInstance();

	this->apply(_mparams /* PARAMS FIRST */, *this->getToModel()->write(sofa::core::VecCoordId::position()),*this->getFromModel()->read(sofa::core::ConstVecCoordId::position()));

	const Data<InVecCoord>& xfree_in = *this->getFromModel()->read(sofa::core::ConstVecCoordId::freePosition());

	const Data<VecCoord>&     x_out = *this->getToModel()->read(sofa::core::VecCoordId::position());
	const Data<VecCoord>& xfree_out = *this->getToModel()->read(sofa::core::VecCoordId::freePosition());

	core::behavior::MechanicalState<TOut>* ms_out =	dynamic_cast<core::behavior::MechanicalState<TOut> *> (this->getToModel());

	if (x_out.getValue().size() != xfree_out.getValue().size())
	{
		ms_out->vInit(_mparams,sofa::core::VecCoordId::freePosition(),sofa::core::ConstVecCoordId::position());
		ms_out->vInit(_mparams,sofa::core::VecDerivId::freeVelocity(),sofa::core::ConstVecDerivId::velocity());
	}

	if (xfree_in.getValue().size() > 0)
	{
		this->apply(_mparams /* PARAMS FIRST */, *this->getToModel()->write(sofa::core::VecCoordId::freePosition()),	*this->getFromModel()->read(sofa::core::ConstVecCoordId::freePosition()));
	}
}



template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::init()
{
        if (m_ircontroller==NULL) {
                ///////// get the Adaptive Interpolation component ///////
                //std::vector<sofa::core::behavior::LinearSolver*> solvers;
                core::objectmodel::BaseContext * c = this->getContext();

                const helper::vector<std::string>& interpolName = m_controlerPath.getValue();
                if (interpolName.empty()) {
                    m_ircontroller = c->get<TInterventionalRadiologyController>(core::objectmodel::BaseContext::Local);
                } else {
                    m_ircontroller = c->get<TInterventionalRadiologyController>(m_controlerPath.getValue()[0]);
                }

                if(m_ircontroller==NULL)
                    serr<<" no Beam Interpolation found !!! the component can not work"<<sendl;
                else
                    sout<<" interpolation named"<<m_ircontroller->getName()<<" found (for "<<this->getName()<<")"<<sendl;
        }

        m_ircontroller->getInstrumentList(m_instrumentList);

        // create a mapping for each instrument
        m_subMappingList.clear();
        for (unsigned int i=0; i<m_instrumentList.size(); i++)
        {
                AdaptiveBeamMapping< TIn, TOut>* newMapping = new AdaptiveBeamMapping< TIn, TOut>(this->fromModel, this->toModel,m_instrumentList[i],true);
                m_subMappingList.push_back(newMapping);
        }


	this->f_listening.setValue(true);

	///////// STEP 3 : get the edgeSet topology and fill it with segments
	this->getContext()->get(_topology);

	if(_topology != NULL && this->f_printLog.getValue() )
		sout<<" FIND topology named "<< _topology->getName()<<sendl;

	this->getContext()->get(_edgeMod);

	if (_edgeMod == NULL)
		serr << "EdgeSetController has no binding EdgeSetTopologyModifier." << sendl;

	// fill topology :
	_topology->clear();
	_topology->cleanup();

	unsigned int numSeg, numLinesInstrument;
	numSeg=0;
	InReal DX=0;
	// we chose the collision parameters of the most discrestized instrument
	for (unsigned int i=0; i<m_instrumentList.size(); i++)
	{
		InReal dx=0;
		m_instrumentList[i]->getNumberOfCollisionSegment(dx, numLinesInstrument);
		if( numSeg < numLinesInstrument ){
			numSeg = numLinesInstrument;
			DX=dx;
		}
	}

	sout<<" aaaaaaaaaaa numSeg found in MultiAdaptiveBeamMapping="<<numSeg<<sendl;


	// add points
	for ( int i=0; i<(int)numSeg+1; i++)
	{
		Real px = i*DX;
		_topology->addPoint( px, 0, 0);
	}
	// add segments
	for (int i=0; i<(int)numSeg; i++)
	{
		_topology->addEdge(i,i+1);
	}
	// resize Mstate
	this->toModel->resize(numSeg+1);

	// resize the internal list of the collision points ( for each point : [x_curv on the global wire , id of the corresponding instrument]
	_xPointList.resize(numSeg+1);
	_idm_instrumentList.resize(numSeg+1);

}


template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::bwdInit()
{

	for (unsigned int i=0; i<m_instrumentList.size(); i++)
	{
		m_subMappingList[i]->setuseCurvAbs(this->useCurvAbs.getValue());
		m_subMappingList[i]->setName(" SubMapping - " + m_instrumentList[i]->getName() );
		m_subMappingList[i]->bwdInit();//////////////////////////////////////////////////
	}
	assignSubMappingFromControllerInfo();

}



template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::draw(const core::visual::VisualParams* vparams)
{
	if (!vparams->displayFlags().getShowMappings()) return;
}





} // namespace mapping

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_MAPPING_MULTIADAPTIVEBEAMMAPPING_INL */
