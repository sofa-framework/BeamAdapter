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
#pragma once

#include <BeamAdapter/component/mapping/MultiAdaptiveBeamMapping.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/core/MechanicalParams.h>
#include <sofa/helper/ScopedAdvancedTimer.h>


namespace sofa::component::mapping
{

using sofa::helper::ScopedAdvancedTimer;


template <class TIn, class TOut>
MultiAdaptiveBeamMapping< TIn, TOut>::MultiAdaptiveBeamMapping(core::State< In >* from, core::State< Out >* to, TInterventionalRadiologyController* _ircontroller)
: Inherit(from, to)
, useCurvAbs(initData(&useCurvAbs,true,"useCurvAbs","true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same"))
, m_controlerPath(initData(&m_controlerPath,"ircontroller", "Path to the ircontroller component on scene"))
, d_parallelMapping(initData(&d_parallelMapping, false, "parallelMapping", "flag to enable parallel internal computation in all the submappings"))
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
, d_parallelMapping(initData(&d_parallelMapping, false, "parallelMapping", "flag to enable parallel internal computation in all the submappings"))
, m_ircontroller(nullptr)
, isBarycentricMapping(false)
{
    this->addAlias(&m_controlerPath, "controller");
}


template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::apply(const core::MechanicalParams* mparams /* PARAMS FIRST */, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    SCOPED_TIMER("MultiAdaptiveBeamMapping_apply");

    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    if(!isBarycentricMapping)
        out.resize(_xPointList.size());
    else if(out.size() != this->getMechTo()[0]->getSize())
        out.resize(this->getMechTo()[0]->getSize());

    for (unsigned int subMap=0; subMap<m_subMappingList.size(); subMap++)
    {
        if (this->f_printLog.getValue()){
            msg_info()<<"apply ";
            m_subMappingList[subMap]->printIstrumentInfo();//ctn_DEV
        }
        m_subMappingList[subMap]->apply(mparams /* PARAMS FIRST */, dOut, dIn);
    }

}

template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::applyJ(const core::MechanicalParams* mparams /* PARAMS FIRST */, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    SCOPED_TIMER("MultiAdaptiveBeamMapping_applyJ");

    for (unsigned int subMap=0; subMap<m_subMappingList.size(); subMap++)
    {
        if (this->f_printLog.getValue()){
            m_subMappingList[subMap]->printIstrumentInfo();//ctn_DEV
        }
        m_subMappingList[subMap]->applyJ(mparams /* PARAMS FIRST */, dOut, dIn);
    }
}

template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::applyJT(const core::MechanicalParams* mparams /* PARAMS FIRST */, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{
    SCOPED_TIMER("MultiAdaptiveBeamMapping_applyJT");

    for (unsigned int subMap=0; subMap<m_subMappingList.size(); subMap++)
    {
        if (this->f_printLog.getValue()){

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
    SCOPED_TIMER("MultiAdaptiveBeamMapping_applyJT");

    for (unsigned int subMap=0; subMap<m_subMappingList.size(); subMap++)
    {
        if (this->f_printLog.getValue()){
            m_subMappingList[subMap]->printIstrumentInfo();
        }
        m_subMappingList[subMap]->applyJT(cparams /* PARAMS FIRST */, dOut, dIn);
    }
}




template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::handleEvent(sofa::core::objectmodel::Event * event)
{
    if (sofa::simulation::AnimateBeginEvent::checkEventType(event))
    {
        assignSubMappingFromControllerInfo();
    }
}

template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::assignSubMappingFromControllerInfo()
{
    sofa::type::vector<int> removeEdgeAtPoint;

    // 1. get the new controls
    m_ircontroller->interventionalRadiologyCollisionControls(_xPointList, _idm_instrumentList, removeEdgeAtPoint);
    
    if(!isBarycentricMapping)
    {
        
        //Case if this is not a barycentric mapping
        // 2. assign each value of xPointList to the corresponding "sub Mapping"
        // i.e = each point from xPointList of the collision model is controlled by a given instrument wich id is provided in _id_m_instrumentList
        for (unsigned int i=0; i<m_subMappingList.size(); i++)
        {
            auto pointList = sofa::helper::getWriteOnlyAccessor(m_subMappingList[i]->d_points);
            pointList.clear();
            m_subMappingList[i]->clearIdPointSubMap();
        }

        for (unsigned int i=0; i< _xPointList.size(); i++)
        {
            unsigned int  id = _idm_instrumentList[i];

            auto pointList = sofa::helper::getWriteOnlyAccessor(m_subMappingList[id]->d_points);
            pointList.wref().emplace_back(  _xPointList[i], 0, 0 );

            m_subMappingList[id]->addIdPointSubMap(i);
        }

        // handle the possible topological change
        sofa::type::vector<sofa::core::topology::BaseMeshTopology::EdgeID> edgeToRemove;

        for (unsigned int i=0; i<removeEdgeAtPoint.size(); i++)
        {
            /// look if the edge is already suppressed (only one edge around the targetted point)
            sofa::type::vector<sofa::core::topology::BaseMeshTopology::PointID> baseEdge(0);
            baseEdge = _topology->getEdgesAroundVertex((sofa::core::topology::BaseMeshTopology::PointID) removeEdgeAtPoint[i]);

            /// need to be removed
            if (baseEdge.size() == 2)
            {
                /// chose the edge to remove
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
            {
                msg_info()<<" ok, the edge is already suppressed";
            }
            else
            {
                msg_error() << "Trying to remove baseEdge which is alreay empty. This case is not supposed to happened.";
            }

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
    if (m_ircontroller==nullptr) 
    {
        ///////// get the Adaptive Interpolation component ///////
        core::objectmodel::BaseContext * c = this->getContext();

        const type::vector<std::string>& interpolName = m_controlerPath.getValue();
        if (interpolName.empty()) {
            m_ircontroller = c->get<TInterventionalRadiologyController>(core::objectmodel::BaseContext::Local);
        } else {
            m_ircontroller = c->get<TInterventionalRadiologyController>(m_controlerPath.getValue()[0]);
        }

        if (m_ircontroller == nullptr) {
            msg_error() << " no Beam Interpolation found !!! the component can not work";
            sofa::core::objectmodel::BaseObject::d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
            return;
        }
        else {
            msg_info() << " interpolation named" << m_ircontroller->getName() << " found (for " << this->getName() << ")";
        }
    }

    m_ircontroller->getInstrumentList(m_instrumentList);

    // create a mapping for each instrument
    m_subMappingList.clear();
    for (unsigned int i=0; i<m_instrumentList.size(); i++)
    {
        typename AdaptiveBeamMapping< TIn, TOut>::SPtr newMapping = sofa::core::objectmodel::New<AdaptiveBeamMapping< TIn, TOut>>(this->fromModel, this->toModel,m_instrumentList[i],true);
        newMapping->d_parallelMapping.setParent(&d_parallelMapping);
        newMapping->d_parallelMapping.update();
        m_subMappingList.push_back(newMapping);
    }


    this->f_listening.setValue(true);

    ///////// STEP 3 : get the edgeSet topology and fill it with segments
    this->getContext()->get(_topology);

    if (_topology != nullptr && this->f_printLog.getValue()) {
        msg_info() << " FIND topology named " << _topology->getName();
    }

    this->getContext()->get(_edgeMod);

    if (_edgeMod == nullptr)
        msg_error() << "EdgeSetController has no binding EdgeSetTopologyModifier.";

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

    msg_info() << "numSeg found in MultiAdaptiveBeamMapping="<< numSeg;


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

    // create edge around vertex array
    _topology->init();

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
        m_subMappingList[i]->setUseCurvAbs(this->useCurvAbs.getValue());
        m_subMappingList[i]->setName(" SubMapping - " + m_instrumentList[i]->getName() );
        m_subMappingList[i]->bwdInit();//////////////////////////////////////////////////
    }
    assignSubMappingFromControllerInfo();

}


template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::setBarycentricMapping()
{
    isBarycentricMapping=true;
    for(unsigned int i=0;i<m_subMappingList.size();i++)
    {
        m_subMappingList[i]->setBarycentricMapping();
    }
}

template <class TIn, class TOut>
int MultiAdaptiveBeamMapping< TIn, TOut>::addBaryPoint(const int& edgeId,const Vec3& _baryCoord,bool isStraight)
{
    int returnId=this->getMechTo()[0]->getSize();
    this->getMechTo()[0]->resize(returnId+1);

    assert(m_ircontroller !=nullptr && isBarycentricMapping);
    const type::vector<type::vector<int> >& id_instrument_curvAbs_table = m_ircontroller->get_id_instrument_curvAbs_table();
    int nbControlledEdge  = static_cast<int>(id_instrument_curvAbs_table.size()) - 1;
    int totalNbEdges = m_ircontroller->getTotalNbEdges();
    int nbUnControlledEdges = totalNbEdges - nbControlledEdge;
    assert(nbUnControlledEdges>=0);

    if (edgeId < (totalNbEdges-nbControlledEdge) )
    {
        //if the edge in question is not under control, dont need to compute for collision
    }
    else
    {
        int controledEdgeId = edgeId-nbUnControlledEdges;
        const sofa::type::vector<int>&  id_instrument_table_on_node = id_instrument_curvAbs_table[controledEdgeId+1];
        sofa::type::vector< sofa::component::fem::WireBeamInterpolation<In>  *> m_instrumentsList;
        m_ircontroller->getInstrumentList(m_instrumentsList);
        Real radius = m_instrumentsList[id_instrument_table_on_node[0]]->getBeamSection(controledEdgeId)._r;
        int idInstrument  = id_instrument_table_on_node[0];
        for(unsigned int i=1;i<id_instrument_table_on_node.size();i++)
        {
            if(radius < m_instrumentsList[id_instrument_table_on_node[i]]->getBeamSection(controledEdgeId)._r)
            {
                radius = m_instrumentsList[id_instrument_table_on_node[i]]->getBeamSection(controledEdgeId)._r;
                idInstrument = id_instrument_table_on_node[i];
            }
        }
        m_subMappingList[idInstrument]->addBaryPoint(controledEdgeId,_baryCoord,isStraight);
        m_subMappingList[idInstrument]->addIdPointSubMap(returnId);
    }
    return returnId;
}

template <class TIn, class TOut>
void MultiAdaptiveBeamMapping< TIn, TOut>::clear(int size)
{
    for(unsigned int i=0;i<m_subMappingList.size();i++)
        m_subMappingList[i]->clear(size);
    this->getMechTo()[0]->resize(0);
}


} // namespace sofa::component::mapping
