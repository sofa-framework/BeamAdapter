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
// C++ Implementation : AdaptiveBeamMapping
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
#include <BeamAdapter/component/mapping/AdaptiveBeamMapping.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <string>
#include <sofa/core/Mapping.inl>
#include <sofa/core/visual/VisualParams.h>

#include <sofa/simulation/AnimateBeginEvent.h>
#include <sofa/helper/ScopedAdvancedTimer.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/core/MechanicalParams.h>

#include <sofa/simulation/ParallelForEach.h>
#include <sofa/simulation/MainTaskSchedulerFactory.h>

namespace sofa::component::mapping
{

namespace _adaptivebeammapping_
{

using sofa::core::State;
using helper::ReadAccessor;
using helper::WriteAccessor;
using sofa::core::ConstVecCoordId;
using sofa::helper::AdvancedTimer;
using sofa::helper::ScopedAdvancedTimer;
using sofa::core::MultiVecCoordId;
using sofa::core::VecCoordId;
using sofa::core::ConstMultiVecCoordId;
using core::MechanicalParams;

template <class TIn, class TOut>
AdaptiveBeamMapping<TIn,TOut>::AdaptiveBeamMapping(State< In >* from, State< Out >* to,
                                                   BeamInterpolation< TIn >* interpolation, bool isSubMapping)
    : Inherit(from, to)
    , d_useCurvAbs(initData(&d_useCurvAbs,true,"useCurvAbs","true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same"))
    , d_points(initData(&d_points, "points", "defines the mapped points along the beam axis (in beam frame local coordinates)"))
    , d_proximity(initData(&d_proximity, 0.0, "proximity", "if positive, the mapping is modified for the constraints to take into account the lever created by the proximity"))
    , d_contactDuplicate(initData(&d_contactDuplicate,false,"contactDuplicate","if true, this mapping is a copy of an input mapping and is used to gather contact points (ContinuousFrictionContact Response)"))
    , d_inputMapName(initData(&d_inputMapName,"nameOfInputMap", "if contactDuplicate==true, it provides the name of the input mapping"))
    , d_nbPointsPerBeam(initData(&d_nbPointsPerBeam, 0.0, "nbPointsPerBeam", "if non zero, we will adapt the points depending on the discretization, with this num of points per beam (compatible with useCurvAbs)"))
    , d_segmentsCurvAbs(initData(&d_segmentsCurvAbs, "segmentsCurvAbs", "the abscissa of each point on the collision model", true, true))
    , d_parallelMapping(initData(&d_parallelMapping, false, "parallelMapping", "flag to enable parallel internal computation"))
    , l_adaptativebeamInterpolation(initLink("interpolation", "Path to the Interpolation component on scene"), interpolation)
    , m_inputMapping(nullptr)
    , m_isSubMapping(isSubMapping)
    , m_isBarycentricMapping(false)
{
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::init()
{
    if (!l_adaptativebeamInterpolation)
        l_adaptativebeamInterpolation.set(dynamic_cast<sofa::core::objectmodel::BaseContext *>(this->getContext())->get<BInterpolation>());

    if (!l_adaptativebeamInterpolation)
    {
        msg_error() <<"No Beam Interpolation found, the component can not work.";

        this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
    }

    if (d_parallelMapping.getValue())
    {
        simulation::MainTaskSchedulerFactory::createInRegistry()->init();
    }
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::reinit()
{
    init();
    computeDistribution();
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::reset()
{
    reinit();
}

template <class TIn, class TOut>
void AdaptiveBeamMapping<TIn,TOut>::printIstrumentInfo() const
{
    if (m_isSubMapping)
    {
        msg_info()<<"Instrument Named "<<l_adaptativebeamInterpolation->getName()<<msgendl
                 <<" MState1:"<<fromModel->getName()<< "  size:"<<fromModel->getSize()<<msgendl
                <<" MState2:"<<toModel->getName()<< "  size:"<<toModel->getSize()<<msgendl
               <<"idPointSubMap."<<m_idPointSubMap.size()<<msgendl
              <<"pointBeamDistribution."<<m_pointBeamDistribution.size();
    }
}


template <class TIn, class TOut>
int AdaptiveBeamMapping< TIn, TOut>::addPoint (const Coord& point, int indexFrom)
{
    SOFA_UNUSED(indexFrom) ;

    auto points = sofa::helper::getWriteAccessor(d_points);

    const auto nbPoints = points.size();
    const auto coord = TOut::getCPos(point);

    points.wref().emplace_back(coord[0], coord[1], coord[2]);

    return static_cast<int>(nbPoints);
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::setBarycentricMapping()
{
    m_isBarycentricMapping=true;

    auto points = sofa::helper::getWriteOnlyAccessor(d_points);
    points.clear();
}


template <class TIn, class TOut>
int AdaptiveBeamMapping< TIn, TOut>::addBaryPoint(const int& beamId, const Vec3& baryCoord, bool straightlineSplineOption)
{
    //TODO(dmarchal 2017-06-01) Please specify who/when this will be done (remove in 1 year)
    //TODO add parameter label for different cases : unactive, linear, spline
    //attention, beamId here is not the edge Id, but the id of a vec_edge_list defined in BeamInterpolation
    SOFA_UNUSED(straightlineSplineOption);

    const auto newPointId = m_pointBeamDistribution.size();
    m_pointBeamDistribution.resize(newPointId+1);
    m_pointBeamDistribution[newPointId].baryPoint=baryCoord;
    m_pointBeamDistribution[newPointId].beamId=beamId;
    return static_cast<int>(newPointId);
}


//void clear(){}; /// CTN_DEV todo for ContactMapper


//clear the mapping in functions of size given
template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::clear(int size)
{
    clearIdPointSubMap();
    m_pointBeamDistribution.clear();
    auto points = sofa::helper::getWriteOnlyAccessor(d_points);
    if(size>0 && !m_isSubMapping)
    {
        m_pointBeamDistribution.reserve(size);
        points->reserve(size);
        this->getMechTo()[0]->resize(size);
    }
    else //case where this clear is call by a Multimapping, all component will be clear to null size
    {
        points->resize(0);
        this->getMechTo()[0]->resize(0);
    }
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::apply(const MechanicalParams* mparams, Data<VecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    SCOPED_TIMER("AdaptiveBeamMapping_Apply");

    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    const InVecCoord& in = dIn.getValue();

    m_isXBufferUsed=false;

    // When using an adaptatif controller, one need to redistribute the points at each time step
    {
        SCOPED_TIMER("computeDistribution");
        if (d_useCurvAbs.getValue() && !d_contactDuplicate.getValue())
            computeDistribution();
    }

    if (!m_isSubMapping)
    {
        SCOPED_TIMER("resizeToModel&Out");
        if (d_nbPointsPerBeam.getValue() > 0)
        {
            this->toModel->resize(m_pointBeamDistribution.size());
            out.resize(m_pointBeamDistribution.size());
        }
        else
        {
            this->toModel->resize(d_points.getValue().size());
            out.resize(d_points.getValue().size());
        }
    }

    MultiVecCoordId x = VecCoordId::position();
    MultiVecCoordId xfree = VecCoordId::freePosition();

    const ConstMultiVecCoordId &xId = mparams->x();
    ConstVecCoordId xtest = xId.getId(this->fromModel);

    if(xtest == xfree.getId(this->fromModel))
    {
        VecCoordId xfreeIn = VecCoordId::freePosition();
        l_adaptativebeamInterpolation->updateBezierPoints(in, xfreeIn);
    }
    else if(xtest == x.getId(this->fromModel))
    {
        VecCoordId positionIn = VecCoordId::position();
        l_adaptativebeamInterpolation->updateBezierPoints(in, positionIn);
    }

    AdvancedTimer::stepBegin("computeNewInterpolation");

    // effective computation for each element in pointBeamDistribution
    auto apply_impl = [&](const auto& range)
    {
        auto elementID = std::distance(m_pointBeamDistribution.begin(), range.start);
        for (auto it = range.start; it != range.end; ++it, ++elementID)
        {
            const auto& pointBeamDistribution = m_pointBeamDistribution[elementID];

            Vec<3, InReal> pos(sofa::type::NOINIT);
            const Vec3 localPos(0., pointBeamDistribution.baryPoint[1], pointBeamDistribution.baryPoint[2]);
            l_adaptativebeamInterpolation->interpolatePointUsingSpline(pointBeamDistribution.beamId, pointBeamDistribution.baryPoint[0], localPos, in, pos, false, xtest);

            if (m_isSubMapping)
            {
                if (!m_idPointSubMap.empty())
                    out[m_idPointSubMap[elementID]] = pos;
            }
            else
                out[elementID] = pos;
        }

    };

    const bool multithreading = d_parallelMapping.getValue();

    const simulation::ForEachExecutionPolicy execution = multithreading ?
        simulation::ForEachExecutionPolicy::PARALLEL :
        simulation::ForEachExecutionPolicy::SEQUENTIAL;

    simulation::TaskScheduler* taskScheduler = simulation::MainTaskSchedulerFactory::createInRegistry();
    assert(taskScheduler);

    simulation::forEachRange(execution, *taskScheduler, m_pointBeamDistribution.begin(), m_pointBeamDistribution.end(), apply_impl);


    AdvancedTimer::stepEnd("computeNewInterpolation");
}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJ(const core::MechanicalParams* mparams, Data<VecDeriv>& dOut, const Data<InVecDeriv>& dIn)
{
    SOFA_UNUSED(mparams);
    SCOPED_TIMER("AdaptiveBeamMapping_applyJ");

    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    const InVecDeriv& in= dIn.getValue();
    
    Data<InVecCoord>& dataInX = *this->getFromModel()->write(VecCoordId::position());
    auto x = sofa::helper::getWriteOnlyAccessor(dataInX);

    if (d_useCurvAbs.getValue() && !d_contactDuplicate.getValue())
        computeDistribution();

    // TODO: check if m_isXBufferUsed is set somewhere else
    // As far as I could see, m_isXBufferUsed is never set to true
    InVecCoord xBuf2{};
    if (m_isXBufferUsed)
    {
        // TODO : solve this problem during constraint motion propagation !!
        xBuf2 = x.ref();
        x.wref() = m_xBuffer;
    }

    // should not be necessary if apply() was called first
    if (!m_isSubMapping)
    {
        if (d_nbPointsPerBeam.getValue() > 0)
        {
            out.resize(m_pointBeamDistribution.size());
        }
        else
        {
            out.resize(d_points.getValue().size());
        }
    }

    // effective computation for each element in pointBeamDistribution
    auto applyJ_impl = [&](const auto& range)
    {
        auto elementID = std::distance(m_pointBeamDistribution.begin(), range.start);
        for (auto it = range.start; it != range.end; ++it, ++elementID)
        {
            const auto& pointBeamDistribution = m_pointBeamDistribution[elementID];

            unsigned int IdxNode0, IdxNode1;
            l_adaptativebeamInterpolation->getNodeIndices(pointBeamDistribution.beamId, IdxNode0, IdxNode1);

            SpatialVector vDOF0{ In::getDRot(in[IdxNode0]), In::getDPos(in[IdxNode0]) };
            SpatialVector vDOF1{ In::getDRot(in[IdxNode1]), In::getDPos(in[IdxNode1]) };

            Deriv vResult(sofa::type::NOINIT);

            applyJonPoint(elementID, vDOF0, vDOF1, vResult, x.ref());

            if (m_isSubMapping)
            {
                if (m_idPointSubMap.size() > 0)
                    out[m_idPointSubMap[elementID]] = vResult;
            }
            else
                out[elementID] = vResult;
        }
    };

    const bool multithreading = d_parallelMapping.getValue();

    const simulation::ForEachExecutionPolicy execution = multithreading ?
        simulation::ForEachExecutionPolicy::PARALLEL :
        simulation::ForEachExecutionPolicy::SEQUENTIAL;

    simulation::TaskScheduler* taskScheduler = simulation::MainTaskSchedulerFactory::createInRegistry();
    assert(taskScheduler);

    simulation::forEachRange(execution, *taskScheduler, m_pointBeamDistribution.begin(), m_pointBeamDistribution.end(), applyJ_impl);

    if(m_isXBufferUsed)
    {
        x.wref() = xBuf2;
        m_isXBufferUsed = false;
    }
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJT(const core::MechanicalParams* mparams, Data<InVecDeriv>& dOut, const Data<VecDeriv>& dIn)
{
    SOFA_UNUSED(mparams);

    SCOPED_TIMER("AdaptiveBeamMapping_applyJT");
    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    const VecDeriv& in= dIn.getValue();

    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(ConstVecCoordId::position());
    const InVecCoord& x = dataInX.getValue();

    for (unsigned int i=0; i<m_pointBeamDistribution.size(); i++)
    {
        PosPointDefinition  ppd = m_pointBeamDistribution[i];
        //1. get the indices
        unsigned int IdxNode0, IdxNode1;
        l_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

        Deriv finput;
        if(m_isSubMapping){
            if(m_idPointSubMap.size()>0)
                finput = in[m_idPointSubMap[i]];}
        else
            finput = in[i];

        SpatialVector FNode0, FNode1;
        applyJTonPoint(i, finput, FNode0, FNode1, x);

        //2. put the result in out vector computes the equivalent forces on nodes + rotate to Global Frame from DOF frame
        In::setDPos(out[IdxNode0], In::getDPos(out[IdxNode0]) + FNode0.getForce());
        In::setDPos(out[IdxNode1], In::getDPos(out[IdxNode1]) + FNode1.getForce());
        In::setDRot(out[IdxNode0], In::getDRot(out[IdxNode0]) + FNode0.getTorque());
        In::setDRot(out[IdxNode1], In::getDRot(out[IdxNode1]) + FNode1.getTorque());
    }
}


/// AdaptiveBeamMapping::applyJT(InMatrixDeriv& out, const OutMatrixDeriv& in)
/// this function propagate the constraint through the Adaptive Beam mapping :
/// if one constraint along (vector n) with a value (v) is applied on the childModel (like collision model)
/// then this constraint is transformed by (Jt.n) with value (v) for the rigid model
/// note : the value v is not propagated through the mapping
template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJT(const core::ConstraintParams* cparams, Data<InMatrixDeriv>& dOut, const Data<OutMatrixDeriv>& dIn)
{
    SOFA_UNUSED(cparams);

    SCOPED_TIMER("AdaptiveBeamMapping_ApplyJT");
    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    const OutMatrixDeriv& in = dIn.getValue();
    const Data<InVecCoord>& dataInX = *this->getFromModel()->read(ConstVecCoordId::position());
    const InVecCoord& x = dataInX.getValue();

    m_isXBufferUsed = false;
    m_xBuffer = x ;

    //////////// What's for ?? it seems not useful//////////
    bool proximity_lever = false;
    if (d_proximity.getValue() > 0.0)
        proximity_lever = true;

    if (proximity_lever && this->f_printLog.getValue() )
        msg_warning() <<"The constraints are contact at the surface of the beam (not at their center)";
    /////////////////////////////////////

    typename Out::MatrixDeriv::RowConstIterator rowItEnd = in.end();
    for (typename Out::MatrixDeriv::RowConstIterator rowIt = in.begin(); rowIt != rowItEnd; ++rowIt)
    {
        typename Out::MatrixDeriv::ColConstIterator colItEnd = rowIt.end();
        for (typename Out::MatrixDeriv::ColConstIterator colIt = rowIt.begin(); colIt != colItEnd; ++colIt)
        {
            typename In::MatrixDeriv::RowIterator o = out.wref().writeLine(rowIt.index());
            unsigned int indexIn = colIt.index();
            const Deriv data = colIt.val();

            if (m_isSubMapping)
            {
                // look if we get the indexIn in the idPointSubMap:
                unsigned int i=0;
                while( i<m_idPointSubMap.size() && m_idPointSubMap[i]!=indexIn)
                    i++;

                if (i<m_idPointSubMap.size())
                    indexIn = i;
                else
                    continue;
            }

            if(indexIn<m_pointBeamDistribution.size()){
                const PosPointDefinition&  ppd = m_pointBeamDistribution[indexIn];
                unsigned int IdxNode0, IdxNode1;
                l_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

                SpatialVector FNode0, FNode1;
                applyJTonPoint(indexIn, data, FNode0, FNode1, x);

                // Compute the mapped Constraint on the beam nodes
                InDeriv direction0;
                In::setDPos(direction0,FNode0.getForce());
                In::setDRot(direction0,FNode0.getTorque());
                InDeriv direction1;
                In::setDPos(direction1,FNode1.getForce());
                In::setDRot(direction1,FNode1.getTorque());

                o.addCol(IdxNode0, direction0);
                o.addCol(IdxNode1, direction1);
            }
            else
            {
                if ( this->f_printLog.getValue() )
                    msg_warning() <<"Wrong index in VecConst in";
                break;
            }
        }
    }
}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::bwdInit()
{
    const auto& pts = d_points.getValue();
    const auto ptsSize = pts.size();

    if (ptsSize == 0)
    {
        helper::ReadAccessor<Data<VecCoord> > xTo = this->toModel->read(sofa::core::ConstVecCoordId::position()) ;

        if(xTo.size()==0)
        {
            msg_warning() <<"No point defined in the AdaptiveBeamMapping ";
            // _problem = true;
        }
        else
        {
            auto pts2 = sofa::helper::getWriteOnlyAccessor(d_points);
            msg_info() <<"No point defined in the AdaptiveBeamMapping - uses positions defined by Mechanical State";
            for(unsigned int i=0; i<xTo.size();i++)
            {
                pts2.wref().emplace_back(xTo[i][0], xTo[i][1], xTo[i][2]);
            }
        }
    }

    bool curvAbs = this->d_useCurvAbs.getValue();
    if (curvAbs)
    {
        int cpt=0;
        std::stringstream tmp;
        for (unsigned int i=0; cpt < 10 && i<pts.size()-1; i++)
        {
            if( pts[i][0]>pts[i+1][0])
            {
                tmp << "- invalid ordering at index: " << i << " : " << pts[i][0] << " > " << pts[i][1] << msgendl ;
                cpt++;
            }
        }
        if(cpt>=10)
        {
            msg_warning() <<"When using useCurvAbs==true, points must be sorted according to their curvAbs: " << msgendl
                          << tmp.str() ;
        }
    }

    l_adaptativebeamInterpolation->bwdInit();
    computeDistribution();
    if (!m_isSubMapping)
    {
        core::Mapping< TIn, TOut>::init();
    }
    if(d_contactDuplicate.getValue()==true)
    {
        this->fromModel->getContext()->get(m_inputMapping, sofa::core::objectmodel::BaseContext::SearchRoot);
        if(m_inputMapping==nullptr)
            msg_warning() <<"Can not found the input  Mapping";
        else
            msg_info()<<"Input Mapping named "<<m_inputMapping->getName()<<" is found";
    }
}

template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::beginAddContactPoint()
{
    auto points = sofa::helper::getWriteOnlyAccessor(d_points);
    points->clear();
    m_pointBeamDistribution.clear();
}


template <class TIn, class TOut>
int AdaptiveBeamMapping< TIn, TOut>::addContactPoint(const Vec3& bary)
{
    unsigned int index = d_points.getValue().size();
    auto points = sofa::helper::getWriteOnlyAccessor(d_points);
    points->push_back(bary);

    if(this->toModel->getSize() <= index)
        this->toModel->resize(index+1);
    return index;
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::computeIdxAndBaryCoordsForAbs(unsigned int &b, Real &xBary, const Real &xAbs )
{
    const InReal xAbsInput = (InReal) xAbs;
    InReal xBaryOutput = (InReal) xBary;

    l_adaptativebeamInterpolation->getBeamAtCurvAbs(xAbsInput,b,xBaryOutput);
    xBary = (Real) xBaryOutput;
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::computeDistribution()
{
    //The normal procedure is givng "points", computeDistribution comppute theses points and put
    //on the baryPoints which are "pointBeamDistribution"
    //If the mapping is a barycentric one, that mean a list of baryPoints are already setted
    //Then no need to recompute this
    if(!m_isBarycentricMapping)
    {
        bool curvAbs = d_useCurvAbs.getValue();
        const type::vector<Vec3>& points = d_points.getValue();
        m_pointBeamDistribution.clear();

        const unsigned int numBeams = l_adaptativebeamInterpolation->getNumBeams();
        if(numBeams==0)
        {
            if (this->f_printLog.getValue())
                msg_info() <<"No beams found in adaptBeamInterpolation in BeamInterpolation named"<<l_adaptativebeamInterpolation->getName();
            return;
        }

        if (curvAbs)
        {
            double ptsPerBeam = d_nbPointsPerBeam.getValue();
            if(ptsPerBeam)
            {	// Recreating the distribution based on the current sampling of the beams
                WriteAccessor<Data<type::vector<Real>>> waSegmentsCurvAbs = d_segmentsCurvAbs;
                waSegmentsCurvAbs.clear();

                double step = 1.0 / ptsPerBeam;
                double posInBeam = 0;
                InReal segStart, segEnd, segLength;
                for(unsigned int b=0; b< numBeams; ++b)
                {
                    l_adaptativebeamInterpolation->getAbsCurvXFromBeam(b, segStart, segEnd);
                    segLength = segEnd - segStart;

                    for (; posInBeam <= 1.0; posInBeam += step)
                    {
                        waSegmentsCurvAbs.push_back(segStart + segLength * posInBeam);
                        m_pointBeamDistribution.emplace_back( b, Vec3{ posInBeam, 0.0, 0.0} );
                    }

                    posInBeam -= 1.0;
                }

                if (numBeams && fabs(posInBeam - step) > 0.01)
                {
                    // Last point
                    waSegmentsCurvAbs.push_back(segEnd);
                    m_pointBeamDistribution.emplace_back( numBeams - 1, Vec3{ 1.0, 0.0, 0.0} );
                }

                BaseMeshTopology* topo = this->getContext()->getMeshTopology();
                if(topo)
                {
                    topo->clear();
                    int nbEdges = static_cast<int>(m_pointBeamDistribution.size()) - 1;
                    for(int i=0; i<nbEdges; ++i)
                        topo->addEdge(i, i+1);
                }
            }
            else
            {
                // We use the points Data
                m_pointBeamDistribution.reserve(points.size());
                for (unsigned int i=0; i<points.size(); i++)
                {
                    unsigned int b=0;
                    Real xAbs = points[i][0];
                    Real xBary = 0.0;

                    computeIdxAndBaryCoordsForAbs(b, xBary,  xAbs );

                    m_pointBeamDistribution.emplace_back( b, Vec3{ xBary, points[i][1], points[i][2]});
                }
            }
        }
        else
        {
            m_pointBeamDistribution.reserve(points.size());
            for (unsigned int i=0; i<points.size(); i++)
            {
                const unsigned int beamId = (int)floor(points[i][0]);
                msg_warning_when(beamId > numBeams-1) << "Points[" << i << "][0] = " << beamId << " is defined outside of the beam length";
                m_pointBeamDistribution.emplace_back( beamId, Vec3{ points[i][0] - floor(points[i][0]), points[i][1], points[i][2]});
            }
        }
    }
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::computeJacobianOnPoint(unsigned int i, const typename In::VecCoord& x)
{
    // TEST : calcul d'une jacobienne:
    Mat3x12 J;
    Mat12x3 Jt;

    for (unsigned int j=0; j<3; j++)
    {
        Deriv Id, Vresult;
        Vec3 Idv(0,0,0);
        Id[j]=1.0;
        SpatialVector v_DOF0, v_DOF1;

        //3 colonnes
        v_DOF0.clear();
        v_DOF0.setLinearVelocity(Idv);
        v_DOF1.clear();
        applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
        J(0,j)=Vresult[0]; J(1,j)=Vresult[1]; J(2,j)=Vresult[2];
        //3 colonnes
        v_DOF0.clear();
        v_DOF0.setAngularVelocity(Idv);
        v_DOF1.clear();
        applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
        J(0,j+3)=Vresult[0]; J(1,j+3)=Vresult[1]; J(2,j+3)=Vresult[2];
        //3 colonnes
        v_DOF0.clear();
        v_DOF1.clear();
        v_DOF1.setLinearVelocity(Idv);

        applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
        J(0,j+6)=Vresult[0]; J(1,j+6)=Vresult[1]; J(2,j+6)=Vresult[2];
        //3 colonnes
        v_DOF0.clear();
        v_DOF1.clear();
        v_DOF1.setAngularVelocity(Idv);
        applyJonPoint(i, v_DOF0, v_DOF1, Vresult, x);
        J(0,j+9)=Vresult[0]; J(1,j+9)=Vresult[1]; J(2,j+9)=Vresult[2];

        SpatialVector F_DOF0, F_DOF1;
        applyJTonPoint(i, Id, F_DOF0, F_DOF1, x);
        Jt(0,j)=F_DOF0.getForce()[0]; Jt(1,j)=F_DOF0.getForce()[1];  Jt(2,j) =F_DOF0.getForce()[2];
        Jt(3,j)=F_DOF0.getTorque()[0];Jt(4,j)=F_DOF0.getTorque()[1]; Jt(5,j) =F_DOF0.getTorque()[2];
        Jt(6,j)=F_DOF1.getForce()[0]; Jt(7,j)=F_DOF1.getForce()[1];  Jt(8,j) =F_DOF1.getForce()[2];
        Jt(9,j)=F_DOF1.getTorque()[0];Jt(10,j)=F_DOF1.getTorque()[1];Jt(11,j)=F_DOF1.getTorque()[2];
    }
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJonPoint(unsigned int i, SpatialVector& VNode0input, SpatialVector& VNode1input, Deriv& vOutput, const InVecCoord& x)
{

    // 1. get the curvilinear abs;
    const PosPointDefinition&  ppd = m_pointBeamDistribution[i];

    // 2. get the indices
    unsigned int IdxNode0, IdxNode1;
    l_adaptativebeamInterpolation->getNodeIndices(ppd.beamId,IdxNode0,IdxNode1);

    // 3. get the transform to DOF in global frame from local frame
    Transform DOF0Global_H_local0, DOF1Global_H_local1;
    l_adaptativebeamInterpolation->getDOFtoLocalTransformInGlobalFrame(ppd.beamId, DOF0Global_H_local0, DOF1Global_H_local1, x);

    // 4. project the velocities in local frame:
    SpatialVector v_local0, v_local1;
    v_local0 = DOF0Global_H_local0.inversed()*VNode0input;
    v_local1 = DOF1Global_H_local1.inversed()*VNode1input;

    // 5. Computes the local velocities of the 4 points of the spline
    Real L = l_adaptativebeamInterpolation->getLength(ppd.beamId);
    Vec3 lever(L/3,0,0);
    Vec3 V0, V1, V2, V3;
    V0 = v_local0.getLinearVelocity();
    V1 = V0 - lever.cross(v_local0.getAngularVelocity());
    V3 = v_local1.getLinearVelocity();
    V2 = V3 + lever.cross(v_local1.getAngularVelocity());

    // 6. Rotate back the vectors in the global frame
    V0 = DOF0Global_H_local0.getOrientation().rotate(V0);
    V1 = DOF0Global_H_local0.getOrientation().rotate(V1);
    V2 = DOF1Global_H_local1.getOrientation().rotate(V2);
    V3 = DOF1Global_H_local1.getOrientation().rotate(V3);

    const Vec3 localPos(0.,ppd.baryPoint[1],ppd.baryPoint[2]);
    if(localPos.norm() > L*1e-10)
    {
        lever = localPos;

        // 7. compute the effect of the angular velocities on the points that are not aligned along the center line
        Vec3 DV0, DV3;
        DV0 = DOF0Global_H_local0.getOrientation().rotate( - lever.cross( v_local0.getAngularVelocity() ) );
        DV3 = DOF1Global_H_local1.getOrientation().rotate( - lever.cross( v_local1.getAngularVelocity() ) );

        // uses spline to interpolate:
        Real bx = ppd.baryPoint[0];
        Real a0=(1-bx)*(1-bx)*(1-bx);
        Real a1=3*bx*(1-bx)*(1-bx);
        Real a2=3*bx*bx*(1-bx);
        Real a3=bx*bx*bx;
        vOutput = V0*a0 + V1*a1 + V2*a2 + V3*a3 + DV0*(a0+a1) + DV3*(a2+a3);

    }
    else
    {
        // uses spline to interpolate:
        Real bx = ppd.baryPoint[0];
        vOutput = V0*(1-bx)*(1-bx)*(1-bx) + V1*3*bx*(1-bx)*(1-bx) + V2*3*bx*bx*(1-bx) + V3*bx*bx*bx;
    }
}


template <class TIn, class TOut>
void AdaptiveBeamMapping< TIn, TOut>::applyJTonPoint(unsigned int i, const Deriv& finput, SpatialVector& FNode0output, SpatialVector& FNode1output, const InVecCoord& x)
{
    const PosPointDefinition&  pointBeamDistribution = m_pointBeamDistribution[i];
    const Vec3 localPos(0.,pointBeamDistribution.baryPoint[1],pointBeamDistribution.baryPoint[2]);
    const Vec3 Fin(finput[0], finput[1], finput[2]);
    l_adaptativebeamInterpolation->MapForceOnNodeUsingSpline(pointBeamDistribution.beamId, pointBeamDistribution.baryPoint[0], localPos, x, Fin, FNode0output, FNode1output );
}

} /// namespace _adaptivebeammapping_

} /// namespace sofa::component::mapping
