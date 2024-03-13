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
// Author: Eulalie Coevoet
//
// Copyright: See COPYING file that comes with this distribution

#pragma once

#include <BeamAdapter/component/mapping/BeamProjectionDifferenceMultiMapping.h>
#include <sofa/core/Multi2Mapping.inl>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/objectmodel/BaseContext.h>
#include <sofa/helper/logging/Message.h>
#include <sofa/type/RGBAColor.h>

#include <Eigen/Dense>

#include <string>

namespace beamadapter::mapping
{
using sofa::core::objectmodel::BaseContext ;
using sofa::helper::WriteAccessor;
using sofa::type::RGBAColor ;
using sofa::core::objectmodel::ComponentState;

template <class TIn1, class TIn2, class TOut>
BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::BeamProjectionDifferenceMultiMapping()
    : d_indices(initData(&d_indices, "indicesInput1", "Indices of model1 to project on model2 (beams)"))
    , d_directions(initData(&d_directions, "directions", "Directions to project (in the local frame)."))
    , d_updateProjectionPosition(initData(&d_updateProjectionPosition, false, "updateProjectionPosition", "Update the projection on the beam at each time step even when direction[0]=1."))
    , d_updateProjectionOrientation(initData(&d_updateProjectionOrientation, false, "updateProjectionOrientation", "Update the projection on the beam at each time step even when direction[0]=1."))
    , d_draw(initData(&d_draw, "draw", "Draw projection points and directions"))
    , d_drawSize(initData(&d_drawSize, Real(3), "drawSize", ""))
    , l_in2Topology(initLink("topologyInput2", "link to input2's topology container (beams to project on)"))
    , l_interpolation(initLink("interpolationInput2", "link to input2's interpolation component (BeamInterpolation)"))
    , m_fromModel1(NULL)
    , m_fromModel2(NULL)
    , m_toModel(NULL)
    , m_updateJ(false)
{
    auto directions = sofa::helper::getWriteAccessor(d_directions);
    directions.resize(OutDeriv::total_size);
}


template <class TIn1, class TIn2, class TOut>
void BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::init()
{
    d_componentState.setValue(ComponentState::Valid);

    m_eigenJacobians.resize( 2 );
    m_eigenJacobians[0] = &m_eigenJacobian1;
    m_eigenJacobians[1] = &m_eigenJacobian2;

    if(this->getFromModels1().empty())
    {
        msg_error() << "Error while initializing, input1 not found." ;
        d_componentState.setValue(ComponentState::Invalid);
        return;
    }
    m_fromModel1 = this->getFromModels1()[0];

    if(this->getFromModels2().empty())
    {
        msg_error() << "Error while initializing, input2 not found." ;
        d_componentState.setValue(ComponentState::Invalid);
        return;
    }
    m_fromModel2 = this->getFromModels2()[0];

    if(this->getToModels().empty())
    {
        msg_error() << "Error while initializing, output not found." ;
        d_componentState.setValue(ComponentState::Invalid);
        return;
    }
    m_toModel = this->getToModels()[0];

    if (l_in2Topology.empty())
    {
        msg_info() << "Link to input2's topology container should be set to ensure right behavior. First Topology found in input2's context will be used.";
        l_in2Topology.set(m_fromModel2->getContext()->getMeshTopologyLink());
        if (l_in2Topology.empty()){
            msg_error() << "Input2's topology container with edges not found.";
            d_componentState.setValue(ComponentState::Invalid);
            return;
        }
    }

    if (l_interpolation.empty())
    {
        msg_error() << "Link to input2's interpolation component is empty. The component cannot work.";
        d_componentState.setValue(ComponentState::Invalid);
        return;
    }

    auto directions = sofa::helper::getWriteAccessor(d_directions);
    if (directions.size() != OutCoord::total_size)
    {
        msg_warning() << "Wrong size for directions, [" << directions << "]. The size should be " << OutCoord::total_size << ".";
        directions.resize(OutCoord::total_size);
    }
}


template <class TIn1, class TIn2, class TOut>
void BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::computeProjection(const In1VecCoord &xFrom, const In2VecCoord &xTo,
                                                                               const bool &updateOrientation)
{
    auto interpolation = l_interpolation.get();
    auto edges = l_in2Topology.get()->getEdges();

    // for each point of P of xFrom
    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    m_mappedPoints.resize(indices.size());

    for (unsigned int i=0; i<indices.size(); i++)
    {
        In1Coord P = xFrom[indices[i]];
        MappedPoint mp;
        bool found = false;
        Real distMin = std::numeric_limits<Real>::max();

        // find the min distance between P and its projection on each edge of xTo
        for (sofa::Size e=0; e<edges.size(); e++)
        {
            const auto& edge = edges[e];
            In2Coord Q1 = xTo[edge[0]];
            In2Coord Q2 = xTo[edge[1]];

            OutDeriv dirAxe = Out::coordDifference(Q2, Q1);
            Real normAxe = In1::getDPos(dirAxe).norm();

            if (std::abs(normAxe)<std::numeric_limits<SReal>::epsilon())  // edge is degenerated, continue with next edge
                continue;

            dirAxe /= normAxe;

            Real r =  In1::getDPos(In1::coordDifference(P, Q1)) * In1::getDPos(dirAxe);
            Real alpha = r / normAxe;
            alpha = (std::abs(alpha)<1e-2)? std::abs(alpha): alpha;

            if (alpha < 0 || alpha > 1)  // not on the edge, continue with next edge
                continue;

            In1Coord proj;
            for (size_t j=0; j<In1Coord::total_size; j++)
                proj[j]= (1 - alpha) * Q1[j] + alpha * Q2[j];
            OutDeriv dirProj = In1::coordDifference(P, proj);
            Real normProj = In1::getDPos(dirProj).norm();

            if (normProj <= distMin)  // closest projection
            {
                distMin = normProj;

                mp.pi1 = edge[0];
                mp.pi2 = edge[1];
                mp.alpha = alpha;
                mp.edgeIndex = e;

                Transform interpolatedTransform;
                interpolation->InterpolateTransformUsingSpline(interpolatedTransform,
                                                               alpha,
                                                               Transform(Q1.getCenter(), Q1.getOrientation()),
                                                               Transform(Q2.getCenter(), Q2.getOrientation()),
                                                               interpolation->getLength(e));
                mp.interpolatedTransform = interpolatedTransform;
                found = true;
            }
        }

        mp.onBeam = found;
        if (!updateOrientation)
        {
            mp.interpolatedTransform = m_mappedPoints[i].interpolatedTransform;
        }
        m_mappedPoints[i] = mp;
    }
}


template <class TIn1, class TIn2, class TOut>
void BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::apply( const sofa::core::MechanicalParams* mparams,
                                                                  const sofa::type::vector<OutDataVecCoord*>& dataVecOutPos,
                                                                  const sofa::type::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
                                                                  const sofa::type::vector<const In2DataVecCoord*>& dataVecIn2Pos)
{
    SOFA_UNUSED(mparams);

    if(d_componentState.getValue() == ComponentState::Invalid)
        return;

    if(dataVecOutPos.empty() || dataVecIn1Pos.empty() || dataVecIn2Pos.empty())
        return;

    const In1VecCoord& in1 = dataVecIn1Pos[0]->getValue();
    const In2VecCoord& in2 = dataVecIn2Pos[0]->getValue();
    OutVecCoord& out = *dataVecOutPos[0]->beginEdit();
    const auto& directions = sofa::helper::getReadAccessor(d_directions);
    auto interpolation = l_interpolation.get();

    if (m_mappedPoints.empty() || !directions[0])
    {
        computeProjection(in1, in2, true);
    }
    else if (d_updateProjectionPosition.getValue() || d_updateProjectionOrientation.getValue())
    {
        computeProjection(in1, in2, d_updateProjectionOrientation.getValue());
    }

    m_updateJ = true;

    sofa::Size sz = m_mappedPoints.size();
    out.resize(sz);
    const auto& indices = sofa::helper::getReadAccessor(d_indices);

    for(sofa::Size i=0; i<sz; i++)
    {
        MappedPoint& mp = m_mappedPoints[i];
        if (mp.onBeam)
        {
            Transform interpolatedTransform;
            interpolation->InterpolateTransformUsingSpline(interpolatedTransform,
                                                           mp.alpha,
                                                           Transform(in2[mp.pi1].getCenter(), in2[mp.pi1].getOrientation()),
                                                           Transform(in2[mp.pi2].getCenter(), in2[mp.pi2].getOrientation()),
                                                           interpolation->getLength(mp.edgeIndex));

            In1Coord projection;
            projection.getCenter() = interpolatedTransform.getOrigin();
            projection.getOrientation() = interpolatedTransform.getOrientation();

            In1Coord p;
            for (size_t j=0; j<In1Coord::total_size; j++)
            {
                p[j] = in1[indices[i]][j] - projection[j];
            }

            In1Coord v;
            v.getCenter() = mp.interpolatedTransform.getRotationMatrix().transposed() * Out::getCPos(p);
            v.getOrientation() = p.getOrientation() * Out::getCRot(p);

            for (size_t j=0; j<In1Coord::total_size; j++)
            {
                out[i][j] = (directions[j])? v[j] : 0;
            }
        }
        else
        {
            out[i] = OutCoord();
        }
    }
    dataVecOutPos[0]->endEdit();
}



template <class TIn1, class TIn2, class TOut>
void BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>:: applyJ(const sofa::core::MechanicalParams* mparams,
                                                                    const sofa::type::vector< OutDataVecDeriv*>& dataVecOutVel,
                                                                    const sofa::type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
                                                                    const sofa::type::vector<const In2DataVecDeriv*>& dataVecIn2Vel)
{
    SOFA_UNUSED(mparams);

    if(d_componentState.getValue() == ComponentState::Invalid)
        return;

    if(dataVecOutVel.empty() || dataVecIn1Vel.empty() || dataVecIn2Vel.empty() )
        return;

    const In1VecDeriv& in1 = dataVecIn1Vel[0]->getValue();
    const In2VecDeriv& in2 = dataVecIn2Vel[0]->getValue();
    OutVecDeriv& outVel = *dataVecOutVel[0]->beginEdit();

    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    const auto& directions = sofa::helper::getReadAccessor(d_directions);

    size_t sz = m_mappedPoints.size();
    outVel.resize(sz);
    for (size_t i = 0 ; i < sz; i++)
    {
        MappedPoint& mp = m_mappedPoints[i];
        if (mp.onBeam)
        {
            In1Deriv vel;
            for (size_t j=0; j<In1Deriv::total_size; j++)
            {
                vel[j] = in1[indices[i]][j] - (1 - mp.alpha) * in2[mp.pi1][j] - mp.alpha * in2[mp.pi2][j];
            }

            In1Deriv v;
            v.getVCenter() = mp.interpolatedTransform.getRotationMatrix().transposed() * Out::getDPos(vel);
            v.getVOrientation() = mp.interpolatedTransform.getRotationMatrix().transposed() * Out::getDRot(vel);

            for (size_t j=0; j<In1Deriv::total_size; j++)
            {
                outVel[i][j] = (directions[j])? v[j] : 0;
            }
        }
    }

    dataVecOutVel[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
void BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::applyJT( const sofa::core::MechanicalParams* mparams,
                                                                    const sofa::type::vector< In1DataVecDeriv*>& dataVecOut1Force,
                                                                    const sofa::type::vector< In2DataVecDeriv*>& dataVecOut2Force,
                                                                    const sofa::type::vector<const OutDataVecDeriv*>& dataVecInForce)
{
    SOFA_UNUSED(mparams);

    if(d_componentState.getValue() == ComponentState::Invalid)
        return;

    if(dataVecOut1Force.empty() || dataVecInForce.empty() || dataVecOut2Force.empty())
        return;

    const OutVecDeriv& in = dataVecInForce[0]->getValue();

    In1VecDeriv& out1 = *dataVecOut1Force[0]->beginEdit();
    In2VecDeriv& out2 = *dataVecOut2Force[0]->beginEdit();

    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    const auto& directions = sofa::helper::getReadAccessor(d_directions);

    size_t sz = m_mappedPoints.size();
    for (size_t i = 0 ; i < sz; i++)
    {
        MappedPoint& mp = m_mappedPoints[i];
        if (mp.onBeam)
        {
            OutDeriv f = in[i];

            for (size_t j = 0; j < OutDeriv::total_size; j++)
            {
                f[j] = (directions[j])? f[j] : 0;
            }

            In1Deriv v;
            v.getVCenter() = mp.interpolatedTransform.getRotationMatrix() * Out::getDPos(f);
            v.getVOrientation() = mp.interpolatedTransform.getRotationMatrix() * Out::getDRot(f);

            for (size_t j = 0; j < OutDeriv::total_size; j++){
                out1[indices[i]][j] += v[j];
                out2[mp.pi1][j] -= (1 - mp.alpha) * v[j];
                out2[mp.pi2][j] -= mp.alpha * v[j];
            }
        }
    }

    dataVecOut1Force[0]->endEdit();
    dataVecOut2Force[0]->endEdit();
}


template <class TIn1, class TIn2, class TOut>
void BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::applyJT( const sofa::core::ConstraintParams* cparams,
                                                                    const sofa::type::vector< In1DataMatrixDeriv*>&  dataMatOut1Const,
                                                                    const sofa::type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
                                                                    const sofa::type::vector<const OutDataMatrixDeriv*>& dataMatInConst)
{
    SOFA_UNUSED(cparams);

    if(d_componentState.getValue() == ComponentState::Invalid)
        return;

    if(dataMatOut1Const.empty() || dataMatOut2Const.empty() || dataMatInConst.empty())
        return;

    In1MatrixDeriv& out1 = *dataMatOut1Const[0]->beginEdit();
    In2MatrixDeriv& out2 = *dataMatOut2Const[0]->beginEdit();
    const OutMatrixDeriv& in = dataMatInConst[0]->getValue();

    const auto& indices = sofa::helper::getReadAccessor(d_indices);
    const auto& directions = sofa::helper::getReadAccessor(d_directions);

    const auto& rowitEnd = in.end();
    for (auto rowIt = in.begin(); rowIt != rowitEnd; ++rowIt)
    {
        const auto& colitEnd = rowIt.end();
        for (auto colIt = rowIt.begin(); colIt != colitEnd; ++colIt)
        {
            typename In1MatrixDeriv::RowIterator o1 = out1.writeLine(rowIt.index());
            typename In2MatrixDeriv::RowIterator o2 = out2.writeLine(rowIt.index());

            MappedPoint mp = m_mappedPoints[colIt.index()];
            if (mp.onBeam)
            {
                OutDeriv h = colIt.val();

                for (size_t j = 0; j < OutDeriv::total_size; j++)
                {
                    h[j] = (directions[j])? h[j] : 0;
                }

                In1Deriv v;
                v.getVCenter() = mp.interpolatedTransform.getRotationMatrix() * Out::getDPos(h);
                v.getVOrientation() = mp.interpolatedTransform.getRotationMatrix() * Out::getDRot(h);
                In1Deriv h1, h2_1, h2_2;

                for (size_t j = 0; j < In1Deriv::total_size; j++)
                {
                    h1[j]   = v[j];
                    h2_1[j] = - (1 - mp.alpha) * v[j];
                    h2_2[j] = - mp.alpha * v[j];
                }

                o1.addCol(indices[colIt.index()], h1);
                o2.addCol(mp.pi1, h2_1);
                o2.addCol(mp.pi2, h2_2);
            }
        }
    }
    dataMatOut1Const[0]->endEdit();
    dataMatOut2Const[0]->endEdit();
}

template <class TIn1, class TIn2, class TOut>
void BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::applyDJT(const sofa::core::MechanicalParams* mparams,
                                                                      sofa::core::MultiVecDerivId inForce,
                                                                      sofa::core::ConstMultiVecDerivId outForce)
{
    SOFA_UNUSED(mparams);
    SOFA_UNUSED(inForce);
    SOFA_UNUSED(outForce);
}

template <class TIn1, class TIn2, class TOut>
const sofa::type::vector<sofa::linearalgebra::BaseMatrix*>* BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::getJs()
{
    const OutVecCoord& out = m_toModel->read(sofa::core::ConstVecCoordId::position())->getValue();
    const In1VecCoord& in1 = m_fromModel1->read(sofa::core::ConstVecCoordId::position())->getValue();
    const In2VecCoord& in2 = m_fromModel2->read(sofa::core::ConstVecCoordId::position())->getValue();

    typename SparseMatrixEigen1::CompressedMatrix& J1 = m_eigenJacobian1.compressedMatrix;
    typename SparseMatrixEigen1::CompressedMatrix& J2 = m_eigenJacobian2.compressedMatrix;
    const auto& indices = sofa::helper::getReadAccessor(d_indices);

    if( m_updateJ || J1.size() == 0 || J2.size() == 0 )
    {
        J1.resize(out.size() * NOut, in1.size() * NIn1);
        J2.resize(out.size() * NOut, in2.size() * NIn2);
        J1.setZero();
        J2.setZero();

        size_t sz = m_mappedPoints.size();
        for (size_t n = 0 ; n < sz; n++)
        {
            MappedPoint& mp = m_mappedPoints[n];
            sofa::type::Mat3x3d base = mp.interpolatedTransform.getRotationMatrix().transposed();

            Eigen::Matrix<Real, OutDeriv::total_size, OutDeriv::total_size> block;
            if (mp.onBeam)
            {
                block.template rightCols<3>() <<
                  base(0,0), base(0,1), base(0,2),
                  base(1,0), base(1,1), base(1,2),
                  base(2,0), base(2,1), base(2,2);
                block.template leftCols<NOut>().setIdentity();
            }

            for(unsigned i = 0; i < NOut; ++i)
            {
                unsigned row = n * NOut + i;

                J1.startVec( row );
                for(unsigned j = 0; j < NIn1; ++j)
                {
                    unsigned col = indices[n] * NIn1 + j;
                    J1.insertBack(row, col) = block(i, j);
                }

                J2.startVec( row );
                for(unsigned j = 0; j < NIn2; ++j)
                {
                    if (mp.onBeam)
                    {
                        unsigned col1 = mp.pi1 * NIn2 + j;
                        unsigned col2 = mp.pi2 * NIn2 + j;
                        J2.insertBack(row, col1) = - (1 - mp.alpha) * block(i, j);
                        J2.insertBack(row, col2) = - mp.alpha * block(i, j);
                    } else
                    {
                        J2.insertBack(row, j) = 0;
                        J2.insertBack(row, j) = 0;
                    }
                }
            }
        }
        J1.finalize();
        J2.finalize();
    }

    m_updateJ = false;
    return &m_eigenJacobians;
}


template <class TIn1, class TIn2, class TOut>
void BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>::draw(const sofa::core::visual::VisualParams* vparams)
{
    if(d_componentState.getValue() == ComponentState::Invalid)
        return;

    auto interpolation = l_interpolation.get();
    Real size = d_drawSize.getValue();
    const auto& x = m_fromModel2->readPositions();

    if (d_draw.getValue()){
        for (MappedPoint& mp : m_mappedPoints)
        {
            if (mp.onBeam)
            {
                Transform interpolatedTransform;
                interpolation->InterpolateTransformUsingSpline(interpolatedTransform,
                                                               mp.alpha,
                                                               Transform(x[mp.pi1].getCenter(), x[mp.pi1].getOrientation()),
                                                               Transform(x[mp.pi2].getCenter(), x[mp.pi2].getOrientation()),
                                                               interpolation->getLength(mp.edgeIndex));
                In1Coord position;
                position.getCenter() = interpolatedTransform.getOrigin();
                sofa::type::Mat3x3 base = interpolatedTransform.getRotationMatrix().transposed();
                vparams->drawTool()->drawArrow(TIn1::getCPos(position),
                                               TIn1::getCPos(position) + sofa::type::Vec3{base(0,0), base(0,1), base(0,2)}*size, size/10,
                                               sofa::type::RGBAColor::red());
                vparams->drawTool()->drawArrow(TIn1::getCPos(position),
                                               TIn1::getCPos(position) + sofa::type::Vec3{base(1,0), base(1,1), base(1,2)}*size, size/10,
                                               sofa::type::RGBAColor::green());
                vparams->drawTool()->drawArrow(TIn1::getCPos(position),
                                               TIn1::getCPos(position) + sofa::type::Vec3{base(2,0), base(2,1), base(2,2)}*size, size/10,
                                               sofa::type::RGBAColor::blue());
            }
        }
    }
}

} // namespace
