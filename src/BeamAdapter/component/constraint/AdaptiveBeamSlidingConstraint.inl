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
#pragma once

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/visual/VisualParams.h>

#include <BeamAdapter/component/constraint/AdaptiveBeamSlidingConstraint.h>
#include <sofa/core/behavior/ConstraintResolution.h>


namespace sofa::component::constraintset::_adaptiveBeamSlidingConstraint_
{

using sofa::core::behavior::ConstraintResolution ;
using sofa::core::ConstVecCoordId;
using sofa::core::ConstraintParams;
using sofa::helper::ReadAccessor;

/*!
 * \class AdaptiveBeamSlidingConstraintResolution
 * @brief AdaptiveBeamSlidingConstraintResolution Class
 */
class AdaptiveBeamSlidingConstraintResolution : public ConstraintResolution
{
public:
    AdaptiveBeamSlidingConstraintResolution(double* sliding = nullptr);

    virtual void init(int line, double** w, double* force);
    virtual void resolution(int line, double** w, double* d, double* force);
    virtual void store(int line, double* force, bool convergence);

    void resolution(int line, double** w, double* d, double* force, double* dfree);

protected:
    double* m_slidingDisp;
    double  m_slidingW;
};


template<class DataTypes>
AdaptiveBeamSlidingConstraint<DataTypes>::AdaptiveBeamSlidingConstraint(TypedMechanicalState* object1, TypedMechanicalState* object2)
    : Inherit(object1, object2)
    , m_interpolation(initLink("interpolation", "link to the WireBeamInterpolation component in the scene"))
{
}

template<class DataTypes>
AdaptiveBeamSlidingConstraint<DataTypes>::AdaptiveBeamSlidingConstraint(TypedMechanicalState* object)
    : AdaptiveBeamSlidingConstraint(object, object)
{
}

template<class DataTypes>
AdaptiveBeamSlidingConstraint<DataTypes>::AdaptiveBeamSlidingConstraint()
    : AdaptiveBeamSlidingConstraint(nullptr, nullptr)
{
}

template<class DataTypes>
void AdaptiveBeamSlidingConstraint<DataTypes>::init()
{
    PairInteractionConstraint<DataTypes>::init();

    if(mstate1 == nullptr || mstate2 == nullptr)
        msg_error() << "This component needs to be linked to two MechanicalObject.";

    if(!m_interpolation.get())
        msg_error() << "This component needs to be linked to a WireBeamInterpolation component.";

    internalInit();
}

template<class DataTypes>
void AdaptiveBeamSlidingConstraint<DataTypes>::reset()
{
    internalInit();
}

template<class DataTypes>
void AdaptiveBeamSlidingConstraint<DataTypes>::internalInit()
{
    // We search for the closest segment, on which to project each point
    // Convention : object1 is the beam model, object2 is the list of point constraints

    ReadAccessor<Data<VecCoord> > x1 = mstate1->read(ConstVecCoordId::position()) ;
    ReadAccessor<Data<VecCoord> > x2 = mstate2->read(ConstVecCoordId::position()) ;

    unsigned int m2Size = x2.size();
    m_previousPositions.clear();
    m_previousPositions.resize(m2Size);
    m_projected.clear();
    m_projected.resize(m2Size);
    m_displacements.clear();
    m_displacements.resize(m2Size);

    WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();

    for(unsigned int i=0; i<m2Size; i++)
    {
        Real r = -1;
        Vec3 pt = x2[i].getCenter();
        bool p = interpolation->getApproximateCurvAbs(pt, x1.ref(), r);

        m_previousPositions[i] = r;
        m_projected[i] = p;
    }
}


template<class DataTypes>
bool AdaptiveBeamSlidingConstraint<DataTypes>::getCurvAbsOfProjection(const Vec3& x_input, const VecCoord& vecX, Real& xcurv_output, const Real& tolerance)
{
    WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();

    // We have put all the code in a new class, because it uses a lot of custom functions and data
    ProjectionSearch<DataTypes> ps(interpolation, x_input, vecX, xcurv_output, tolerance);
    return ps.doSearch(xcurv_output);
}


template<class DataTypes>
void AdaptiveBeamSlidingConstraint<DataTypes>::buildConstraintMatrix(const ConstraintParams * cParams,
                                                              DataMatrixDeriv &c1_d, DataMatrixDeriv &c2_d,
                                                              unsigned int &constraintId,
                                                              const DataVecCoord &x1_d, const DataVecCoord &x2_d)
{
    SOFA_UNUSED(cParams);

    m_violations.clear();
    m_nbConstraints = 0;
    m_cid = constraintId;

    Transform Tnode0, Tnode1, Tresult;
    Real baryCoord;
    unsigned int beam = 0;

    ReadAccessor<Data<VecCoord> > x1free=mstate1->read(ConstVecCoordId::freePosition()) ;
    ReadAccessor<Data<VecCoord> > x2free=mstate2->read(ConstVecCoordId::freePosition()) ;

    unsigned int m2 = x2free.size();
    WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();
    auto c1 = sofa::helper::getWriteOnlyAccessor(c1_d);
    auto c2 = sofa::helper::getWriteOnlyAccessor(c2_d);
    const VecCoord& x1= x1_d.getValue();
    const VecCoord& x2= x2_d.getValue();


    for(unsigned int i=0; i<m2 ; i++)
    {
        if(!m_projected[i])
            continue;

        // Get new projection on the curve
        m_previousPositions[i] += (Real) m_displacements[i];
        if(!this->getCurvAbsOfProjection(x2[i].getCenter(), x1, m_previousPositions[i], 1e-5))
        {
            m_projected[i] = false;
            continue;
        }

        unsigned int node0, node1;
        interpolation->getNodeIndices(beam, node0, node1);

        // Position and frame on the curve
        interpolation->getBeamAtCurvAbs(m_previousPositions[i], beam, baryCoord);
        interpolation->computeTransform(beam, node0, node1, Tnode0, Tnode1, x1free.ref());
        interpolation->InterpolateTransformUsingSpline(Tresult, baryCoord, Tnode0, Tnode1, interpolation->getLength(beam));
        Pos p = Tresult.getOrigin();
        Pos dir0, dir1, dir2;
        Rot rot = Tresult.getOrientation();
        dir0 = rot.rotate(Pos(1,0,0));
        dir1 = rot.rotate(Pos(0,1,0));
        dir2 = rot.rotate(Pos(0,0,1));

        // Compute violations
        Pos violation = p - x2free[i].getCenter();
        m_violations.push_back(violation * dir1);
        m_violations.push_back(violation * dir2);
        m_violations.push_back(violation * dir0);

        // Define the constraint
        SpatialVector sv0, sv1;
        Vec3 nullRot(0,0,0);

        MatrixDerivRowIterator c1_it = c1.wref().writeLine(m_cid + m_nbConstraints);
        MatrixDerivRowIterator c2_it = c2.wref().writeLine(m_cid + m_nbConstraints);
        interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir1, sv0, sv1);
        c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
        c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
        c2_it.addCol(i, Deriv(-dir1, nullRot));
        m_nbConstraints++;

        c1_it = c1.wref().writeLine(m_cid + m_nbConstraints);
        c2_it = c2.wref().writeLine(m_cid + m_nbConstraints);
        interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir2, sv0, sv1);
        c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
        c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
        c2_it.addCol(i, Deriv(-dir2, nullRot));
        m_nbConstraints++;

        c1_it = c1.wref().writeLine(m_cid + m_nbConstraints);
        c2_it = c2.wref().writeLine(m_cid + m_nbConstraints);
        interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir0, sv0, sv1);
        c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
        c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
        c2_it.addCol(i, Deriv(-dir0, nullRot));
        m_nbConstraints++;
    }

    constraintId += m_nbConstraints;
}


template<class DataTypes>
void AdaptiveBeamSlidingConstraint<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                               BaseVector *v, const
                                                               DataVecCoord &x1, const DataVecCoord &x2,
                                                               const DataVecDeriv &v1, const DataVecDeriv &v2)
{
    SOFA_UNUSED(cParams);
    SOFA_UNUSED(x1);
    SOFA_UNUSED(x2);
    SOFA_UNUSED(v1);
    SOFA_UNUSED(v2);

    unsigned int nb = m_violations.size();
    for(unsigned int i=0; i<nb; i++)
        v->set(m_cid+i, m_violations[i]);
}


template<class DataTypes>
void AdaptiveBeamSlidingConstraint<DataTypes>::getConstraintResolution(const ConstraintParams* cParams,
                                                                std::vector<ConstraintResolution*>& resTab,
                                                                unsigned int& offset)
{
    SOFA_UNUSED(cParams);

    unsigned int nb = mstate2->getSize();
    for(unsigned int i=0; i<nb; i++)
    {
        if(!m_projected[i]) continue;

        resTab[offset] = new AdaptiveBeamSlidingConstraintResolution(&m_displacements[i]);
        offset += 3;
    }
}


template<class DataTypes>
void AdaptiveBeamSlidingConstraint<DataTypes>::draw(const VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowInteractionForceFields()) return;
    
    vparams->drawTool()->saveLastState();

    ReadAccessor<Data<VecCoord> > x = mstate2->read(ConstVecCoordId::position());
    sofa::type::Vec3 point;
    std::vector< sofa::type::Vec3 > points;
    std::vector< sofa::type::RGBAColor> colors;

    points.reserve(x.size());
    colors.reserve(x.size());

    for (sofa::Size i = 0; i < x.size(); i++)
    {
        point = DataTypes::getCPos(x[i]);
        points.push_back(point);
        
        if (m_projected[i])
            colors.emplace_back(0.0f, 1.0f, 0.0f, 1.0f);
        else
            colors.emplace_back(0.0f, 1.0f, 1.0f, 1.0f);
    }
    vparams->drawTool()->drawPoints(points, 10.f, colors);

    vparams->drawTool()->restoreLastState();
}


template<class DataTypes>
bool ProjectionSearch<DataTypes>::doSearch(Real& result)
{
    initSearch(m_e);

    // Do one pass of the newton method
    newtonMethod();

    Real dist = computeDistAtCurvAbs(m_e);

    // If the new estimate is good, go with it
    if (m_found)
    {
        result = m_e;
        return testForProjection(result);
    }

    // Look if the new estimate is closer to the target than the previous estimate
    if (dist < m_de)
    {
        // If it is, continue with the newton method, it should converge fast
        while (m_totalIterations < m_interpolation->getNumBeams() + 1)
        {
            m_de = dist;
            newtonMethod();
            m_totalIterations++;

            if (m_found)
            {
                // Go back to global curv abs
                result = m_beamStart + (m_beamEnd - m_beamStart) * m_le;
                return testForProjection(result);
            }

            dist = computeDistAtCurvAbs(m_e);
            if (dist > m_de)	// Diverging
                break;

            if (m_le < 0 || m_le > 1)	// We have to change beam, use the dichotomic method instead
                break;
        }
    }

    // If the estimate is outside the beam, or is further than the previous one, change beam
    //  and use a dichotomic search until we find a solution
    if (!m_found)
    {
        m_totalIterations = 0;

        while (m_totalIterations < m_interpolation->getNumBeams() + 10 && m_dichotomicIterations < 10)
        {
            m_totalIterations++;
            m_dichotomicIterations++;

            // We will compute 10 samples
            range = m_segEnd - m_segStart;
            rangeSampling = range / s_sampling;
            for (unsigned int i = 0; i <= s_sampling; i++)
            {
                Real curvAbs = m_segStart + rangeSampling * (Real)i;
                m_distTab[i] = computeDistAtCurvAbs(curvAbs);
            }
            unsigned int minIndex = std::min_element(m_distTab, m_distTab + s_sampling + 1) - m_distTab;	// Where is the minimum
            m_e = m_segStart + rangeSampling * minIndex;
            if (testForProjection(m_e))
            {
                result = m_e;
                return true;
            }

            // If the minium is at one extremity, change beam
            if (minIndex == 0)
                changeCurrentBeam(m_beamIndex - 1);
            else if (minIndex == s_sampling)
                changeCurrentBeam(m_beamIndex + 1);
            else
            {	// We continue the search with a smaller interval (keeping only 3 points)
                m_segStart = m_e - rangeSampling;
                m_segEnd = m_e + rangeSampling;
            }
        }
    }

    result = m_e;
    return testForProjection(result);
}

template<class DataTypes>
void ProjectionSearch<DataTypes>::initSearch(Real curvAbs)
{
    m_e = curvAbs;
    if (m_e < 0)
    {
        m_beamIndex = 0;
        m_le = 0;
    }
    else if (m_e > m_interpolation->getRestTotalLength())
    {
        m_beamIndex = m_interpolation->getNumBeams() - 1;
        m_le = 1;
    }
    else
        m_interpolation->getBeamAtCurvAbs(m_e, m_beamIndex, m_le);

    m_searchDirection = 0;
    m_interpolation->getAbsCurvXFromBeam(m_beamIndex, m_beamStart, m_beamEnd);
    m_segStart = m_beamStart;
    m_segEnd = m_beamEnd;
    m_interpolation->getSplinePoints(m_beamIndex, m_x, P0, P1, P2, P3);
    m_de = computeDistAtCurvAbs(m_e);
}

template<class DataTypes>
void ProjectionSearch<DataTypes>::newtonMethod()
{
    Real bx = m_le, bx2 = bx * bx, bx3 = bx2 * bx;
    Real obx = 1 - bx, obx2 = obx * obx, obx3 = obx2 * obx;
    Vec3 P = P0 * obx3 + P1 * 3 * bx * obx2 + P2 * 3 * bx2 * obx + P3 * bx3;
    Vec3 dP = P0 * (-3 * obx2) + P1 * (3 - 12 * bx + 9 * bx2) + P2 * (6 * bx - 9 * bx2) + P3 * (3 * bx2);
    Real f_x = dot((m_target - P), dP);

    if (f_x == 0.0 || fabs(f_x) / dP.norm() < m_tolerance) // reach convergence
    {
        m_interpolation->getCurvAbsAtBeam(m_beamIndex, m_le, m_e);
        m_found = true;
        return;
    }
    Vec3 d2P = P0 * 6 * (1 - bx) + P1 * (-12 + 18 * bx) + P2 * (6 - 18 * bx) + P3 * 6 * bx;

    Real df_x = dot(-dP, dP) + dot((m_target - P), d2P);

    if (fabs(df_x) < 1e-5 * m_tolerance)
    {
        return;
    }

    Real d_bx = -f_x / df_x;
    m_le += d_bx;

    m_e = m_beamStart + (m_beamEnd - m_beamStart) * m_le;

    // NOTE : bx+d_bx-1.0 ne donne pas une estimation correcte de la position dans l'autre beam, puisque sa longueur peut �tre diff�rente !
}

template<class DataTypes>
bool ProjectionSearch<DataTypes>::changeCurrentBeam(int index)
{
    // If at the end of the thread
    if (index < 0)
    {
        m_segEnd = m_segStart + rangeSampling;
        return false;
    }
    else if (index > static_cast<int>(m_interpolation->getNumBeams()) - 1)
    {
        m_segStart = m_segEnd - rangeSampling;
        return false;
    }

    int dir = index - m_beamIndex;
    if (m_searchDirection * dir < 0)
    {	// Changing the direction of search means we are looking for a point near an extremity
        // We know we are looking for a point inside the interval [0.9;1.1]
        // but the ranges for each beam can be different
        Real nStart, nEnd;
        m_interpolation->getAbsCurvXFromBeam(index, nStart, nEnd);
        Real nRangeSampling = (nEnd - nStart) / s_sampling;

        if (dir < 0)
        {
            m_segStart = m_beamStart - nRangeSampling;
            m_segEnd = m_beamStart + rangeSampling;
        }
        else
        {
            m_segStart = m_beamEnd - rangeSampling;
            m_segEnd = m_beamEnd + nRangeSampling;
        }
        return false;
    }


    if (dir < 0)
        m_searchDirection = -1;
    else if (dir > 0)
        m_searchDirection = 1;

    // Really changing beam
    m_beamIndex = index;
    m_interpolation->getAbsCurvXFromBeam(m_beamIndex, m_beamStart, m_beamEnd);
    m_segStart = m_beamStart;
    m_segEnd = m_beamEnd;
    m_interpolation->getSplinePoints(m_beamIndex, m_x, P0, P1, P2, P3);
    m_dichotomicIterations = 0;

    return true;
}

template<class DataTypes>
typename ProjectionSearch<DataTypes>::Real ProjectionSearch<DataTypes>::computeDistAtCurvAbs(Real curvAbs)
{
    if (curvAbs >= m_beamStart && curvAbs <= m_beamEnd)
    {	// We can use the control points we saved
        Real bx = (curvAbs - m_beamStart) / (m_beamEnd - m_beamStart), bx2 = bx * bx, bx3 = bx2 * bx;
        Real obx = 1 - bx, obx2 = obx * obx, obx3 = obx2 * obx;
        Vec3 P = P0 * obx3 + P1 * 3 * bx * obx2 + P2 * 3 * bx2 * obx + P3 * bx3;

        return (m_target - P).norm();
    }
    else
    {
        // TODO(dmarchal 2017-05-17) Please specify who/when this will be done
        // TODO : save all the control points so we don't have to compute them again
        Real bx;
        unsigned int index;
        Vec3 tP0, tP1, tP2, tP3;
        m_interpolation->getBeamAtCurvAbs(curvAbs, index, bx);
        m_interpolation->getSplinePoints(index, m_x, tP0, tP1, tP2, tP3);
        Real bx2 = bx * bx, bx3 = bx2 * bx;
        Real obx = 1 - bx, obx2 = obx * obx, obx3 = obx2 * obx;
        Vec3 P = tP0 * obx3 + tP1 * 3 * bx * obx2 + tP2 * 3 * bx2 * obx + tP3 * bx3;

        return (m_target - P).norm();
    }
}

template<class DataTypes>
bool ProjectionSearch<DataTypes>::testForProjection(Real curvAbs)
{
    Real bx;
    unsigned int index;
    Vec3 tP0, tP1, tP2, tP3;
    m_interpolation->getBeamAtCurvAbs(curvAbs, index, bx);
    m_interpolation->getSplinePoints(index, m_x, tP0, tP1, tP2, tP3);
    Real bx2 = bx * bx, bx3 = bx2 * bx;
    Real obx = 1 - bx, obx2 = obx * obx, obx3 = obx2 * obx;
    Vec3 P = tP0 * obx3 + tP1 * 3 * bx * obx2 + tP2 * 3 * bx2 * obx + tP3 * bx3;
    Vec3 dP = tP0 * (-3 * obx2) + tP1 * (3 - 12 * bx + 9 * bx2) + tP2 * (6 * bx - 9 * bx2) + tP3 * (3 * bx2);
    Real f_x = dot((m_target - P), dP);

    if (fabs(f_x) / dP.norm() < m_tolerance)
        return true;

    return false;
}


} // namespace sofa::component::constraintset::_adaptiveBeamSlidingConstraint_


