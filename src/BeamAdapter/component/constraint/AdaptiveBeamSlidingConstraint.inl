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

namespace sofa::component::constraintset
{

namespace _adaptiveBeamSlidingConstraint_
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
    unsigned int beam;

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
        if(!interpolation->getCurvAbsOfProjection(x2[i].getCenter(), x1, m_previousPositions[i], 1e-5))
        {
            m_projected[i] = false;
            continue;
        }

        // Position and frame on the curve
        interpolation->getBeamAtCurvAbs(m_previousPositions[i], beam, baryCoord);
        interpolation->computeTransform2(beam, Tnode0, Tnode1, x1free.ref());
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
        unsigned int node0, node1;
        SpatialVector sv0, sv1;
        Vec3 nullRot(0,0,0);
        interpolation->getNodeIndices(beam, node0, node1);

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

    for (auto i = 0; i < x.size(); i++)
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

} // namespace _adaptiveBeamSlidingConstraint_

} // namespace sofa::component::constraintset
