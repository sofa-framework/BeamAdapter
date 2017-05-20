/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 RC 1        *
*                (c) 2006-2011 MGH, INRIA, USTL, UJF, CNRS                    *
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
#ifndef SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMCONSTRAINT_INL
#define SOFA_COMPONENT_CONSTRAINTSET_ADAPTIVEBEAMCONSTRAINT_INL


//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/core/visual/VisualParams.h>
#include <sofa/helper/gl/template.h>

#include "AdaptiveBeamConstraint.h"

namespace sofa
{

namespace component
{

namespace constraintset
{

namespace _adaptivebeamconstraint_
{

using std::vector;
using sofa::core::behavior::ConstraintResolution ;
using sofa::core::ConstVecCoordId;
using sofa::core::ConstraintParams;
using sofa::helper::ReadAccessor;

/*!
 * \class AdaptiveBeamConstraintResolution
 * @brief AdaptiveBeamConstraintResolution Class
 */
class AdaptiveBeamConstraintResolution : public ConstraintResolution
{
public:
    AdaptiveBeamConstraintResolution(double* sliding = nullptr);

    virtual void init(int line, double** w, double* force);
    virtual void resolution(int line, double** w, double* d, double* force);
    virtual void store(int line, double* force, bool convergence);

    void resolution(int line, double** w, double* d, double* force, double* dfree);

protected:
    double* m_slidingDisp;
    double  m_slidingW;
};


template<class DataTypes>
AdaptiveBeamConstraint<DataTypes>::AdaptiveBeamConstraint(TypedMechanicalState* object1, TypedMechanicalState* object2)
    : Inherit(object1, object2)
    , m_interpolation(initLink("interpolation", "link to the interpolation component in the scene"))
{
}

template<class DataTypes>
AdaptiveBeamConstraint<DataTypes>::AdaptiveBeamConstraint(TypedMechanicalState* object)
    : AdaptiveBeamConstraint(object, object)
{
}

template<class DataTypes>
AdaptiveBeamConstraint<DataTypes>::AdaptiveBeamConstraint()
    : AdaptiveBeamConstraint(nullptr, nullptr)
{
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::init()
{
    assert(mstate1);
    assert(mstate2);
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::reset()
{
    internalInit();
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::internalInit()
{
    // TODO(eulalie): I'm not sure that this comment belongs here?
    // We search for the closest segment, on which to project each point
    // Convention : object1 is the beam model, object2 is the list of point constraints

    if(!m_interpolation.get())
        return;

    ReadAccessor<Data<VecCoord> > x1 = mstate1->read(ConstVecCoordId::position()) ;
    ReadAccessor<Data<VecCoord> > x2 = mstate2->read(ConstVecCoordId::position()) ;

    unsigned int m2 = x2.size();
    m_previousPositions.clear();
    m_previousPositions.resize(m2);
    m_projected.clear();
    m_projected.resize(m2);
    m_displacements.clear();
    m_displacements.resize(m2);

    WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();

    for(unsigned int i=0; i<m2; i++)
    {
        Real r = -1;
        Vec3 pt = x2[i].getCenter();
        bool p = interpolation->getApproximateCurvAbs(pt, x1.ref(), r);

        m_previousPositions[i] = r;
        m_projected[i] = p;
    }
}

template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::buildConstraintMatrix(const ConstraintParams * cParams,
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
    ReadAccessor<Data<VecCoord> > x2free=mstate1->read(ConstVecCoordId::freePosition()) ;

    unsigned int m2 = x2free.size();
    WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();
    MatrixDeriv& c1 = *c1_d.beginEdit();
    MatrixDeriv& c2 = *c2_d.beginEdit();
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
        Pos dir, dir1, dir2;
        Rot rot = Tresult.getOrientation();
        dir = rot.rotate(Pos(1,0,0));
        dir1 = rot.rotate(Pos(0,1,0));
        dir2 = rot.rotate(Pos(0,0,1));

        // Compute violations
        Pos violation = p - x2free[i].getCenter();
        m_violations.push_back(violation * dir1);
        m_violations.push_back(violation * dir2);
        m_violations.push_back(violation * dir);

        // Define the constraint
        unsigned int node0, node1;
        SpatialVector sv0, sv1;
        Vec3 nullRot(0,0,0);
        interpolation->getNodeIndices(beam, node0, node1);

        MatrixDerivRowIterator c1_it = c1.writeLine(m_cid + m_nbConstraints);
        MatrixDerivRowIterator c2_it = c2.writeLine(m_cid + m_nbConstraints);
        interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir1, sv0, sv1);
        c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
        c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
        c2_it.addCol(i, Deriv(-dir1, nullRot));
        m_nbConstraints++;

        c1_it = c1.writeLine(m_cid + m_nbConstraints);
        c2_it = c2.writeLine(m_cid + m_nbConstraints);
        interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir2, sv0, sv1);
        c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
        c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
        c2_it.addCol(i, Deriv(-dir2, nullRot));
        m_nbConstraints++;

        c1_it = c1.writeLine(m_cid + m_nbConstraints);
        c2_it = c2.writeLine(m_cid + m_nbConstraints);
        interpolation->MapForceOnNodeUsingSpline(beam, baryCoord, Pos(0,0,0), x1, dir, sv0, sv1);
        c1_it.addCol(node0, Vec6(sv0.getForce(), sv0.getTorque()));
        c1_it.addCol(node1, Vec6(sv1.getForce(), sv1.getTorque()));
        c2_it.addCol(i, Deriv(-dir, nullRot));
        m_nbConstraints++;
    }

    constraintId += m_nbConstraints;
    c1_d.endEdit();
    c2_d.endEdit();
}


template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::getConstraintViolation(const ConstraintParams* cParams,
                                                               defaulttype::BaseVector *v, const
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
void AdaptiveBeamConstraint<DataTypes>::getConstraintResolution(vector<ConstraintResolution*>& resTab,
                                                                unsigned int& offset)
{
    unsigned int nb = mstate2->getSize();
    for(unsigned int i=0; i<nb; i++)
    {
        if(!m_projected[i]) continue;

        resTab[offset] = new AdaptiveBeamConstraintResolution(&m_displacements[i]);
        offset += 3;
    }
}


template<class DataTypes>
void AdaptiveBeamConstraint<DataTypes>::draw(const VisualParams* vparams)
{
#ifndef SOFA_NO_OPENGL
    if(!vparams->displayFlags().getShowInteractionForceFields()) return;

    glDisable(GL_LIGHTING);
    glPointSize(10);
    glBegin(GL_POINTS);

    ReadAccessor<Data<VecCoord> > x = mstate2->read(ConstVecCoordId::position()) ;
    unsigned int m = x.size();

    for(unsigned int i=0; i<m; i++)
    {
        glColor4f(0.0f,1.0f,m_projected[i]?1:0.0f,1.0f);
        helper::gl::glVertexT(x[i]);
    }

    glEnd();
    glPointSize(1);
#endif // SOFA_NO_OPENGL
}

} // namespace _adaptivebeamconstraint_

} // namespace constraintset

} // namespace component

} // namespace sofa

#endif
