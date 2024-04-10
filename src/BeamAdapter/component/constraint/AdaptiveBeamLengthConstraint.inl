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
#include <sofa/type/Vec.h>
#include <sofa/linearalgebra/BaseVector.h>
#include <sofa/helper/visual/DrawTool.h>
#include <sofa/core/behavior/ConstraintResolution.h>
#include <BeamAdapter/component/constraint/AdaptiveBeamLengthConstraint.h>


namespace sofa::component::constraintset::_adaptivebeamlengthconstraint_
{

using helper::ReadAccessor;
using sofa::core::ConstVecCoordId;
using std::stringstream;
using sofa::core::ConstraintParams;
using sofa::linearalgebra::BaseVector;
using sofa::core::visual::VisualParams;


class AdaptiveBeamLengthConstraintResolution : public ConstraintResolution
{
public:
    AdaptiveBeamLengthConstraintResolution(SReal* initF=nullptr, bool* active=nullptr) : ConstraintResolution(1) ,m_initF(initF), m_active(active)
    {
    }
    virtual void init(int line, SReal** w, SReal* force) override;
    virtual void store(int line, SReal* force, bool convergence) override;

    using sofa::core::behavior::ConstraintResolution::resolution;
    void resolution(int line, SReal** w, SReal* d, SReal* force);
    
protected:
    SReal*    m_initF;
    bool*      m_active;
};


template<class DataTypes>
AdaptiveBeamLengthConstraint<DataTypes>::AdaptiveBeamLengthConstraint(TypedMechanicalState* object)
    : Inherit(object)
    , m_alarmLength(initData(&m_alarmLength, (Real)1.02, "alarmLength", "Elongation before creating a constraint (default=1.02)"))
    , m_constrainedLength(initData(&m_constrainedLength, (Real)1.05, "constrainedLength", "Allowed elongation of a beam (default=1.05"))
    , m_maxBendingAngle(initData(&m_maxBendingAngle,  (Real)0.1, "maxBendingAngle", "max bending criterion (in rad) for one constraint interval (default=0.1)"))
    , m_interpolation(initLink("interpolation", "link to the interpolation component in the scene"))
{
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::init()
{
    this->mstate= dynamic_cast< MechanicalState<DataTypes> *> (this->getContext()->getMechanicalState());
    assert(this->mstate);
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::reset()
{
    internalInit();
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::internalInit()
{
    /// We search for the closest segment, on which to project each point
    /// Convention : object1 is the beam model, object2 is the list of point constraints
    if(!m_interpolation.get())
        return;
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::detectElongation(const VecCoord& x, const VecCoord& xfree)
{
    Vec3 P0,P1,P2,P3;
    Real length, rest_length, rest_length_interval;
    Real alarmLength = m_alarmLength.getValue();
    bool prev_stretch = false;

    fem::WireBeamInterpolation<DataTypes>* interpolation = m_interpolation.get();

    /// storage of the length (and the rest_length)  of the interval being stretched
    ///	length_interval=0.0; //commented to remove compilation warning
    rest_length_interval=0.0;

    /// storage of the interval information
    IntervalDefinition<Real> intervalDef;
    Real angleInterval=0.0;

    for (unsigned int b=0; b<interpolation->getNumBeams(); b++)
    {
        /// 1. compute the actual length and the rest length of the beams
        interpolation->getSplinePoints(b,x,P0,P1,P2,P3);

        //TODO(dmarchal 2017) Please specify who/when this will be done
        /// TODO : optimization: finally it is not necessary to compute all the spline points
        length=(P0-P3).norm();
        rest_length = interpolation->getLength(b);

        unsigned n0, n1;
        interpolation->getNodeIndices(b, n0, n1);

        /// 2. compute the bending angle
        Transform Tnode0, Tnode1;
        interpolation->computeTransform(b, n0, n1, Tnode0,Tnode1,x);
        Real angleBeam = interpolation->ComputeTotalBendingRotationAngle(rest_length/10.0, Tnode0, Tnode1,rest_length , 0.0, 1.0);

        /// 3. treatment of the different case..
        bool case1a = (n0==n1);

        if(prev_stretch)
        {
            ////CASE 1: previous beam was stretched:
            /// find the case that necessitates to stop the interval:
            /// (a) rigidification
            /// (b) current beam not stretched
            /// (c) too large bending angle...
            /// (d) last beam ! [ see after the loop "for" ]
            //// => store the information : index [begin end], Pos [begin end]

            bool case1b = (rest_length*alarmLength > length);
            bool case1c = (angleBeam >= 0.99*m_maxBendingAngle.getValue());
            case1c=false;
            case1b=false;

            /// CASE 1 (a) + (b) + (c)
            if (case1a || case1b || case1c)
            {

                dmsg_info() <<" beam "<<b<<" case 1 detected (a):"<<n0<<" == ?"<<n1<<"  (b) :"<<
                               rest_length*alarmLength<<" > ?"<<length<<"  (c) : "<<angleBeam+angleInterval<< " > ?" <<m_maxBendingAngle.getValue() ;

                /// Stop the rigidification
                /// the interval ends at the beginning of the current beam
                intervalDef.posEnd = P0; /// store the position
                intervalDef.IdxEnd = n0; /// store the index of the dof

                Transform DOF0_H_local0, DOF1_H_local1;
                interpolation->getDOFtoLocalTransform(b, DOF0_H_local0,  DOF1_H_local1);

                intervalDef.dof_H_end = DOF0_H_local0; /// store the transform from dof to pos

                interpolation->getSplinePoints(b,xfree,P0,P1,P2,P3);

                Transform global_H_local0_free, global_H_local1_free;
                interpolation->computeTransform(b, n0, n1, global_H_local0_free,  global_H_local1_free, xfree);
                intervalDef.posFreeEnd  = global_H_local0_free.getOrigin(); /// store the free position

                intervalDef.rest_length=rest_length_interval; /// store the rest_length

                dmsg_info() << " rest_length_interval ="<<rest_length_interval ;

                /// ends the interval
                prev_stretch=false;
                angleInterval=0.0;
                rest_length_interval=0.0;

                /// verify that the interval length is not null:
                if ((intervalDef.posBegin-intervalDef.posEnd).norm() < 0.0000001)
                    msg_warning() <<"interval of size = 0 detected => 1 beam has a bendingAngle > m_maxBendingAngle" ;
                else
                    m_constraintIntervals.push_back(intervalDef);

                /// isolate the case 1 (c)
                if (!case1a && !case1b && case1c) //
                {
                    dmsg_info() <<" isolate case 1(c)" ;

                    /// create a new interval:
                    prev_stretch=true;
                    angleInterval = angleBeam;
                    rest_length_interval = rest_length;

                    /// store the information of the beginning of the new interval
                    /// -> it corresponds to the end of the previous interval
                    intervalDef.dof_H_begin = DOF0_H_local0;
                    intervalDef.IdxBegin = n0;
                    intervalDef.posBegin = intervalDef.posEnd;
                    intervalDef.posFreeBegin = intervalDef.posFreeEnd;
                }
            }
            else
            {
                /// Continue the rigidification
                angleInterval+=angleBeam;
                rest_length_interval += rest_length;
            }
        }
        else
        {
            //// CASE 2: previous beam was not stretched:
            /// (a) the current beam is stretched: start the interval => store index_begin Pos_begin
            ///     veriy the bending angle (in case the interval= the beam)
            /// (b) the current beam is not stretched=> nothing to do !
            bool case2a = (rest_length*alarmLength < length);
            case2a=true;

            if (  case2a && !case1a ) /// CASE 2 (a)
            {
                dmsg_info() << " beam "<<b<<" case 2 (a) detected "<<rest_length*alarmLength<<" < ?"<<length ;

                /// create a new interval:
                prev_stretch=true;
                angleInterval = angleBeam;
                rest_length_interval = rest_length;

                /// store the information of the beginning of the new interval
                /// -> it corresponds to the position of node 0 of the beam
                Transform DOF0_H_local0, DOF1_H_local1;
                interpolation->getDOFtoLocalTransform(b, DOF0_H_local0,  DOF1_H_local1);

                Transform global_H_local0, global_H_local1;
                interpolation->computeTransform(b, n0, n1, global_H_local0,  global_H_local1, x);

                Transform global_H_local0_free, global_H_local1_free;
                interpolation->computeTransform(b, n0, n1, global_H_local0_free,  global_H_local1_free, xfree);

                intervalDef.dof_H_begin = DOF0_H_local0;
                intervalDef.IdxBegin = n0;
                intervalDef.posBegin = global_H_local0.getOrigin();
                intervalDef.posFreeBegin = global_H_local0_free.getOrigin();
            }
            else
            {
                dmsg_info() <<" beam "<<b<<" case 2 (b) detected "<<rest_length*alarmLength<<" > ?"<<length<<" or n0="<<n0<<" ==? "<<"n1="<<n1 ;
            }
        }
    }

    unsigned int b = interpolation->getNumBeams()-1;
    unsigned n0, n1;
    interpolation->getNodeIndices(b, n0, n1);
    if(prev_stretch) /// case 1(d)
    {
        dmsg_info() <<" case 1 (d) detected on the last beam" ;

        /// store the information of the beginning of the new interval
        /// -> it corresponds to the position of node 0 of the beam
        Transform DOF0_H_local0, DOF1_H_local1;
        interpolation->getDOFtoLocalTransform(b, DOF0_H_local0,  DOF1_H_local1);

        Transform global_H_local0, global_H_local1;

        interpolation->computeTransform(b, n0, n1, global_H_local0,  global_H_local1, x);


        Transform global_H_local0_free, global_H_local1_free;
        interpolation->computeTransform(b, n0, n1, global_H_local0_free,  global_H_local1_free, xfree);


        intervalDef.dof_H_end = DOF1_H_local1;
        intervalDef.IdxEnd = n1;
        intervalDef.posEnd = global_H_local1.getOrigin();
        intervalDef.posFreeEnd = global_H_local1_free.getOrigin();
        intervalDef.rest_length=rest_length_interval; /// store the rest_length

        dmsg_info() <<" rest_length_interval ="<<rest_length_interval ;

        m_constraintIntervals.push_back(intervalDef);
    }
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::buildConstraintMatrix(const ConstraintParams* cParams, DataMatrixDeriv &c_d,
                                                                    unsigned int &constraintId, const DataVecCoord & x_d)
{
    SOFA_UNUSED(cParams) ;
    SOFA_UNUSED(x_d) ;

    m_violations.clear();
    Real constrainedLength = m_constrainedLength.getValue();

    m_nbConstraints = 0;
    m_cid = constraintId;

    ReadAccessor<Data<VecCoord> > x = this->mstate->read(ConstVecCoordId::position()) ;
    ReadAccessor<Data<VecCoord> > xfree = this->mstate->read(ConstVecCoordId::freePosition()) ;

    auto c = sofa::helper::getWriteOnlyAccessor(c_d);

    m_constraintIntervals.clear();
    detectElongation( x.ref(), xfree.ref());

    if( this->f_printLog.getValue())
    {
        stringstream tmp;
        for (unsigned int i=0; i<m_constraintIntervals.size();i++)
        {
            tmp <<"constraint["<<i<<"] between pos: "<<m_constraintIntervals[i].posBegin<<" and pos: "<<m_constraintIntervals[i].posEnd << msgendl ;
        }
        msg_info() << tmp.str() ;
    }

    /// for each constraint Interval, a constraint is created
    for (unsigned int i=0; i<m_constraintIntervals.size();i++)
    {
        /// get the indices of the dofs involved
        unsigned int n0,n1;
        n0=m_constraintIntervals[i].IdxBegin;
        n1=m_constraintIntervals[i].IdxEnd;

        /// compute the violation and the local direction of the constraint
        Vec3 dir = m_constraintIntervals[i].posEnd - m_constraintIntervals[i].posBegin;
        Vec3 PbPeFree= m_constraintIntervals[i].posFreeEnd - m_constraintIntervals[i].posFreeBegin;

        dir.normalize();
        Real length_free= dot(dir, PbPeFree);

        /// put the violation in a buffer
        m_violations.push_back(m_constraintIntervals[i].rest_length * constrainedLength - length_free);

        /// project the  direction of the constraint (considered as a force) in the DOF frame
        Vec3 lever0 = x[n0].getOrientation().rotate( m_constraintIntervals[i].dof_H_begin.getOrigin() );
        Vec3 lever1 = x[n1].getOrientation().rotate( m_constraintIntervals[i].dof_H_end.getOrigin()  );

        dmsg_info_when(lever0.norm() > 0.0001) << " lever0 ="<<lever0<<" dir ="<<dir ;
        dmsg_info_when(lever1.norm() > 0.0001) << " lever1 ="<<lever1<<" -dir ="<<-dir ;

        MatrixDerivRowIterator c_it = c.wref().writeLine(m_cid + m_nbConstraints);

        c_it.addCol(n0, Vec6(dir, cross(lever0, dir) ) );
        c_it.addCol(n1, Vec6(-dir, cross(lever1, -dir)  ) );

        m_nbConstraints++;

    }
    constraintId +=m_nbConstraints;
    dmsg_info() << "end buildConstraintMatrix " ;
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::getConstraintViolation(const ConstraintParams*, BaseVector *v,
                                                                     const DataVecCoord &, const DataVecDeriv &)
{
    unsigned int nb = m_violations.size();
    for(unsigned int i=0; i<nb; i++)
    {
        v->set(m_cid+i, m_violations[i]);
    }
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::getConstraintResolution(const ConstraintParams*,std::vector<ConstraintResolution*>& resTab,
                                                                      unsigned int& offset)
{
    unsigned int nb = m_violations.size();
    for(unsigned int i=0; i<nb; i++)
    {
        resTab[offset] = new AdaptiveBeamLengthConstraintResolution(nullptr, &m_constraintIntervals[i].active);
        offset++;
    }
}

template<class DataTypes>
void AdaptiveBeamLengthConstraint<DataTypes>::draw(const VisualParams* vparams)
{
    if (!vparams->displayFlags().getShowInteractionForceFields())
        return;

    if(m_constraintIntervals.size()==0)
        return;

    vparams->drawTool()->saveLastState();
    vparams->drawTool()->setLightingEnabled(false);

    std::vector< sofa::type::Vec3 > points;
    std::vector< sofa::type::RGBAColor> colors;
    points.reserve(m_constraintIntervals.size()*2);
    colors.reserve(m_constraintIntervals.size());

    for (unsigned int i = 0; i < m_constraintIntervals.size(); i++)
    {
        points.emplace_back(m_constraintIntervals[i].posBegin);
        points.emplace_back(m_constraintIntervals[i].posEnd);

        if (m_constraintIntervals[i].active)
            colors.emplace_back(1.0f, 0.5f, 0.0f, 1.0f);
        else
            colors.emplace_back(0.0f, 1.0f, 0.0f, 1.0f);
    }
    vparams->drawTool()->drawPoints(points, 10, sofa::type::RGBAColor::green());
    vparams->drawTool()->drawLines(points, 1.0f, colors);

    vparams->drawTool()->restoreLastState();
}

} // namespace sofa::component::constraintset::_adaptivebeamlengthconstraint_
