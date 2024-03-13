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


#include <BeamAdapter/component/constraint/ImplicitSurfaceAdaptiveConstraint.h>

#include <sofa/helper/system/thread/CTime.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/DataTypeInfo.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/core/ConstraintParams.h>
#include <sofa/component/constraint/lagrangian/model/UnilateralInteractionConstraint.h>


double TimeResolution = 0.0;
double TimeCount = 0.0;
double TimeProjection = 0.0;
double TimeProjection2= 0.0;

//eulalie: can I remove these?
//#define DEBUG
//#define DEBUG_TOPO_DATA
//#define DEBUG_PROJECTION
//#define DEBUG_LAST_CONSTRAINT_ONLY
//#define DEBUG_DFREE_COMPUTATION


namespace sofa::component::constraint
{

using sofa::core::VecCoordId;
using sofa::core::MultiVecCoordId;
using sofa::core::ConstMultiVecCoordId;
using sofa::core::ConstVecCoordId;

using sofa::component::constraint::lagrangian::model::UnilateralConstraintResolution;
using sofa::component::constraint::lagrangian::model::UnilateralConstraintResolutionWithFriction;

template<class DataTypes>
ImplicitSurfaceAdaptiveConstraint<DataTypes>::ImplicitSurfaceAdaptiveConstraint(MechanicalState* object1, MechanicalState* object2)
    : l_wireBinterpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , d_visualization(initData(&d_visualization, false, "visualization", "visualization of the implicit surface potential"))
    , d_posMin(initData(&d_posMin, Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , d_posMax(initData(&d_posMax, Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , d_initDomain(initData(&d_initDomain, 0, "domainId", "optional1: give an initial id for the domain"))
    , d_alarmDistance(initData(&d_alarmDistance, (Real) 1.0, "alarmDistance", "kind of alarm distance computed on implicit surface"))
    , d_radius(initData(&d_radius, (Real)1.0, "radius", "radius: this value is used to include the thickness in the collision response"))
    , d_frictionCoef(initData(&d_frictionCoef, (Real)0.0, "frictionCoef", "coefficient of friction (Coulomb's law)"))
    , d_listBeams(initData(&d_listBeams, "listBeams", "list of beams used by Interpolation that are activated for contact (all by default)"))
{
    SOFA_UNUSED(object1);
    SOFA_UNUSED(object2);
    m_isoValue=0.0;
}

template<class DataTypes>
ImplicitSurfaceAdaptiveConstraint<DataTypes>::ImplicitSurfaceAdaptiveConstraint(MechanicalState* object)
    :  l_wireBinterpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , d_visualization(initData(&d_visualization, false, "visualization", "visualization of the implicit surface potential"))
    , d_posMin(initData(&d_posMin, Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , d_posMax(initData(&d_posMax, Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , d_initDomain(initData(&d_initDomain, 0, "domainId", "optional2: give an initial id for the domain"))
    , d_alarmDistance(initData(&d_alarmDistance, (Real) 0.0, "alarmDistance", "the point for which the potential is more than the threshold are desactivated"))
    , d_radius(initData(&d_radius, (Real)1.0, "radius", "radius: this value is used to include the thickness in the collision response"))
    , d_frictionCoef(initData(&d_frictionCoef, (Real)0.0, "frictionCoef", "coefficient of friction (Coulomb's law)"))
    , d_listBeams(initData(&d_listBeams, "listBeams", "list of beams used by Interpolation that are activated for contact (all by default)"))
{
    SOFA_UNUSED(object);
    m_isoValue=0.0;
}

template<class DataTypes>
ImplicitSurfaceAdaptiveConstraint<DataTypes>::ImplicitSurfaceAdaptiveConstraint()
    :  l_wireBinterpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , d_visualization(initData(&d_visualization, false, "visualization", "visualization of the implicit surface potential"))
    , d_posMin(initData(&d_posMin, Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , d_posMax(initData(&d_posMax, Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , d_initDomain(initData(&d_initDomain, 0, "domainId", "optional3: give an initial id for the domain"))
    , d_alarmDistance(initData(&d_alarmDistance, (Real) 0.0, "alarmDistance", "the point for which the potential is more than the threshold are desactivated"))
    , d_radius(initData(&d_radius, (Real)1.0, "radius", "radius: this value is used to include the thickness in the collision response"))
    , d_frictionCoef(initData(&d_frictionCoef, (Real)0.0, "frictionCoef", "coefficient of friction (Coulomb's law)"))
    , d_listBeams(initData(&d_listBeams, "listBeams", "list of beams used by Interpolation that are activated for contact (all by default)"))
{
    m_isoValue=0.0;
}


template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::init()
{
    Inherit::init();
}


template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::clear()
{
    internalInit();
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::reset()
{
    internalInit();
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::bwdInit()
{
    internalInit();
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::internalInit()
{
    // STEP 1 : get the implicit surface
    if(mstate1)
        mstate1->getContext()->get(m_contactSurface);
    else
        msg_error() <<"No mstate1 found.";

    if (m_contactSurface == nullptr)
    {
        msg_error() <<"No surface found for contact.";
        return;
    }
    else
    {
        msg_info() <<"Contact with surface name:"<<m_contactSurface->getName();
    }

    // STEP 2: verify that we have the mechanical state of the Adaptive Beams
    if (!this->mstate2)
        msg_warning() <<"No mstate2 found";

    // STEP 3: given the coef of Friction, set the bool friction [friction or frictionless contact]
    if(d_frictionCoef.getValue()>0.0)
        m_friction=true;
    else
        m_friction=false;


    // STEP 4: verify if listBeams is void => all beams are activated
    const type::vector<int> &list_B= d_listBeams.getValue();
    if( list_B.size() == 0)
        m_allActivated=true;
    else
        m_allActivated=false;

    // STEP 5: init domain
    m_domainSample.clear();
    m_domainSample.push_back(-1);

    // STEP 6: init marching cube
#ifdef SOFAEVE
    mc = new sofaeve::implicit::MarchingCube();
#endif


}


template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getOrthogonalVectors(const Vec3& dir, Vec3& vec1, Vec3& vec2)
{
    Vec3 temp; // Any vector such as temp != dir
    temp[0] = dir[1];
    temp[1] = dir[2];
    temp[2] = dir[0];

    if(temp == dir) // x = y = z
        temp = Vec3(1,0,0);

    vec1 = cross(dir, temp);
    vec1.normalize();

    vec2 = cross(dir, vec1);
    vec2.normalize();
}


/// function detectPotentialContactOnImplicitSurface()
/// computation of the distance with the implicit surface  (and verify that it is < alarmDistance)
/// "filters" contacts (the normal of contact is supposed to be perp to the tangent of the curve except on the tip point)
/// set  normal and tangential directions
template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::detectPotentialContactOnImplicitSurface(const ConstVecCoordId &vecXId, type::vector<int>& listBeam)
{
    unsigned int numBeams= l_wireBinterpolation->getNumBeams();

    VectorIntIterator it;
    while( m_domainSample.size() < m_posSample.size() )
    {
        it=m_domainSample.begin();
        m_domainSample.insert(it,-1);
    }

    while ( m_domainSample.size() > m_posSample.size() )
    {
        it=m_domainSample.begin();
        m_domainSample.erase(it);
    }

#ifdef DEBUG
    if(m_domainSample.size()!= m_posSample.size() )
    {
        msg_error()<<"Domain and pos samples do not have the same size : cancel detection";
        return;
    }
#endif
    for (unsigned int p=0; p<m_posSample.size(); p++)
    {
        Vec3d pos=(Vec3d)m_posSample[p];
        m_domainSample[p] = m_contactSurface->getDomain(pos, m_domainSample[p]);
        Real value = m_contactSurface->getValue(pos, m_domainSample[p]);
        Vec3 grad = (Vec3)m_contactSurface->getGradient(pos, m_domainSample[p]);

        Real gradNorm= grad.norm();
        if (gradNorm > 1e-12 )
        {
            Real d=value/gradNorm;

            if (d<this->d_alarmDistance.getValue())
            {
                potentialContact pt;
                unsigned int bi = p / 10;
                pt.beamId = listBeam[bi];
                pt.posSampleId = p;
                Real bc = (p+1-(10*bi))/10.0;
                pt.baryCoord=Vec3(bc,0,0); // todo: change to account for the radius
                grad.normalize();
                pt.n=grad;
                pt.d = d;

                // get the tangent to the curve at this point...
                Vec3 t;
                l_wireBinterpolation->getTangentUsingSplinePoints( pt.beamId, bc, vecXId, t );
                t.normalize();
                pt.t = t;

                // for all contact (except the one on the tip of the last beam), the normal is modified to be perp to the tangent of the curve
                double dnt=dot(pt.n, pt.t);
                if(pt.beamId ==(numBeams-1) && pt.baryCoord[0]>0.999)
                {
                    dmsg_info() <<" Last contact point (no modif on the normal) id:"<<pt.beamId <<"  bc"<<  pt.baryCoord;

                    // build a "frame" for friction contact:
                    Vec3 newT;
                    if (fabs(dnt)<0.999)
                    {
                        newT = pt.t - pt.n*dnt;
                        newT.normalize();
                        pt.t=newT;
                    }
                    else
                    {   // in case pt.t and pt.n are aligned:
                        getOrthogonalVectors(pt.n, pt.t, pt.s);
                    }
                }
                else
                {
                    if(fabs(dnt) < 0.5)
                    {
                        Vec3 newN = pt.n - pt.t*dnt;
                        newN.normalize();
                        pt.n = newN;
                    }
                    else
                    {
                        continue;
                    }

                }
                pt.s = cross(pt.n, pt.t);

                m_vecPotentialContact.push_back(pt);
            }
        }
    }
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams* cParams,
                                                                         DataMatrixDeriv &c1,
                                                                         DataMatrixDeriv &c2,
                                                                         unsigned int &cIndex,
                                                                         const DataVecCoord &x1,
                                                                         const DataVecCoord &x2)
{
    SOFA_UNUSED(c1);
    SOFA_UNUSED(x1);

#ifdef DEBUG
    dmsg_info() <<"Entering buildConstraintMatrix ";
#endif

    //1st step: computation of sample point position (10 / beam)
    // if all activated => all beam considered otherwize, use listBeams...
    type::vector<int> list_B;
    unsigned int numBeams= l_wireBinterpolation->getNumBeams();
    m_allActivated=false;
    if (m_allActivated)
    {
        list_B.clear();
        for ( unsigned int i=0; i<numBeams; i++)
            list_B.push_back(i);
    }
    else
    {
        list_B= d_listBeams.getValue();
    }
    numBeams = list_B.size();
    m_posSample.clear();
    m_vecPotentialContact.clear();

    dmsg_info() <<" list_B = "<<list_B;

    if (list_B.size()==0)
        return;

    m_posSample.resize(10*numBeams);

    //PosSample corresponds to points that are placed along the spline which position is stored in sofa::core::VecCoordId::position()
    // the following code verify that the position is correct:
    dmsg_info() <<" in buildConstraintMatrix: cParams->x()="<<cParams->x();

    MultiVecCoordId x = VecCoordId::position(); 
    ConstVecCoordId xtest = cParams->x()[this->mstate2.get()];
    //ConstVecCoordId xtest = cParams->readX(this->mstate2);

    if(xtest != x.getId(this->mstate2))
    {
        msg_warning() <<"In buildConstraintMatrix, cParams->x() != sofa::core::VecCoordId::position()";
    }

    // update the position of the Bezier Points (if necessary)
    const VecCoord &x2buf = x2.getValue();
    sofa::core::VecCoordId x_in = sofa::core::VecCoordId::position();

    auto x2test = sofa::helper::getReadAccessor(*(this->mstate2->read(core::ConstVecCoordId::position())));
    l_wireBinterpolation->updateBezierPoints(x2test, x_in );

    // interpolates the position of the sample points
    for (unsigned int bi=0; bi<numBeams; bi++)
    {
        unsigned int b= list_B[bi];
#ifdef DEBUG
        dmsg_info() <<" interpolate point on beam b="<<b<<" bary Coord.";
#endif
        for (unsigned int p=0;p<10; p++)
        {
            Vec3 localPos(0,0,0);
            Real baryCoord = (p+1)/10.0;
#ifdef DEBUG
        dmsg_info() <<" ("<<baryCoord<<") ";
#endif

            l_wireBinterpolation->interpolatePointUsingSpline(b, baryCoord, localPos, x2buf, m_posSample[10*bi+p], false, x_in);
        }
#ifdef DEBUG
        dmsg_info() <<" ";
#endif
    }
#ifdef DEBUG
    dmsg_info() <<" * 1st step ok: m_posSample = "<<m_posSample;
#endif
    /////////////////
    // 2d step: computation of the distance with the implicit surface  (and verify that it is < alarmDistance)
    this->detectPotentialContactOnImplicitSurface(x_in, list_B);

#ifdef DEBUG
    dmsg_info() <<" * 2d step ok:  "<<m_VecPotentialContact.size()<<" potential(s) contact detected";
#endif


    /////////////////
    // 3d step: setting the constraint direction
    auto c2w = sofa::helper::getWriteOnlyAccessor(c2);

    m_nbConstraints = 0;
    m_cid = cIndex;
    for (unsigned int i=0; i<m_vecPotentialContact.size(); i++)
    {

        potentialContact &pt= m_vecPotentialContact[i];



        // computes the mapping of the force on the DOFs' space => N_node0 and N_node1 (6dofs each)
        SpatialVector N_node0, N_node1;
        l_wireBinterpolation->MapForceOnNodeUsingSpline(pt.beamId, pt.baryCoord[0], Vec3(0.0, pt.baryCoord[1],pt.baryCoord[2]),
                x2buf, pt.n, N_node0, N_node1);

        unsigned int node0Idx, node1Idx;
        l_wireBinterpolation->getNodeIndices( pt.beamId,  node0Idx, node1Idx );


        // put the normal direction of contact in the Matrix of constraint directions

        MatrixDerivRowIterator c2_it = c2w.writeLine(m_cid + m_nbConstraints);

        c2_it.addCol(node0Idx, Deriv(N_node0.getForce(), N_node0.getTorque() ) );
        c2_it.addCol(node1Idx, Deriv(N_node1.getForce(), N_node1.getTorque() ) );
        m_nbConstraints++;

        if (m_friction)
        {

            SpatialVector T_node0, T_node1, S_node0, S_node1 ;

            // put the first tangential direction of contact in the Matrix of constraint directions
            l_wireBinterpolation->MapForceOnNodeUsingSpline(pt.beamId, pt.baryCoord[0], Vec3(0.0, pt.baryCoord[1],pt.baryCoord[2]),
                    x2buf, pt.t, T_node0, T_node1);

            c2_it = c2w.writeLine(m_cid + m_nbConstraints);
            c2_it.addCol(node0Idx, Deriv(T_node0.getForce(), T_node0.getTorque() ) );
            c2_it.addCol(node1Idx, Deriv(T_node1.getForce(), T_node1.getTorque() ) );
            m_nbConstraints++;


            // put the second tangential direction of contact in the Matrix of constraint directions
            l_wireBinterpolation->MapForceOnNodeUsingSpline(pt.beamId, pt.baryCoord[0], Vec3(0.0, pt.baryCoord[1],pt.baryCoord[2]),
                    x2buf, pt.s, S_node0, S_node1);

            c2_it = c2w.writeLine(m_cid + m_nbConstraints);
            c2_it.addCol(node0Idx, Deriv(S_node0.getForce(), S_node0.getTorque() ) );
            c2_it.addCol(node1Idx, Deriv(S_node1.getForce(), S_node1.getTorque() ) );
            m_nbConstraints++;
        }

    }

#ifdef DEBUG
    dmsg_info() <<" * 3d step ok: constraints are set...\n done *************";
#endif

    cIndex+=m_nbConstraints;


}



template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::computeTangentialViolation(const Vec3 &Pos, const Vec3 &freePos, const Vec3 &t, const Vec3 &s,
                                                                              const Real& d, const Real& dfree, Real &dfree_t, Real &dfree_s )
{

    //  Real dfree = pt.d + dot(pt.n, freePos - m_posSample[pt.posSampleId]) - radius.getValue();

    Real dTest= d-d_radius.getValue();

    if ( fabs(dTest) <=0.01*d_alarmDistance.getValue())  // actual position is already in contact
    {
#ifdef DEBUG_DFREE_COMPUTATION
        std::cout<<" case 1"<<std::endl;
#endif
        dfree_t = dot(t, freePos-Pos);
        dfree_s = dot(s, freePos-Pos);
    }
    else if(fabs(dTest-dfree) > 0.001*dTest) // significant variation of d between actual position and free position
    {

        if (dfree>0)
        {
            // not in contact yet...
            dfree_t =0.0;
            dfree_s =0.0;
#ifdef DEBUG_DFREE_COMPUTATION
            std::cout<<" case 2"<<std::endl;
#endif
        }
        else if (dTest<=0.0)
        {
            // already in contact
            dfree_t = dot(t, freePos-Pos);
            dfree_s = dot(s, freePos-Pos);
#ifdef DEBUG_DFREE_COMPUTATION
            std::cout<<" case 3"<<std::endl;
#endif
        }
        else
        {
            // dfree < 0 and d > 0 => Try to find the "time" and position of collision
            Real dt = dTest/(dTest-dfree);
            dt=std::max((Real)0.0,dt);
            dt=std::min((Real)1.0,dt);

            Vec3 Pt= Pos * (1-dt) + freePos*dt;

            dfree_t = dot(t, freePos-Pt);
            dfree_s = dot(s, freePos-Pt);
#ifdef DEBUG_DFREE_COMPUTATION
            std::cout<<" case 4"<<std::endl;
#endif
        }

    }
    else
    {
        // no significant variation of d between actual position and free position
        dfree_t =0.0;
        dfree_s =0.0;
#ifdef DEBUG_DFREE_COMPUTATION
        std::cout<<" case 5"<<std::endl;
#endif
    }

}


template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getConstraintViolation(const core::ConstraintParams* cParams,
                                                                          linearalgebra::BaseVector *v,
                                                                          const DataVecCoord &x1,
                                                                          const DataVecCoord &x2,
                                                                          const DataVecDeriv &v1,
                                                                          const DataVecDeriv &v2)
{
    SOFA_UNUSED(cParams);
    SOFA_UNUSED(x1);
    SOFA_UNUSED(v1);
    SOFA_UNUSED(v2);

#ifdef DEBUG
    std::cout<<" entering getConstraintViolation: cParams->x()="<<cParams->x()<<std::endl;
#endif

    // update the position of the Bezier Points (if necessary)
    const VecCoord &x2buf = x2.getValue();
    sofa::core::VecCoordId x_in = sofa::core::VecCoordId::freePosition();
    l_wireBinterpolation->updateBezierPoints(x2buf, x_in );
    dmsg_info()<<" updateBezierPoints ok";

    dmsg_info()<<" m_VecPotentialContact size = "<<m_vecPotentialContact.size();
    // computes the free position of the potential contact points
    unsigned int itConstraints=0;
    for (unsigned int i=0; i<m_vecPotentialContact.size(); i++)
    {
        potentialContact &pt= m_vecPotentialContact[i];
        Vec3 freePos;
        l_wireBinterpolation->interpolatePointUsingSpline(pt.beamId, pt.baryCoord[0], Vec3(0,0,0), x2buf, freePos, false, x_in);

        Real dfree = pt.d + dot(pt.n, freePos - m_posSample[pt.posSampleId]) - d_radius.getValue();
        v->set(m_cid+itConstraints,dfree);
        itConstraints++;

        if(m_friction)
        {
            Real dfree_t=0.0;
            Real dfree_s=0.0;

            computeTangentialViolation(m_posSample[pt.posSampleId], freePos, pt.t, pt.s, pt.d, dfree, dfree_t, dfree_s);

            // along direction t
            v->set(m_cid+itConstraints, dfree_t);
            itConstraints++;

            // along direction s
            v->set(m_cid+itConstraints, dfree_s);
            itConstraints++;

        }


    }
#ifdef DEBUG
    std::cout<<" leaving getConstraintViolation"<<std::endl;
#endif
}


template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution* > & resTab, unsigned int& offset)
{

#ifdef DEBUG
    std::cout<<" entering getConstraintResolution"<<std::endl;
#endif

    for (unsigned int i=0; i<m_vecPotentialContact.size(); i++)
    {
        if(m_friction)
        {

            resTab[offset] = new UnilateralConstraintResolutionWithFriction(d_frictionCoef.getValue());
            offset+=3;
        }
        else
        {
            resTab[offset] = new UnilateralConstraintResolution();
            offset +=1;
        }
    }
#ifdef DEBUG
    dmsg_info()<<" leaving getConstraintResolution";
#endif
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifdef DEBUG
    dmsg_info()<<" entering draw";
#endif

    if (!vparams->displayFlags().getShowInteractionForceFields()) return;
    TimeCount =0;
    TimeResolution =0;
    TimeProjection =0;
    TimeProjection2=0;


#ifdef DEBUG_LAST_CONSTRAINT_ONLY
    const VecCoord& x		= (*this->mstate2->getX());
    unsigned int _nbNodes = x.size();
#endif

    const auto stateLifeCycle = vparams->drawTool()->makeStateLifeCycle();
    vparams->drawTool()->disableLighting();

    std::vector<sofa::type::Vec3> pointsToDraw;
    std::vector<type::RGBAColor> colors;

    typename sofa::type::vector<potentialContact>::iterator it = m_vecPotentialContact.begin();
    for (unsigned int p=0; p<m_posSample.size(); p++)
    {
        if(m_vecPotentialContact.size()>0 && (*it).posSampleId==p) // in potiential contact
        {
            colors.push_back(type::RGBAColor::red());
            pointsToDraw.push_back(m_posSample[p]);

            if (it!=m_vecPotentialContact.end())
                it++;
        }
        else
        {
            colors.push_back(type::RGBAColor::cyan());
            pointsToDraw.push_back(m_posSample[p]);
        }

    }

    vparams->drawTool()->drawPoints(pointsToDraw, 2, colors);

    colors.clear();
    std::vector<sofa::type::Vec3> linesToDraw;

    for (unsigned int i=0; i<m_vecPotentialContact.size(); i++)
    {
        potentialContact &pt= m_vecPotentialContact[i];
        Vec3 Pos = m_posSample[pt.posSampleId];

        colors.push_back(type::RGBAColor::red());
        linesToDraw.push_back(Pos);
        linesToDraw.push_back(Pos + pt.n);

        if(m_friction)
        {
            colors.push_back(type::RGBAColor::green());
            linesToDraw.push_back(Pos);
            linesToDraw.push_back(Pos + pt.t);

            colors.push_back(type::RGBAColor::blue());
            linesToDraw.push_back(Pos);
            linesToDraw.push_back(Pos + pt.s);
        }

    }
    vparams->drawTool()->drawLines(linesToDraw, 1, colors);

    vparams->drawTool()->enableLighting();

    if (d_visualization.getValue() && m_posSample.size()>0)
    {

        dmsg_info()<<" Visualization ";
        dmsg_info()<<" m_domainSample"<<m_domainSample;

        double dx = (d_posMax.getValue()[0] -d_posMin.getValue()[0])/50;
        double dy = (d_posMax.getValue()[1] -d_posMin.getValue()[1])/50;
        double dz = (d_posMax.getValue()[2] -d_posMin.getValue()[2])/50;
        double z = d_posMin.getValue()[2];
        int idx =0;

        double mmin = 150;
        double mmax = -150;

        type::Vec3d PosLastPoint=m_posSample[m_posSample.size()-1];

        while(z<=d_posMax.getValue()[2])
        {
            double y = d_posMin.getValue()[1];
            while(y<=d_posMax.getValue()[1])
            {
                double x = d_posMin.getValue()[0];

                while(x<=d_posMax.getValue()[0])
                {
                    Vec3d Pos(x,y,z);
                    Pos += PosLastPoint;
                    //Pos += PP;
                    double test = m_contactSurface->getValue(Pos, m_domainSample[m_posSample.size()-1]);
                    if(test > mmax)
                        mmax = test;
                    if(test < mmin)
                        mmin = test;
                    m_cData[idx] = test;
                    idx++;

                    x = x+dx;
                }
                y = y+dy;

            }
            z = z+dz;
        }

        vparams->drawTool()->setPolygonMode(1, true);

#ifdef SOFAEVE
        mc->buildMesh(mc_data, 50, 50, 50, _isoValue);
#endif

        vparams->drawTool()->setPolygonMode(0, false);


        Vec3d Pos(d_posMin.getValue()[0],d_posMin.getValue()[1],d_posMin.getValue()[2]);
        Pos += PosLastPoint;
#ifdef SOFAEVE
        mc->draw(Pos[0], Pos[1], Pos[2], 0.2, 0.2, 0.2);
#endif

    }

#ifdef DEBUG
    dmsg_info() <<" leaving draw";
#endif

}


////////////////// IMPLICITSURFACEADAPTIVECONSTRAINTRESOLUTION //////////////////////

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraintResolution<DataTypes>::resolution(int line, double** w, double* d, double* force)
{
    // dummy resolution
    d[line]=w[line][line]*force[line];
    force[line]=0.0;
}

} // namespace sofa::component::constraint
