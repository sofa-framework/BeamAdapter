#ifndef SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_INL
#define SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_INL


#include "ImplicitSurfaceAdaptiveConstraint.h"

#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/gl/template.h>
#include <sofa/helper/vector.h>

#include <sofa/helper/io/MassSpringLoader.h>
#include <sofa/helper/gl/template.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/DataTypeInfo.h>
#include <BaseTopology/RegularGridTopology.h>
#include <BaseMechanics/AddMToMatrixFunctor.h>
#include <BaseTopology/TopologyData.inl>

#include <sofa/helper/system/thread/CTime.h>
#include <Constraint/constraintset/UnilateralInteractionConstraint.h>

#include <sofa/core/visual/VisualParams.h>


double TimeResolution = 0.0;
double TimeCount = 0.0;
double TimeProjection = 0.0;
double TimeProjection2= 0.0;
//#define DEBUG
//#define DEBUG_TOPO_DATA
//#define DEBUG_PROJECTION
//#define DEBUG_LAST_CONSTRAINT_ONLY


//#define DEBUG_DFREE_COMPUTATION


namespace sofa
{

namespace component
{

namespace constraint
{



//////////////////////////////
// ImplicitSurfaceAdaptiveConstraint
//////////////////////////////



template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::init()
{
std::cout<<"ImplicitSurfaceAdaptiveConstraint<DataTypes>::init()"<<std::endl;
    assert(this->mstate1);
    assert(this->mstate2);
}




template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::clear()
{
    std::cout<<"ImplicitSurfaceAdaptiveConstraint<DataTypes>::clear()"<<std::endl;
    internalInit();
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::reset()
{
    std::cout<<"ImplicitSurfaceAdaptiveConstraint<DataTypes>::reset()"<<std::endl;
    internalInit();
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::bwdInit()
{
    std::cout<<"ImplicitSurfaceAdaptiveConstraint<DataTypes>::bwdInit()"<<std::endl;
    internalInit();
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::internalInit()
{
    std::cout<<"ImplicitSurfaceAdaptiveConstraint<DataTypes>::internalInit()"<<std::endl;

    //  STEP 1 : get the implicit surface
    if (this->mstate1)
        this->mstate1->getContext()->get(_contact_surface);
    else
        serr<<" no mstate1 found: WARNING !!"<<sendl;

    if (_contact_surface == NULL){
        serr<<"oooooooooo\n oooERROR: no surface found for contact"<<sendl;
        return;
    }
    else
    {
        std::cout<<"CONTACT with surface name:"<<_contact_surface->getName()<<std::endl;
    }


    // STEP 2: verify that we have the mechanical state of the Adaptive Beams

    if (!this->mstate2)
        serr<<" no mstate2 found: WARNING !!"<<sendl;

    // STEP 3: given the coef of Friction, set the bool friction [friction or frictionless contact]
    if(frictionCoef.getValue()>0.0)
        friction=true;
    else
        friction=false;


    // STEP 4: verify if listBeams is void => all beams are activated
    const sofa::helper::vector<int> &list_B= listBeams.getValue();
    if( list_B.size() == 0)
        all_activated=true;
    else
        all_activated=false;

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
    Vec3 temp;	// Any vector such as temp != dir
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




/*

template <class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::handleTopologyChange()
{

    #ifdef DEBUG_TOPO_DATA
        std::cout<<"handleTopologyChange before : ";

        int m2 = this->mstate2->getSize();
        for (int i=0; i<m2-1; i++)
        {
            std::cout<<"listProbe["<<i<<"].pos :"<<listProbe[i].pos()<<std::endl;
            std::cout<<"                          ";
        }
        std::cout<<""<<std::endl;
    #endif

    core::topology::BaseMeshTopology* topology = this->mstate2->getContext()->getMeshTopology();
    std::list<const core::topology::TopologyChange *>::const_iterator itBegin=topology->beginChange();
    std::list<const core::topology::TopologyChange *>::const_iterator itEnd=topology->endChange();
    std::list<const core::topology::TopologyChange *>::const_iterator it;
    listProbe.handleTopologyEvents(itBegin,itEnd);

    #ifdef DEBUG_TOPO_DATA
         m2 = this->mstate2->getSize();
        for (int i=0; i<m2; i++)
        {
            std::cout<<"                     after  : listProbe["<<i<<"].pos :"<<listProbe[i].pos()<<std::endl;
        }
        std::cout<<""<<std::endl;


    #endif

}

*/



/// function detectPotentialContactOnImplicitSurface()
/// computation of the distance with the implicit surface  (and verify that it is < alarmDistance)
/// "filters" contacts (the normal of contact is supposed to be perp to the tangent of the curve except on the tip point)
/// set  normal and tangential directions
template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::detectPotentialContactOnImplicitSurface(const sofa::core::ConstVecCoordId &vecXId, sofa::helper::vector<int>& listBeam)
{

    unsigned int numBeams= m_wireBinterpolation->getNumBeams();


    vectorIntIterator it;
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
        serr<<" domain and pos samples do not have the same size : cancel detection"<<send:
        return;
    }
#endif
    for (unsigned int p=0; p<m_posSample.size(); p++)
    {
        Vec3d pos=(Vec3d)m_posSample[p];
        m_domainSample[p] = _contact_surface->getDomain(pos, m_domainSample[p]);
        Real value = _contact_surface->getValue(pos, m_domainSample[p]);
        Vec3 grad = (Vec3)_contact_surface->getGradient(pos, m_domainSample[p]);

        Real gradNorm= grad.norm();
        if (gradNorm > 1e-12 )
        {
            Real d=value/gradNorm;

            if (d<this->alarmDistance.getValue())
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
                m_wireBinterpolation->getTangentUsingSplinePoints( pt.beamId, bc, vecXId, t );
                t.normalize();
                pt.t = t;


                ///////// force n perp t:
                // for all contact (except the one on the tip of the last beam), the normal is modified to be perp to the tangent of the curve
                double dnt=dot(pt.n, pt.t);
                if(pt.beamId ==(numBeams-1) && pt.baryCoord[0]>0.999)
                {
                     std::cout<<" ***********\n last contact point (no modif on the normal) id:"<<pt.beamId <<"  bc"<<  pt.baryCoord<<std::endl;


                     // build a "frame" for friction contact:
                     Vec3 newT;
                     if (fabs(dnt)<0.999)
                     {
                           newT = pt.t - pt.n*dnt;
                           newT.normalize();
                           pt.t=newT;
                     }
                     else
                     {      // in case pt.t and pt.n are aligned:
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
                       // std::cout<<" warning n and t are far from being perp => contact cancelled"<<std::endl;
                        continue;
                    }

                }
                pt.s = cross(pt.n, pt.t);

                m_VecPotentialContact.push_back(pt);

            }
        }

    }

}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::buildConstraintMatrix(const sofa::core::ConstraintParams* cParams , DataMatrixDeriv &/*c1*/, DataMatrixDeriv &c2, unsigned int &cIndex
                                                                         , const DataVecCoord &/*x1*/, const DataVecCoord &x2)
{

#ifdef DEBUG
    std::cout<<"Entering buildConstraintMatrix "<<std::endl;
#endif

    /////////////////
    //1st step: computation of sample point position (10 / beam)

    // if all activated => all beam considered otherwize, use listBeams...
    sofa::helper::vector<int> list_B;
    unsigned int numBeams= m_wireBinterpolation->getNumBeams();
    all_activated=false;
    if (all_activated)
    {
        list_B.clear();
        for ( unsigned int i=0; i<numBeams; i++)
            list_B.push_back(i);

    }
    else
    {
        list_B= listBeams.getValue();
    }
    numBeams = list_B.size();
    m_posSample.clear();
    m_VecPotentialContact.clear();

    std::cout<<" list_B = "<<list_B<<std::endl;

    if (list_B.size()==0)
        return;

    m_posSample.resize(10*numBeams);

    //PosSample corresponds to points that are placed along the spline which position is stored in sofa::core::VecCoordId::position()
    // the following code verify that the position is correct:
    std::cout<<" in buildConstraintMatrix: cParams->x()="<<cParams->x()<<std::endl;

    sofa::core::MultiVecCoordId x = sofa::core::VecCoordId::position();
    const sofa::core::ConstMultiVecCoordId &xId = cParams->x();
    sofa::core::ConstVecCoordId xtest = xId.getId(this->mstate2);

    if(  xtest != x.getId(this->mstate2))
    {
        serr<<" WARNING in buildConstraintMatrix, cParams->x() != sofa::core::VecCoordId::position()"<<sendl;
    }

    // update the position of the Bezier Points (if necessary)
    const VecCoord &x2buf = x2.getValue();
    sofa::core::VecCoordId x_in = sofa::core::VecCoordId::position();

    const VecCoord &x2test=*this->mstate2->getX();
    m_wireBinterpolation->updateBezierPoints(x2test, x_in );

    // interpolates the position of the sample points
    for (unsigned int bi=0; bi<numBeams; bi++)
    {
        unsigned int b= list_B[bi];
#ifdef DEBUG
        std::cout<<" interpolate point on beam b="<<b<<" bary Coord:";
#endif
        for (unsigned int p=0;p<10; p++)
        {
            Vec3 localPos(0,0,0);
            Real baryCoord = (p+1)/10.0;
#ifdef DEBUG
            std::cout<<" ("<<baryCoord<<") ";
#endif

            m_wireBinterpolation->interpolatePointUsingSpline(b, baryCoord, localPos, x2buf, m_posSample[10*bi+p], false, x_in);
        }
 #ifdef DEBUG
        std::cout<<" "<<std::endl;
#endif
    }
#ifdef DEBUG
    std::cout<<" * 1st step ok: m_posSample = "<<m_posSample<<std::endl;
#endif
    /////////////////
    // 2d step: computation of the distance with the implicit surface  (and verify that it is < alarmDistance)
    this->detectPotentialContactOnImplicitSurface(x_in, list_B);

#ifdef DEBUG
    std::cout<<" * 2d step ok:  "<<m_VecPotentialContact.size()<<" potential(s) contact detected"<<std::endl;
#endif


    /////////////////
    // 3d step: setting the constraint direction
    MatrixDeriv& c2w = *c2.beginEdit();
    _nbConstraints = 0;
    cid = cIndex;
    for (unsigned int i=0; i<m_VecPotentialContact.size(); i++)
    {

        potentialContact &pt= m_VecPotentialContact[i];



        // computes the mapping of the force on the DOFs' space => N_node0 and N_node1 (6dofs each)
        SpatialVector N_node0, N_node1;
        m_wireBinterpolation->MapForceOnNodeUsingSpline(pt.beamId, pt.baryCoord[0], Vec3(0.0, pt.baryCoord[1],pt.baryCoord[2]),
                                                        x2buf, pt.n, N_node0, N_node1);

        unsigned int node0Idx, node1Idx;
        m_wireBinterpolation->getNodeIndices( pt.beamId,  node0Idx, node1Idx );


        // put the normal direction of contact in the Matrix of constraint directions

        MatrixDerivRowIterator c2_it = c2w.writeLine(cid + _nbConstraints);

        c2_it.addCol(node0Idx, Deriv(N_node0.getForce(), N_node0.getTorque() ) );
        c2_it.addCol(node1Idx, Deriv(N_node1.getForce(), N_node1.getTorque() ) );
        _nbConstraints++;

        if (friction)
        {

            SpatialVector T_node0, T_node1, S_node0, S_node1 ;

            // put the first tangential direction of contact in the Matrix of constraint directions
            m_wireBinterpolation->MapForceOnNodeUsingSpline(pt.beamId, pt.baryCoord[0], Vec3(0.0, pt.baryCoord[1],pt.baryCoord[2]),
                                                            x2buf, pt.t, T_node0, T_node1);

            c2_it = c2w.writeLine(cid + _nbConstraints);
            c2_it.addCol(node0Idx, Deriv(T_node0.getForce(), T_node0.getTorque() ) );
            c2_it.addCol(node1Idx, Deriv(T_node1.getForce(), T_node1.getTorque() ) );
            _nbConstraints++;


            // put the second tangential direction of contact in the Matrix of constraint directions
            m_wireBinterpolation->MapForceOnNodeUsingSpline(pt.beamId, pt.baryCoord[0], Vec3(0.0, pt.baryCoord[1],pt.baryCoord[2]),
                                                            x2buf, pt.s, S_node0, S_node1);

            c2_it = c2w.writeLine(cid + _nbConstraints);
            c2_it.addCol(node0Idx, Deriv(S_node0.getForce(), S_node0.getTorque() ) );
            c2_it.addCol(node1Idx, Deriv(S_node1.getForce(), S_node1.getTorque() ) );
            _nbConstraints++;
        }

    }
    c2.endEdit();

#ifdef DEBUG
    std::cout<<" * 3d step ok: constraints are set...\n done *************"<<std::endl;
#endif

    cIndex+=_nbConstraints;


}





template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::computeTangentialViolation(const Vec3 &Pos, const Vec3 &freePos, const Vec3 &t, const Vec3 &s,
                                                                              const Real& d, const Real& dfree, Real &dfree_t, Real &dfree_s )
{

    //  Real dfree = pt.d + dot(pt.n, freePos - m_posSample[pt.posSampleId]) - radius.getValue();

    Real dTest= d-radius.getValue();

    if ( fabs(dTest) <=0.01*alarmDistance.getValue())  // actual position is already in contact
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
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getConstraintViolation(const sofa::core::ConstraintParams* /*cParams*/, defaulttype::BaseVector *v, const DataVecCoord &/*x1*/, const DataVecCoord &x2
                                                                               , const DataVecDeriv &/*v1*/, const DataVecDeriv &/*v2*/)
{
#ifdef DEBUG
    std::cout<<" entering getConstraintViolation: cParams->x()="<<cParams->x()<<std::endl;
#endif

    // update the position of the Bezier Points (if necessary)
    const VecCoord &x2buf = x2.getValue();
    sofa::core::VecCoordId x_in = sofa::core::VecCoordId::freePosition();
    m_wireBinterpolation->updateBezierPoints(x2buf, x_in );
    std::cout<<" updateBezierPoints ok"<<std::endl;

    std::cout<<" m_VecPotentialContact size = "<<m_VecPotentialContact.size()<<std::endl;
    // computes the free position of the potential contact points
    unsigned int itConstraints=0;
    for (unsigned int i=0; i<m_VecPotentialContact.size(); i++)
    {
        potentialContact &pt= m_VecPotentialContact[i];
        Vec3 freePos;
        m_wireBinterpolation->interpolatePointUsingSpline(pt.beamId, pt.baryCoord[0], Vec3(0,0,0), x2buf, freePos, false, x_in);

        Real dfree = pt.d + dot(pt.n, freePos - m_posSample[pt.posSampleId]) - radius.getValue();
         v->set(cid+itConstraints,dfree);
         itConstraints++;

         if(friction)
         {
             Real dfree_t=0.0;
             Real dfree_s=0.0;

             computeTangentialViolation(m_posSample[pt.posSampleId], freePos, pt.t, pt.s, pt.d, dfree, dfree_t, dfree_s);

             // along direction t
             v->set(cid+itConstraints, dfree_t);
             itConstraints++;

             // along direction s
             v->set(cid+itConstraints, dfree_s);
             itConstraints++;

         }


    }
 #ifdef DEBUG
    std::cout<<" leaving getConstraintViolation"<<std::endl;
#endif
}


template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{

#ifdef DEBUG
    std::cout<<" entering getConstraintResolution"<<std::endl;
#endif

    for (unsigned int i=0; i<m_VecPotentialContact.size(); i++)
    {
        if(friction)
        {

            resTab[offset] = new sofa::component::constraintset::UnilateralConstraintResolutionWithFriction(frictionCoef.getValue());
            offset+=3;
        }
        else
        {
            resTab[offset] = new sofa::component::constraintset::UnilateralConstraintResolution();
            offset +=1;
        }
    }
#ifdef DEBUG
    std::cout<<" leaving getConstraintResolution"<<std::endl;
#endif
}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{
#ifdef DEBUG
    std::cout<<" entering draw"<<std::endl;
#endif

    if (!vparams->displayFlags().getShowInteractionForceFields()) return;

    //double timeScale = 1.0 / (double)sofa::helper::system::thread::CTime::getRefTicksPerSec();
    //std::cout<<" time for TimeResolution : "<< TimeResolution *timeScale<<" s"<<std::endl;
    //std::cout<<" time for TimeCount      : "<< TimeCount      *timeScale<<" s"<<std::endl;
    //std::cout<<" time for TimeProjection : "<< TimeProjection *timeScale<<" s"<<std::endl;
    //std::cout<<" time for TimeProjection2: "<< TimeProjection2*timeScale<<" s"<<std::endl;
    TimeCount =0;
    TimeResolution =0;
    TimeProjection =0;
    TimeProjection2=0;


#ifdef DEBUG_LAST_CONSTRAINT_ONLY
    const VecCoord& x		= (*this->mstate2->getX());
    unsigned int _nbNodes = x.size();
    //std::cout << "Drawing centered on " << x << std::endl;
#endif

    glDisable(GL_LIGHTING);
    glPointSize(2);

    glBegin(GL_POINTS);


    typename sofa::helper::vector<potentialContact>::iterator it = m_VecPotentialContact.begin();
    for (unsigned int p=0; p<m_posSample.size(); p++)
    {


        if(m_VecPotentialContact.size()>0 && (*it).posSampleId==p) // in potiential contact
        {
            glColor4f(1.0f,0.0f,0.0f,1.0f);
            helper::gl::glVertexT(m_posSample[p]);

            if (it!=m_VecPotentialContact.end())
                it++;
        }
        else
        {
            glColor4f(0.0f,1.0f,1.0f,1.0f);
            helper::gl::glVertexT(m_posSample[p]);
        }

    }
    glEnd();
    glPointSize(1);

    glBegin(GL_LINES);
    for (unsigned int i=0; i<m_VecPotentialContact.size(); i++)
    {
        glColor4f(1.0f,0.0f,0.0f,1.0f);
        potentialContact &pt= m_VecPotentialContact[i];
        Vec3 Pos = m_posSample[pt.posSampleId];
        helper::gl::glVertexT(Pos);
        helper::gl::glVertexT(Pos + pt.n);

        if(friction)
        {
            glColor4f(0.0f,1.0f,0.0f,1.0f);
            helper::gl::glVertexT(Pos);
            helper::gl::glVertexT(Pos + pt.t);

            glColor4f(0.0f,0.0f,1.0f,1.0f);
            helper::gl::glVertexT(Pos);
            helper::gl::glVertexT(Pos + pt.s);
        }

    }
    glEnd();



    glEnable(GL_LIGHTING);
    // glDisable(GL_BLEND);

    if (visualization.getValue() && m_posSample.size()>0)
    {

        std::cout<<" Visualization "<<std::endl;
        std::cout<<" m_domainSample"<<m_domainSample<<std::endl;

        double dx = (PosMax.getValue()[0] -PosMin.getValue()[0])/50;
        double dy = (PosMax.getValue()[1] -PosMin.getValue()[1])/50;
        double dz = (PosMax.getValue()[2] -PosMin.getValue()[2])/50;
        double z = PosMin.getValue()[2];
        int idx =0;

        double mmin = 150;
        double mmax = -150;

        defaulttype::Vec3d PosLastPoint=m_posSample[m_posSample.size()-1];

        //std::cout<<" mc_data =";
        while(z<=PosMax.getValue()[2])
        {
            double y = PosMin.getValue()[1];
            while(y<=PosMax.getValue()[1])
            {
                double x = PosMin.getValue()[0];


                while(x<=PosMax.getValue()[0])
                {
                    defaulttype::Vec3d Pos(x,y,z);
                    Pos += PosLastPoint;
                    //Pos += PP;
                    double test = _contact_surface->getValue(Pos, m_domainSample[m_posSample.size()-1]);
                    if(test > mmax)
                        mmax = test;
                    if(test < mmin)
                        mmin = test;
                    mc_data[idx] = test;
                    //std::cout << test << " ";
                    idx++;

                    x = x+dx;
                }
                y = y+dy;

            }
            z = z+dz;
        }

        glPolygonMode(GL_FRONT, GL_LINE);


#ifdef SOFAEVE
        mc->buildMesh(mc_data, 50, 50, 50, _isoValue);
#endif

        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);


        defaulttype::Vec3d Pos(PosMin.getValue()[0],PosMin.getValue()[1],PosMin.getValue()[2]);
        Pos += PosLastPoint;
#ifdef SOFAEVE
        mc->draw(Pos[0], Pos[1], Pos[2], 0.2, 0.2, 0.2);
#endif



    }



#ifdef DEBUG
    std::cout<<" leaving draw"<<std::endl;
#endif


}







////////////////// FONCTORS RESOLUTION ///////////////////////////////////////////


template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraintResolution<DataTypes>::resolution(int line, double** w, double* d, double* force)
{

    // dummy resolution
    d[line]=w[line][line]*force[line];
    force[line]=0.0;

    /*
    sofa::helper::system::thread::CTime *timer;
    timer = new sofa::helper::system::thread::CTime();
    double time = (double) timer->getTime();


    // TODO : faire l'init en associant le probe

    helper::vector<ImplSurfContact>& lp = *(_implSurfConstraint->listProbe.beginEdit());



    ImplSurfContact* probe = &lp[_id_probe];
    _implSurfConstraint->listProbe.endEdit();

    // the point that is targeted by the resolution is set in _id_probe
    if (_first)
    {
        // place the value of w in the local value of the probe
        probe->setW(w[line][line], w[line][line+1], w[line][line+2], w[line+1][line+1], w[line+1][line+2], w[line+2][line+2]);
        // provide the initial guess
#ifdef DEBUG
                std::cout<<"W = ["<<w[line][line]<<" "<<w[line][line+1]<<" "<<w[line][line+2]<<"\n "
                        <<w[line+1][line]<<" "<<w[line+1][line+1]<<" "<<w[line+1][line+2]<<"\n "
                        <<w[line+2][line]<<" "<<w[line+2][line+1]<<" "<<w[line+2][line+2]<<"]"<<std::endl;
#endif
        force[line  ] = probe->saveForce()[0];
        force[line+1] = probe->saveForce()[1];
        force[line+2] = probe->saveForce()[2];
        _first=false;
        probe->initStep();
        probe->storeBeginInfo();

    }
    //else
    //{

        // computation of the current position of the probe //
        defaulttype::Vec3d pos;
        pos = probe->freePos();

        pos[0] += d[line]; pos[1] += d[line+1]; pos[2] += d[line+2];
                //defaulttype::Vec3d saveForce = probe->saveForce() ;
                //pos += probe->computeDisp(saveForce);   //> NOT NECESSARY ANY MORE !!
        probe->setPos(pos);


#ifdef DEBUG
                std::cout<<"ImplicitSurfaceAdaptiveConstraintResolution : pos = "<<pos<<std::endl;
#endif

        // resolution //
        probe->solveConstraint2(_mu);

        probe->storeEndInfo();



        // resulting force;
        force[line  ] = probe->force()[0];
        force[line+1] = probe->force()[1];
        force[line+2] = probe->force()[2];

        // keep the value of the force in memory
        probe->updateForce();
    //}



    TimeResolution += (double) timer->getTime() - time;

    */

}








} // namespace constraint

} // namespace component

} // namespace sofa

#endif
