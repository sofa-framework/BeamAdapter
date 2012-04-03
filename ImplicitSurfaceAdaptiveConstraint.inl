#ifndef SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_INL
#define SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_INL


#include "ImplicitSurfaceAdaptiveConstraint.h"

#include <sofa/defaulttype/Vec.h>
#include <sofa/helper/gl/template.h>


#include <sofa/helper/io/MassSpringLoader.h>
#include <sofa/helper/gl/template.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/DataTypeInfo.h>
#include <sofa/component/topology/RegularGridTopology.h>
#include <sofa/component/mass/AddMToMatrixFunctor.h>
#include <sofa/component/topology/TopologyData.inl>

#include <sofa/helper/system/thread/CTime.h>


double TimeResolution = 0.0;
double TimeCount = 0.0;
double TimeProjection = 0.0;
double TimeProjection2= 0.0;
//#define DEBUG
//#define DEBUG_TOPO_DATA
//#define DEBUG_PROJECTION
//#define DEBUG_LAST_CONSTRAINT_ONLY


namespace sofa
{

namespace component
{

namespace constraint
{


defaulttype::Vec3d Probe::computeDisp(defaulttype::Vec3d &f)
{
    defaulttype::Vec3d result(0.0,0.0,0.0);

    for (int i=0; i<3; i++)
    {
        for (int j=0; j<3; j++)
            result[i] += _w[i][j]*f[j];
    }
    return result;

}



ImplSurfContact::ImplSurfContact( sofa::component::container::ImplicitSurface *ImplicitSurface)
{
    _inContact = false;
    _saveForce.clear();
    _implSurf = ImplicitSurface;
    _domain = -1;
    _friction = true;
}


void ImplSurfContact::noContact()
{

    _force.clear();
    _saveForce.clear();
    _fn = 0.0;
    _ft = 0.0;
    _fs = 0.0;
    _inContact = false;
    _savePos = _pos;
#ifdef DEBUG_PROJECTION
    std::cout<<"no Contact is called at pos:"<<_pos<<std::endl;
#endif
}


void ImplSurfContact::getOrthogonalVectors(const defaulttype::Vec3d& dir, defaulttype::Vec3d& vec1, defaulttype::Vec3d& vec2)
{
    defaulttype::Vec3d temp;	// Any vector such as temp != dir
    temp[0] = dir[1];
    temp[1] = dir[2];
    temp[2] = dir[0];

    if(temp == dir) // x = y = z
        temp = defaulttype::Vec3d(1,0,0);

    vec1 = cross(dir, temp);
    vec1.normalize();

    vec2 = cross(dir, vec1);
    vec2.normalize();
}

void ImplSurfContact::rotateFrame(const defaulttype::Vec3d& n_old, const defaulttype::Vec3d& n_new, defaulttype::Vec3d& t, defaulttype::Vec3d& s)
{
    defaulttype::Vec3d axis = cross(n_old, n_new);
    double dot_prod = dot(n_old, n_new);


    if (axis.norm() > 0.0001 && dot_prod < 0.9999)
    {

        double angle = acos(dot_prod);
        //std::cout<<"quat axis :"<<axis<<"- angle :"<< angle <<" dot(n_old, n_new) :"<<dot(n_old, n_new)<<std::endl;
        axis.normalize();

        helper::Quater<double> q(axis, angle);
        //n_new = q.rotate(n_old);
        t = q.rotate(t);
        s = q.rotate(s);

        //std::cout<<"quat t :"<<t<<" - s :"<<s<<std::endl;
    }

    //std::cout<<"cross n_new :"<<n_new<<"- n_old :"<< n_old <<" - s :"<<s<<std::endl;
    t = - cross(n_new, s);
    t.normalize();

    s = cross(n_new, t);
    s.normalize();

    // debug
    if (t.norm()>1.00001 || t.norm()<0.9999 || s.norm()>1.00001 || s.norm()<0.9999)
    //	std::cout<<"ERROR in rotate Frame : t"<<t<< " s "<<s<<std::endl;

    return;



}

void ImplSurfContact::changeContactFrame(defaulttype::Vec3d& posSurf)
{

    defaulttype::Vec3d n_old = _n;

    defaulttype::Vec3d grad = _implSurf->getGradient(posSurf, _domain);
    grad.normalize();
    _n = grad;
    _n.normalize();

    rotateFrame(n_old, _n, _t, _s);
}

bool ImplSurfContact::detectNewContact()
{
    //std::cout << "Detection _pos =" << _pos << "Domain = " << _domain << std::endl;
    double valueBuf = _implSurf->getValue(_pos, _domain);
    //std::cout << "Value : " << valueBuf << std::endl;
    // the probe was not in contact
    if(valueBuf > 0)			// TODO : definir un epsilon / choix exterieur-interieur
    {	// case 1: was not in contact and is still not in contact
#ifdef DEBUG
        std::cout<<"case1 :  value ="<<valueBuf<<std::endl;
#endif
        this->noContact();

        return false;
    }
    else
    {	// case 2: was not in contact but is now colliding the implicit surface
#ifdef DEBUG
        std::cout<<"case2 : ";
#endif
        // we try to find the colinding point on the surface
        _inContact = true;
        // debug //
        if(_implSurf->getValue(_savePos, _domain) < 0) // TODO : choix exterieur-interieur
        {
            std::cout<<"Problem : the probe was not supposed to be in contact during the previous step"<<std::endl;

            return false;
        } // end debug //

        if(!_implSurf->computeSegIntersection(_pos, _savePos, _posSurf, _domain))
            std::cout<<"Problem: no intersection found"<<std::endl;

        if (_posSurf.norm() < 0.001)
        {
            std::cout<<"Problem with _posSurf   : _pos "<<_pos << " _savePos "<< _savePos <<std::endl;
        }



        _n = _implSurf->getGradient( _posSurf , _domain);  // TODO : definir choix exterieur-interieur _n = - _implSurf->getGradient( _posSurf ); si on fait des contact à l'interieur de la surface implicite
        _n.normalize();


#ifdef DEBUG
        std::cout<<"Pos intersection = "<<_posSurf<< "normale"<< _n<< std::endl;
#endif
        // TODO _friction : computation of _t and _s
        if (_friction){
            getOrthogonalVectors(_n, _t, _s);
            _posBegin=_posSurf;
        }


        return true;

    }


}


void ImplSurfContact::storeBeginInfo()
{
    if(_inContact)
    {
        changeContactFrame(_posSurf);
        _posBegin = _posSurf;
#ifdef DEBUG_TOPO_DATA
        std::cout<<"begin info: in contact"<<std::endl;
#endif
/*		_nBegin = _n;
        _tBegin = _t;
        _sBegin = _s;*/
        _domain = _implSurf->getDomain(_posSurf, _domain);
    }
    else
    {
#ifdef DEBUG_TOPO_DATA
        std::cout<<"begin info: no contact"<<std::endl;
#endif
        _domain = _implSurf->getDomain(_pos, _domain);
    }

#ifdef DEBUG
    std::cout<<"Begin info : pos:"<<_pos<<" - _posSurf"<<_posSurf<<"_n"<<_n <<"_t"<<_t<< "_s"<<_s<<std::endl;
#endif

    /*if(_domain < 0)
    {
        std::cout<<"WARNING problem : Domain ="<< _domain<<"  at pos "<<_pos<<std::endl;
    }*/



}





// Gauss-Seidel iterative call
bool ImplSurfContact::solveConstraint2(double &mu)
{
#ifdef DEBUG
    std::cout<<"solveConstraint2 domain"<<this->_domain<<std::endl;
#endif

    int numItMax = 1;

    if (mu>0.0)
        _friction = true;
    else
        _friction = false;

// DEBUG:
    defaulttype::Vec3d posWithNoForce;


    // test: Was the probe already in contact or not ?
    if(!_inContact )
    {
        if(!detectNewContact())
        {
            //std::cout<<" no contact - domain ="<<this->_domain<<std::endl;
            return false;


        }

        posWithNoForce = _pos;

    }
    else
    {
        // we were in contact...
        posWithNoForce = _pos - computeDisp(_force); //
    }




    _wnn = compute_wdd(_n);
    _wtt = compute_wdd(_t);
    _wss = compute_wdd(_s);


    //std::cout << "solveConstraint2 - posWithNoForce:" << posWithNoForce<<" - _wnn ="<<_wnn<<" - _wtt ="<<_wtt<<" - _wss ="<<_wss<<std::endl;
    int it;
    defaulttype::Vec3d posBuf;
    defaulttype::Vec3d pos_posBuf;
    for (it=0; it<numItMax ; it++)
    {
        //
        posBuf = _pos;
        //
        defaulttype::Vec3d DPos = _pos - _posSurf;
        double dn = dot(_n,DPos);
        _fn -= dn / _wnn;

        if (_fn < 0)
        {
            _force.clear();
            _fn = 0.0;
            _ft = 0.0;
            _fs = 0.0;
            _pos = posWithNoForce;
            return false;
        }
        // contact force:
        _force = _n * _fn +  _t * _ft + _s * _fs;

        if (_friction)
        {
            //
            _pos = posWithNoForce + computeDisp(_force);
            //
            DPos =_pos - _posBegin;
            double dt = dot(_t,DPos);
            double ds = dot(_s,DPos);


            _ft -= 2*dt/(_wtt+_wss);
            _fs -= 2*ds/(_wtt+_wss);
            double normForceTg = sqrt(_ft*_ft+_fs*_fs);
            if( normForceTg > mu * _fn )
            {
                // glissement
                _ft *= mu * abs(_fn)/normForceTg;
                _fs *= mu * abs(_fn)/normForceTg;
            }

            // contact force:
            _force = _n * _fn +  _t * _ft + _s * _fs;
        }

        _pos = posWithNoForce + computeDisp(_force);

        pos_posBuf = _pos-posBuf;

        if(pos_posBuf.norm()<0.001)
            break;
    }

    /*
    if (it==numItMax && pos_posBuf.norm()>0.01)
    {
        std::cout<<"\n no convergence : error  = "<<pos_posBuf.norm() ;
    }
    */

    // on finit toujours par une correction selon n pour être sûr de ne pas avoir de pb de projection
    defaulttype::Vec3d DPos = _pos - _posSurf;
    double dn = dot(_n,DPos);
    _fn -= dn / _wnn;
    _force = _n * _fn +  _t * _ft + _s * _fs;
    _pos = posWithNoForce + computeDisp(_force);



    return true;



}



void
ImplSurfContact::storeEndInfo()
{

    sofa::helper::system::thread::CTime *timer;
    timer = new sofa::helper::system::thread::CTime();
    double time = (double) timer->getTime();

    // dans le cas où on ne le fait pas avant...
    projectSurfacePos();
    TimeCount += (double) timer->getTime() - time;

    if(_inContact)
    {
#ifdef DEBUG_TOPO_DATA
        std::cout<<"end info: in contact"<<std::endl;
#endif
        //changeContactFrame(_posSurf);
    }
#ifdef DEBUG_TOPO_DATA
    else
        std::cout<<"end info: no contact"<<std::endl;
#endif



    /* on remet à jour la normale */
    //changeContactFrame(_posSurf);

}


void
ImplSurfContact::projectSurfacePos()
{

    sofa::helper::system::thread::CTime *timer;
    timer = new sofa::helper::system::thread::CTime();
    sofa::helper::system::thread::CTime *timer2;
    timer2 = new sofa::helper::system::thread::CTime();

    defaulttype::Vec3d newPos;
    double value;

    #ifdef DEBUG
    std::cout<<"projectSurfacePos : _pos:"<< _pos<<"  - value:" << value<<" - _n:"<< _n<<" - _posSurf:"<<_posSurf<<std::endl;
    #endif

    // if the point is not in contact (case 1), it is not necessary to project it !
    if(!_inContact)
    {
#ifdef DEBUG_PROJECTION
        if(value< 0.0)
            std::cout<<"WARNING : point in collision is not activated"<<std::endl;
#endif
        return;
    }


    value= _implSurf->getValue(_pos, _domain);
    // _inContact is activated (case 2 - 3 - 4)
    if(_fn>0.0)
    {
        _inContact = true;

/* Ce probe est maintenu en contact, on sépare deux cas:
    1. point à l'intérieur (value <0 )
        on est dans un cas concave (ou le gauss-seidel a laissé une erreur)
        le point sur la surface est trouvé en utilisant la normale à la surface


    2. point à l'extérieur (value >0)
        on est dans un cas convexe (ou alors en adherence et le gauss seidel a bien convergé)
        on regarde  le point obtenu:
            - si le point est proche de la surface, on garde la contrainte active et on projette le point à l'aide de la normale à la surface
            - si le point est loin de la surface, on supprime la contrainte
        on regarde  le point projeté:

*/

        if (value<0)
        {
#ifdef DEBUG_PROJECTION
            std::cout<<"concave surface case"<<std::endl;
#endif
            // cas concave //
            // on place _posSurf hors de l'objet (à une distance de 0.1, par exemple)
            double dist_out = 0.001;
            if (! _implSurf->projectPointOutOfSurface(_posSurf,_domain, _n, dist_out))
            {

                    std::cout<<"\n projectPointOutOfSurface does not work !!";
                    return;
            }


            if (_implSurf->getValue(_posSurf, _domain)<0)
                std::cout<<"\n WARNING : _posSurf is not placed out of the surface ! in projectSurfacePos";

            // on cherche le point d'intersection entre ce nouveau _posSurf et la position actuelle, c'est notre nouveau point de reference

            _implSurf->computeSegIntersection(_pos,_posSurf,newPos, _domain);
            _posSurf = newPos;

            _pos = newPos;	// utile ?


        }
        else
        {
#ifdef DEBUG_PROJECTION
            std::cout<<"convex surface case or stick point"<<std::endl;
#endif
            // cas convexe
            if (value >1.0)
            {
                std::cout<<"WARNING: Very convex case "<<std::endl;
                // cas très convexe, évaluation du gradient difficile
                noContact();
                //_saveForce.reset();// ici, il faut mettre saveForce à zero (utilisé dans le calcul addSaveForce(force, torque) au début du calcul
                return;
            }
            else
            {


                /*

                double time2 = (double) timer2->getTime();
                // evaluation du gradient à la position _pos
                defaulttype::Vec3d dir = _implSurf->getGradient(_pos,_domain); // on remonte à la surface en utilisant le gradient
                dir.normalize();
                newPos = _pos;



                if(!_implSurf->projectPointonSurface2(newPos, _domain, dir)) // on projette la position actuelle sur la surface
                {
                    std::cout<<"\n WARNING no projection found for smooth convex case";
                    noContact();
                    //_saveForce.reset();// ici, il faut mettre saveForce à zero (utilisé dans le calcul addSaveForce(force, torque) qui vient après storeBeginInfo)
                    return;

                }
                TimeProjection2 += (double) timer2->getTime() - time2;

                */



                double time = (double) timer->getTime();
                newPos = _pos;
                // cas "fixe"
                defaulttype::Vec3d Disp = _pos - _posSurf;

                if(Disp.norm() > 0.001)    //TODO : mettre en paramètre
                {
                    _implSurf->projectPointonSurface(newPos, _domain);
                    //_implSurf->projectPointonSurface2(newPos, _domain, dir);
                }

                _posSurf = newPos;
                _pos = newPos;	// utile ?
                TimeProjection += (double) timer->getTime() - time;


            }

        }






    }
    else
    {

/* Le probe n'est pas en contact, on sépare deux cas:
    1. le point est malgré tout en collision
        - le point n'était pas dans la liste des contacts "potentiels" (pas normal)
        - cas concave : la direction de la contrainte a évolué

    2. le point n'est pas en collision (rien à faire)

*/


        if(value < 0)
        {
            std::cout<<" WARNING the actual _pos is colliding and no collision was detected"<<std::endl;

            if (!_inContact)
            {
                defaulttype::Vec3d depl = _savePos - _pos;
                std::cout<<" WARNING the actual _pos is colliding and no collision was detected: Depl = "<< depl.norm();
            }

            double dist_out = 0.1;
            if(!_implSurf->projectPointOutOfSurface(_posSurf,_domain, _n, dist_out))
            {
                std::cout<<" pb for non contact case \n";
            }
            else
            {
                defaulttype::Vec3d newPos;
                _implSurf->computeSegIntersection(_pos,_posSurf,newPos, _domain);
                _pos = newPos;
            }

        }
        //else
        //{

        //}
#ifdef DEBUG
        std::cout<<"simple no contact case"<<std::endl;
#endif
        noContact();
        //_inContact = false;

    }



}




void
ImplSurfContact::draw()
{


    if(_inContact){
    glBegin(GL_LINES);
        helper::gl::glVertexT(_pos);
        helper::gl::glVertexT(_posSurf);
    glEnd();


    // draw force direction //

    glColor4f(1.0f,1.0f,0.0f,1.0f);



    defaulttype::Vec3d n_force = _force;
    n_force.normalize();
    glBegin(GL_LINES);
        helper::gl::glVertexT(_pos);
        helper::gl::glVertexT(_pos + n_force);
    glEnd();

    }


    //
}






//////////////////////////////
// ImplicitSurfaceAdaptiveConstraint
//////////////////////////////



template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::init()
{
std::cout<<"ImplicitSurfaceAdaptiveConstraint<DataTypes>::init()"<<std::endl;
    assert(this->object1);
    assert(this->object2);
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
    helper::vector<ImplSurfContact>& lp = *(listProbe.beginEdit());
    lp.clear();
    listProbe.endEdit();

    // get the implicit surface
    this->object1->getContext()->get(_contact_surface);

    if (_contact_surface == NULL){
        serr<<"oooooooooo\n oooERROR: no surface found for contact"<<sendl;
        return;
    }
    else
    {
        std::cout<<"CONTACT with surface name:"<<_contact_surface->getName()<<std::endl;
    }


    // we create three constraint at each point of object2
    // even in the case of _frictionless contact, in order to have the value of the compliance in the three global direction
    int m2 = this->object2->getSize();
    const VecCoord& x		= (*this->object2->getX());
#ifdef DEBUG_LAST_CONSTRAINT_ONLY
    _nbConstraints = 1;
#else
    _nbConstraints = m2;
#endif
    std::cout<<"NUM PROBES = "<<m2<<std::endl;

    //VecCoord x = (*this->object2->getX());
    lp = *(listProbe.beginEdit());
    lp.resize(m2);
    active_list.clear();

    for (int p=0; p<m2; p++){

        //ImplSurfContact* probe = new ImplSurfContact(contact_surface);
        //listProbe.push_back(probe);
        lp[p].setImplicitSurface(_contact_surface);

        CPos P = DataTypes::getCPos(x[p]);
        lp[p]._domain = _contact_surface->getDomain(P, -1); //
        lp[p].setPos(P);
        //std::cerr<<"domain found at position : "<< P<<" is "<<lp[p]._domain<<std::endl;

        init_domain.setValue(lp[p]._domain);

        active_list.push_back(false);

        //std::cout << "LP domain " << lp[p]._domain <<" pos["<<p<<"] ="<< lp[p]._pos<< std::endl;
    }

    listProbe.endEdit();
    //_radiusSphere = 1.0;
}

/*
template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getOrthogonalVectors(const Deriv& dir, Deriv& vec1, Deriv& vec2)
{

}
*/



/*

template <class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::handleTopologyChange()
{

    #ifdef DEBUG_TOPO_DATA
        std::cout<<"handleTopologyChange before : ";

        int m2 = this->object2->getSize();
        for (int i=0; i<m2-1; i++)
        {
            std::cout<<"listProbe["<<i<<"].pos :"<<listProbe[i].pos()<<std::endl;
            std::cout<<"                          ";
        }
        std::cout<<""<<std::endl;
    #endif

    core::topology::BaseMeshTopology* topology = this->object2->getContext()->getMeshTopology();
    std::list<const core::topology::TopologyChange *>::const_iterator itBegin=topology->beginChange();
    std::list<const core::topology::TopologyChange *>::const_iterator itEnd=topology->endChange();
    std::list<const core::topology::TopologyChange *>::const_iterator it;
    listProbe.handleTopologyEvents(itBegin,itEnd);

    #ifdef DEBUG_TOPO_DATA
         m2 = this->object2->getSize();
        for (int i=0; i<m2; i++)
        {
            std::cout<<"                     after  : listProbe["<<i<<"].pos :"<<listProbe[i].pos()<<std::endl;
        }
        std::cout<<""<<std::endl;


    #endif

}

*/

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST =ConstraintParams::defaultInstance()*/, core::MultiMatrixDerivId cId, unsigned int &constraintId)
{

    std::cout<<"ImplicitSurfaceAdaptiveConstraint<DataTypes>::buildConstraintMatrix()"<<std::endl;


    if (!cParams)
    {
        serr<<" WARNING cParams not defined can not build constraint matrix"<<sendl;
        return;
    }

    DataMatrixDeriv &cData = *cId[object2].write();
    MatrixDeriv& c2 = *cData.beginEdit();


    cid = constraintId;
    const VecCoord& x  = (*this->object2->getX());
    const VecCoord& xfree	= (*this->object2->getXfree());
    //MatrixDeriv& c1 = *this->object1->getC();
    //MatrixDeriv& c2 = *this->object2->getC();



    unsigned int _nbNodes = x.size();

    // option 1:

#ifdef DEBUG_LAST_CONSTRAINT_ONLY
    _nbConstraints = 1;
#else
    _nbConstraints = _nbNodes;
#endif


    // each point creates a constraint (but the constraint can be not active)
    //_nbConstraints = x.size();
    Deriv vx,vy,vz;
    DataTypes::setDPos(vx, DPos(1.0,0.0,0.0));
    DataTypes::setDPos(vy, DPos(0.0,1.0,0.0));
    DataTypes::setDPos(vz, DPos(0.0,0.0,1.0));

    helper::vector<ImplSurfContact>& lp = *(listProbe.beginEdit());

    std::cout<<" _nbConstraints ="<<_nbConstraints<<std::endl;

    active_list.clear();
#ifdef DEBUG_LAST_CONSTRAINT_ONLY
    for (int p=	_nbNodes-1 ; p<_nbNodes; p++)
    {
#else
    for (int p=0; p<_nbConstraints; p++)
     {
#endif
        // initialize with domain = init_domain;
        if (lp[p]._domain<0)
            lp[p]._domain = init_domain.getValue();
        ////////////// update the probes ////////////
        CPos Pfree = DataTypes::getCPos(xfree[p]);
        lp[p].setPosFree(Pfree);
        CPos P = DataTypes::getCPos(x[p]);
        lp[p].setPos(P);



        if(cancelPointThreshold.getValue()>0.0)
        {
            if (_contact_surface->getValue(P, lp[p]._domain) > cancelPointThreshold.getValue())
            {
                //std::cout<<" point "<<p<<" is not active : value = "<<_contact_surface->getValue(P, lp[p]._domain) <<std::endl;
                active_list.push_back(false);
                continue;
            }
            else
                active_list.push_back(true);
        }
        else
        {

            active_list.push_back(true);

        }

        lp[p].setImplicitSurface(_contact_surface);


        ////////////// create 3 constraints along each axis x, y and z ///////
        //SparseVecDeriv svd1;
        //SparseVecDeriv svd2;
        //x
//#ifdef DEBUG_LAST_CONSTRAINT_ONLY
//        this->object2->setConstraintId(cid );
//#else
//        this->object2->setConstraintId(cid + 3*p);
//#endif
        //svd2.add(p, vx);
        //c2.push_back(svd2);

        MatrixDerivRowIterator xConstIt = c2.writeLine(constraintId);
        xConstIt.setCol(p, vx);
        //y
//#ifdef DEBUG_LAST_CONSTRAINT_ONLY
//        this->object2->setConstraintId(cid + 1);
//#else
//        this->object2->setConstraintId(cid + 3*p+1);
//#endif
        //svd2.set(p, vy);
        //c2.push_back(svd2);
        MatrixDerivRowIterator yConstIt = c2.writeLine(constraintId+1);
        yConstIt.setCol(p, vy);
        //z
//#ifdef DEBUG_LAST_CONSTRAINT_ONLY
//        this->object2->setConstraintId(cid + 2);
//#else
//        this->object2->setConstraintId(cid + 3*p+2);
//#endif
        //svd2.set(p, vz);
        //c2.push_back(svd2);
        MatrixDerivRowIterator zConstIt = c2.writeLine(constraintId+2);
        zConstIt.setCol(p, vz);


        constraintId += 3;

   }

    listProbe.endEdit();
    cData.endEdit();

}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getConstraintValue(double* v)
{
    std::cout<<"\n \n ************* is it really used ? \n \n **************"<<std::endl;
    for(int i=0; i<_nbConstraints; i++)
    {
        if(active_list[i])
            v[cid+i] = 0.0;
    }
}
/*
template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getConstraintId(long* , unsigned int &)
{

    // not very useful here


}
*/

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset)
{


    //
#ifdef DEBUG_LAST_CONSTRAINT_ONLY
    const VecCoord& x		= (*this->object2->getX());
    unsigned int _nbNodes = x.size();

    for(int i=_nbNodes-1; i<_nbNodes; i++){
        int it = 0;

#else

      //  std::cout << "Nb Constraints :" << _nbConstraints << " nbNodes = " << _nbNodes << std::endl;
    for(int i=0; i<_nbConstraints; i++)	{
        int it=i;
#endif
        //std::cout << "General case" << std::endl;

        //printf("Tab = 0x%x\n", resTab.size());


        if(this->active_list[it])
        {
            resTab[offset] = new ImplicitSurfaceAdaptiveConstraintResolution<DataTypes>(mu.getValue(), this, i);

            offset+=3;
        }
    }
    //std::cerr<<"ImplicitSurfaceAdaptiveConstraint<DataTypes>::getConstraintResolution() _nbConstraints ="<<_nbConstraints<<std::endl;


}

template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraint<DataTypes>::draw(const core::visual::VisualParams* vparams)
{

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
    const VecCoord& x		= (*this->object2->getX());
    unsigned int _nbNodes = x.size();
    //std::cout << "Drawing centered on " << x << std::endl;
#endif


    if (visualization.getValue() && _nbConstraints>0)
    {

        glDisable(GL_LIGHTING);
        glPointSize(2);
        glBegin(GL_POINTS);
        double dx = (PosMax.getValue()[0] -PosMin.getValue()[0])/10;
        double dy = (PosMax.getValue()[1] -PosMin.getValue()[1])/10;
        double dz = (PosMax.getValue()[2] -PosMin.getValue()[2])/10;
        double x = PosMin.getValue()[0];



        helper::vector<ImplSurfContact>& lp = *(listProbe.beginEdit());



        while(x<=PosMax.getValue()[0])
        {
            double y = PosMin.getValue()[1];
            while(y<=PosMax.getValue()[1])
            {
                double z = PosMin.getValue()[2];


                while(z<=PosMax.getValue()[2])
                {
                    defaulttype::Vec3d Pos(x,y,z);
#ifdef DEBUG_LAST_CONSTRAINT_ONLY
                    Pos += lp[_nbNodes-1].pos();
                    if(_contact_surface->getValue(Pos,lp[_nbNodes-1]._domain ) > 0){
#else
                    Pos += lp[_nbConstraints-1].pos();
                    if(_contact_surface->getValue(Pos,lp[_nbConstraints-1]._domain ) > 0){
#endif
                        glColor4f(1.0f,0.0f,0.0f,1.0f);
                        helper::gl::glVertexT(Pos);
                        }
                    else{
                        glColor4f(0.0f,1.0f,1.0f,1.0f);
                        helper::gl::glVertexT(Pos);
                        }

                    z = z+dz;
                }
                y = y+dy;

            }
            x = x+dx;
        }
        glEnd();
        glPointSize(1);
        listProbe.endEdit();

    }

    //if (!context->getShowInteractionForceFields()) return;

    glDisable(GL_LIGHTING);
    glPointSize(10);
    glBegin(GL_POINTS);

    /*
    int m = this->object1->getSize();
    glColor4f(0,0,1,1);
    for(int i=0; i<m; i++)
        helper::gl::glVertexT((*this->object1->getX())[i]);
    */

    int m = this->object2->getSize();

    helper::vector<ImplSurfContact>& lp = *(listProbe.beginEdit());


    for(int i=0; i<m && lp.size(); i++)
    {
        if (lp[i].inContact())
            glColor4f(1.0f,0.0f,0.0f,1.0f);
        else
            glColor4f(0.0f,0.0f,1.0f,1.0f);
        helper::gl::glVertexT((*this->object2->getX())[i]);
    }

    glEnd();
    glPointSize(1);

    // line between _pos and _posSurf

    for(int i=0; i<_nbConstraints; i++)
    {

        if(this->active_list[i])
            lp[i].draw();
    }
    listProbe.endEdit();

    /*
    glBegin(GL_LINES);

    m = this->object1->getSize();
    for(int i=1; i<m; i++)
    {
        float c = 1.0f/((float)(m-1)*(i-1));
        glColor4f(c,0.0f,1.0f-c,1.0f);
        helper::gl::glVertexT((*this->object1->getX())[i-1]);
        c = 1.0f/((float)(m-1)*i);
        glColor4f(c,0.0f,1.0f-c,1.0f);
        helper::gl::glVertexT((*this->object1->getX())[i]);
    }
    glEnd();
    */


}







////////////////// FONCTORS RESOLUTION ///////////////////////////////////////////


template<class DataTypes>
void ImplicitSurfaceAdaptiveConstraintResolution<DataTypes>::resolution(int line, double** w, double* d, double* force)
{

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

        /* si on change de normale à chaque iteration */
        probe->storeEndInfo();



        // resulting force;
        force[line  ] = probe->force()[0];
        force[line+1] = probe->force()[1];
        force[line+2] = probe->force()[2];

        // keep the value of the force in memory
        probe->updateForce();
    //}



    TimeResolution += (double) timer->getTime() - time;

}








} // namespace constraint

} // namespace component

} // namespace sofa

#endif
