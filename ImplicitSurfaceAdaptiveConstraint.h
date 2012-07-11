#ifndef SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_H

#include <sofa/core/behavior/BaseInteractionConstraint.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/visual/VisualModel.h>
#include <sofa/helper/gl/template.h>
#include <iostream>
#include <sofa/component/container/ImplicitSurfaceContainer.h>

#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/helper/vector.h>
#include <sofa/component/topology/TopologyData.h>

namespace sofa
{

namespace component
{

namespace constraint
{

/*!
 * \class Probe
 * \brief Parent class for all Probe type
 */
class Probe
{
public:

        Probe(){}
        virtual ~Probe(){}

       // get-set-update functions//
        inline const defaulttype::Vec3d& freePos() const { return _freePos; }
        inline const defaulttype::Vec3d& pos() const { return _pos; }
        inline const defaulttype::Vec3d& savepos() const { return _savePos; }
        inline const defaulttype::Vec3d& saveForce() const { return _saveForce; }
        inline const defaulttype::Vec3d& force() const { return _force; }

        inline void setPos(defaulttype::Vec3d &p)		{ _pos = p; }
        inline void setPosFree(defaulttype::Vec3d &p)	{ _freePos = p; }
        inline void setSavePos(defaulttype::Vec3d &p) { _savePos = p; }

        inline void setW(double a, double b, double c, double d, double e, double f){ _w[0][0]=a; _w[0][1]=b; _w[0][2]=c;
                 _w[1][0]=b; _w[1][1]=d; _w[1][2]=e;
                  _w[2][0]=c; _w[2][1]=e; _w[2][2]=f;}

        double compute_wdd(defaulttype::Vec3d& dir){
            double result = dir[0]*(_w[0][0]*dir[0] + 2*_w[0][1]*dir[1] + 2*_w[0][2]*dir[2]) +
                            dir[1]*(				    _w[1][1]*dir[1] + 2*_w[1][2]*dir[2]) +
                            dir[2]*(									+ _w[2][2]*dir[2]);
            return result;
        }


        //!< update the free position of the probe using the 6D position of the tool (beginning of Gauss-Seidel it)
        //inline void updateFreePos(defaulttype::Vec3d& center_free, Quaternion& orientation_free) {_freePos = center_free + orientation_free.rotate(_refPos) ;}

        //inline void updateLever(Quaternion& orientation) {_lever = orientation.rotate(_refPos);}
        //inline void updatePos(defaulttype::Vec3d& center_depl, defaulttype::Vec3d& orientation_depl) {_pos = _freePos + center_depl - cross(_lever,orientation_depl);}
        /*inline void updateForce(defaulttype::Vec3d &center_force, defaulttype::Vec3d &orientation_force) {df = _force - _saveForce;
                                                                                center_force += df;
                                                                                orientation_force += cross(_lever,df);
                                                                                _saveForce = _force;}*/
        inline void updateForce() {_saveForce = _force;}
        //!< compute the local compliance matrix on the Probe
        //void updateCompliance(Quaternion& orientation, Vec3 &compTrans, Vec3 &compRot);

        //!< compute the local displacement due to a given force
        defaulttype::Vec3d computeDisp(defaulttype::Vec3d &f);

        // virtual function: their implementation depends on the type of constraint
       virtual void initStep()=0;

       virtual void draw()=0;





public:

        //!< Note: these references do not change after contruction
        defaulttype::Vec3d _pos;			//!< position in global coordinates (updated during Gauss-Seidel)
        defaulttype::Vec3d _freePos;	//!< free position in global coordinate (computed at the beginning of each time-step)
        defaulttype::Vec3d _savePos;	    //!< Save position in global coordinates (updated when Gauss-Seidel algo had converged) : it is used as the position from the previous time-step
        defaulttype::Vec3d _force;		//!< Force reaction on the probe
        defaulttype::Vec3d _saveForce;	//!< Save force reaction (updated at the end of "updateForce")

        double _w[3][3];		//!< Local compliance of the probe


};





/*!
 * \class ImplSurfContact
 * \brief Probe class for contact (with friction) with an implicit surface
 */
class ImplSurfContact	:	public Probe
{
public:
    ImplSurfContact( sofa::component::container::ImplicitSurface *ImplicitSurface = NULL);
    ~ImplSurfContact(){}

    void setImplicitSurface(sofa::component::container::ImplicitSurface *s) { _implSurf = s; }

    inline bool inContact(){return _inContact;}
    void noContact();

    void getOrthogonalVectors(const defaulttype::Vec3d& dir, defaulttype::Vec3d& vec1, defaulttype::Vec3d& vec2);
    void rotateFrame(const defaulttype::Vec3d& n_old, const defaulttype::Vec3d& n_new, defaulttype::Vec3d& t, defaulttype::Vec3d& s);
    void changeContactFrame(defaulttype::Vec3d& posSurf);


    void initStep(){_posBegin = _posSurf;}

    bool detectNewContact();
    //bool solveConstraint();

    void storeBeginInfo();
    bool solveConstraint2(double &mu);	// new procedure

    void storeEndInfo();
    void projectSurfacePos();

    void draw();

    /// Output stream
    inline friend std::ostream& operator<< ( std::ostream& os, const ImplSurfContact& /*isc*/ )
    {
        return os;
    }

    /// Input stream
    inline friend std::istream& operator>> ( std::istream& in, ImplSurfContact& /*isc*/ )
    {
            return in;
    }

    int _domain;
private:
    sofa::component::container::ImplicitSurface *_implSurf;
    defaulttype::Vec3d	_posSurf; // the ref position on the surface
    defaulttype::Vec3d  _posBegin; // the ref position on the surface at the beginning of the time step (useful for friction)
    defaulttype::Vec3d	_n, _t, _s; // contact normal and 2 tg
    double _fn, _ft, _fs;
    double _wnn, _wtt, _wss;
    bool _inContact;
    bool _friction;

};





/*!
 * \class ImplicitSurfaceAdaptiveConstraint
 * \brief ImplicitSurfaceAdaptiveConstraint Class
 */
template<class DataTypes>
class ImplicitSurfaceAdaptiveConstraint : public core::behavior::BaseInteractionConstraint
{
public:        
    SOFA_CLASS(SOFA_TEMPLATE(ImplicitSurfaceAdaptiveConstraint,DataTypes),sofa::core::behavior::BaseInteractionConstraint);

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::CPos CPos;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::DPos DPos;
    typedef typename Coord::value_type Real;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;

    typedef core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;




protected:
    MechanicalState* object1; ///< MechanicalState that drives the implicit surface : for now, not used//
    MechanicalState* object2; ///< MechanicalState of the points constrained by the implicit surface //
    bool yetIntegrated;

    unsigned int cid;
    Data<double> mu;

    // For visualization of the potential //
    Data<bool> visualization;
    Data<defaulttype::Vec3d> PosMin;
    Data<defaulttype::Vec3d> PosMax;

    // domain
    Data<int> init_domain;

    // all points are active or not
    Data<Real> cancelPointThreshold;
    sofa::helper::vector<bool> active_list;

    // pointer to the implicit surface //
    sofa::component::container::ImplicitSurface* _contact_surface;

    //void getOrthogonalVectors(const Deriv& dir, Deriv& vec1, Deriv& vec2);

    /// this constraint is NOT holonomic
    bool isHolonomic() {return false;}
    int _nbConstraints;

public:

    ImplicitSurfaceAdaptiveConstraint(MechanicalState* object1, MechanicalState* object2)
    : object1(object1), object2(object2), yetIntegrated(false)
    , mu(initData(&mu, 0.0, "friction", "friction coefficient"))
    , visualization(initData(&visualization, false, "visualization", "visualization of the implicit surface potential"))
    , PosMin(initData(&PosMin, defaulttype::Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , PosMax(initData(&PosMax, defaulttype::Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , init_domain(initData(&init_domain, 0, "domainId", "optional1: give an initial id for the domain"))
    , cancelPointThreshold(initData(&cancelPointThreshold, (Real) 0.0, "cancelPointThreshold", "the point for which the potential is more than the threshold are desactivated"))
    , listProbe(initData(&listProbe, "listProbe", "listProbe"))
    {
    }

    ImplicitSurfaceAdaptiveConstraint(MechanicalState* object)
    : object1(object), object2(object), yetIntegrated(false)
    , mu(initData(&mu, 0.0, "friction", "friction coefficient"))
    , visualization(initData(&visualization, false, "visualization", "visualization of the implicit surface potential"))
    , PosMin(initData(&PosMin, defaulttype::Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , PosMax(initData(&PosMax, defaulttype::Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , init_domain(initData(&init_domain, 0, "domainId", "optional2: give an initial id for the domain"))
    , cancelPointThreshold(initData(&cancelPointThreshold, (Real) 0.0, "cancelPointThreshold", "the point for which the potential is more than the threshold are desactivated"))
    , listProbe(initData(&listProbe, "listProbe", "listProbe"))
    {
    }

    ImplicitSurfaceAdaptiveConstraint()
    : object1(NULL), object2(NULL), yetIntegrated(false)
    , mu(initData(&mu, 0.0, "friction", "friction coefficient"))
    , visualization(initData(&visualization, false, "visualization", "visualization of the implicit surface potential"))
    , PosMin(initData(&PosMin, defaulttype::Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , PosMax(initData(&PosMax, defaulttype::Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , init_domain(initData(&init_domain, 0, "domainId", "optional3: give an initial id for the domain"))
    , cancelPointThreshold(initData(&cancelPointThreshold, (Real) 0.0, "cancelPointThreshold", "the point for which the potential is more than the threshold are desactivated"))
    , listProbe(initData(&listProbe, "listProbe", "listProbe"))
    {
    }

     ~ImplicitSurfaceAdaptiveConstraint()
    {
    }

    MechanicalState* getObject1() { return object1; }
    MechanicalState* getObject2() { return object2; }
    core::behavior::BaseMechanicalState* getMechModel1() { return object1; }
    core::behavior::BaseMechanicalState* getMechModel2() { return object2; }



    /////////////////// to be modified ////////////
    // implicit function of a sphere
    double  getValue(Coord& Pos);
    defaulttype::Vec3d getGradient(Coord& Pos) ;

    Coord _posSphere;
    ///////////////////////////////////////////////
    sofa::component::topology::PointData<sofa::helper::vector< ImplSurfContact> > listProbe;

    //std::vector<ImplSurfContact*> listProbe ;

    void init();

    void clear();

    void internalInit();

    void bwdInit();

    void reset();

    // void handleTopologyChange();

    //void applyConstraint(unsigned int &constraintId);

    void getConstraintValue(double* v /*, unsigned int &numContacts */);

    void getConstraintId(long* , unsigned int &) { } // not used //

    void getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset);

    void projectJacobianMatrix() { }

    void buildConstraintMatrix(unsigned int& constraintId, sofa::core::VecId) {
        std::cout<<" THIS << buildConstraintMatrix >> SHOULD NOT BE CALLED "<<std::endl;
        //applyConstraint(constraintId);
    }


    /////////////// TODO: use the standart interface in the constraint //////////////

     void buildConstraintMatrix(const core::ConstraintParams* cParams /* PARAMS FIRST =ConstraintParams::defaultInstance()*/, core::MultiMatrixDerivId cId, unsigned int &constraintId);


     void getConstraintViolation(const sofa::core::ConstraintParams * /* PARAMS FIRST */, sofa::defaulttype::BaseVector *v) {}


    // Previous Constraint Interface
    void projectResponse(){}
    void projectVelocity(){}
    void projectPosition(){}
    void projectFreeVelocity(){}
    void projectFreePosition(){}

	//TODO : BaseInteractionConstraint
//	void buildConstraintMatrix(const sofa::core::ConstraintParams * /* PARAMS FIRST */, sofa::core::MultiMatrixDerivId,unsigned int &) {}
//	void getConstraintViolation(const sofa::core::ConstraintParams * /* PARAMS FIRST */, sofa::defaulttype::BaseVector *) {}

    /// Pre-construction check method called by ObjectFactory.
    template<class T>
    static bool canCreate(T*& obj, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        if (arg->getAttribute("object1") || arg->getAttribute("object2"))
        {
        std::cout << "TESTING object1 " << arg->getAttribute("object1","..") << std::endl;
            if (dynamic_cast<MechanicalState*>(arg->findObject(arg->getAttribute("object1",".."))) == NULL)
                return false;
        std::cout << "TESTING object2 " << arg->getAttribute("object2","..") << std::endl;
            if (dynamic_cast<MechanicalState*>(arg->findObject(arg->getAttribute("object2",".."))) == NULL)
                return false;
        std::cout << "TESTING DONE." << std::endl;
        }
        else
        {
            if (dynamic_cast<MechanicalState*>(context->getMechanicalState()) == NULL)
                return false;
        }
        return core::behavior::BaseInteractionConstraint::canCreate(obj, context, arg);
    }

    /// Construction method called by ObjectFactory.
    template<class T>
    static typename T::SPtr create(T* p0, core::objectmodel::BaseContext* context, core::objectmodel::BaseObjectDescription* arg)
    {
        typename T::SPtr obj = core::behavior::BaseInteractionConstraint::create(p0, context, arg);
        if (arg && (arg->getAttribute("object1") || arg->getAttribute("object2")))
        {
            obj->object1 = dynamic_cast<MechanicalState*>(arg->findObject(arg->getAttribute("object1","..")));
            obj->object2 = dynamic_cast<MechanicalState*>(arg->findObject(arg->getAttribute("object2","..")));
        }
        else if (context)
        {
            obj->object1 =
            obj->object2 =
                dynamic_cast<MechanicalState*>(context->getMechanicalState());
        }
        return obj;
    }

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const ImplicitSurfaceAdaptiveConstraint<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    void draw(const core::visual::VisualParams* vparams);
};


/////////////////////////////////////////////////////////
/*!
 * \class ImplicitSurfaceAdaptiveConstraintResolution
 * \brief ImplicitSurfaceAdaptiveConstraintResolution Class
 */
template<class DataTypes>
class ImplicitSurfaceAdaptiveConstraintResolution : public core::behavior::ConstraintResolution
{
public:
    ImplicitSurfaceAdaptiveConstraintResolution(double frictionCoef, ImplicitSurfaceAdaptiveConstraint<DataTypes>* constraintInput, int id_probe)
    : _implSurfConstraint(constraintInput), _mu(frictionCoef), _first(true), _id_probe(id_probe) { nbLines = 3; }
    virtual void resolution(int line, double** w, double* d, double* force);


private:
    ImplicitSurfaceAdaptiveConstraint<DataTypes>* _implSurfConstraint;
    double _mu;  /*, _fConstraint*/
    bool _first;
    int _id_probe;
};


} // namespace constraint

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINT_BILATERALINTERACTIONCONSTRAINT_H
