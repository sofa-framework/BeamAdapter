#ifndef SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_H

#include <sofa/core/behavior/BaseInteractionConstraint.h>
#include <sofa/core/behavior/PairInteractionConstraint.h>
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

#include "WireBeamInterpolation.h"


#include <../applications/plugins/SofaEVE/Implicit/Isosurface.h>


namespace sofa
{

namespace component
{

namespace constraint
{
using namespace sofa::defaulttype;


/////////////////////////////////////////////////////////
/*!
 * \class ImplicitSurfaceAdaptiveConstraintResolution
 * \brief ImplicitSurfaceAdaptiveConstraintResolution Class
 */
template<class DataTypes>
class ImplicitSurfaceAdaptiveConstraintResolution : public core::behavior::ConstraintResolution
{
public:
    ImplicitSurfaceAdaptiveConstraintResolution(double frictionCoef, int line, sofa::component::fem::WireBeamInterpolation<DataTypes>* wireInterpol)
    : _wireInterpolation(wireInterpol), _mu(frictionCoef), _line(line) { nbLines = 3; }
    virtual void resolution(int line, double** w, double* d, double* force);


private:
    sofa::component::fem::WireBeamInterpolation<DataTypes>* _wireInterpolation;
    double _mu;  /*, _fConstraint*/
    int _line;
};


/*!
 * \class ImplicitSurfaceAdaptiveConstraint
 * \brief ImplicitSurfaceAdaptiveConstraint Class
 */
template<class DataTypes>
class ImplicitSurfaceAdaptiveConstraint : public core::behavior::PairInteractionConstraint<DataTypes>
{
public:        
    SOFA_CLASS(SOFA_TEMPLATE(ImplicitSurfaceAdaptiveConstraint,DataTypes),SOFA_TEMPLATE(sofa::core::behavior::PairInteractionConstraint,DataTypes));


    typedef typename core::behavior::PairInteractionConstraint<DataTypes> Inherit;

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::MatrixDeriv MatrixDeriv;
    typedef typename MatrixDeriv::RowIterator MatrixDerivRowIterator;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::CPos CPos;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::DPos DPos;
    typedef typename Coord::value_type Real;
    typedef  Vec<3,Real> Vec3;
    typedef Vec<3,double> Vec3d;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef  sofa::component::fem::WireBeamInterpolation<DataTypes> WBInterpolation;

    typedef core::objectmodel::Data<VecCoord>		DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv>		DataVecDeriv;
    typedef core::objectmodel::Data<MatrixDeriv>    DataMatrixDeriv;

    typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
    typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

    typedef typename sofa::helper::vector<int>::iterator vectorIntIterator;

protected :
    /// pointer to the interpolation
    SingleLink<ImplicitSurfaceAdaptiveConstraint<DataTypes>, WBInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_wireBinterpolation;



    unsigned int cid;

    // For visualization of the potential //
    Data<bool> visualization;
    Data<defaulttype::Vec3d> PosMin;
    Data<defaulttype::Vec3d> PosMax;
    sofaeve::implicit::MarchingCube * mc;
    double mc_data[100000000];
    Real _isoValue;

    // domain
    Data<int> init_domain;

    // all points are active or not
    Data<Real> alarmDistance;
    Data<Real> radius;
    sofa::helper::vector<bool> active_list;

    // friction coef
    Data<Real> frictionCoef;
    bool friction;

    // list of beams to be considered for collision
    Data< sofa::helper::vector<int> > listBeams;
    bool all_activated;

    // pointer to the implicit surface //
    sofa::component::container::ImplicitSurface* _contact_surface;

    //void getOrthogonalVectors(const Deriv& dir, Deriv& vec1, Deriv& vec2);

    /// this constraint is NOT holonomic
    bool isHolonomic() {return false;}
    int _nbConstraints;

public:

    ImplicitSurfaceAdaptiveConstraint(MechanicalState* /*object1*/, MechanicalState* /*object2*/)
    : m_wireBinterpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , visualization(initData(&visualization, false, "visualization", "visualization of the implicit surface potential"))
    , PosMin(initData(&PosMin, defaulttype::Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , PosMax(initData(&PosMax, defaulttype::Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , init_domain(initData(&init_domain, 0, "domainId", "optional1: give an initial id for the domain"))
    , alarmDistance(initData(&alarmDistance, (Real) 1.0, "alarmDistance", "kind of alarm distance computed on implicit surface"))
    , radius(initData(&radius, (Real)1.0, "radius", "radius: this value is used to include the thickness in the collision response"))
    , frictionCoef(initData(&frictionCoef, (Real)0.0, "frictionCoef", "coefficient of friction (Coulomb's law)"))
    , listBeams(initData(&listBeams, "listBeams", "list of beams used by Interpolation that are activated for contact (all by default)"))
    {
        _isoValue=0.0;
    }

    ImplicitSurfaceAdaptiveConstraint(MechanicalState* /*object*/)
    :  m_wireBinterpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , visualization(initData(&visualization, false, "visualization", "visualization of the implicit surface potential"))
    , PosMin(initData(&PosMin, defaulttype::Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , PosMax(initData(&PosMax, defaulttype::Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , init_domain(initData(&init_domain, 0, "domainId", "optional2: give an initial id for the domain"))
    , alarmDistance(initData(&alarmDistance, (Real) 0.0, "alarmDistance", "the point for which the potential is more than the threshold are desactivated"))
    , radius(initData(&radius, (Real)1.0, "radius", "radius: this value is used to include the thickness in the collision response"))
    , frictionCoef(initData(&frictionCoef, (Real)0.0, "frictionCoef", "coefficient of friction (Coulomb's law)"))
    , listBeams(initData(&listBeams, "listBeams", "list of beams used by Interpolation that are activated for contact (all by default)"))
    {
        _isoValue=0.0;
    }

    ImplicitSurfaceAdaptiveConstraint()
    :  m_wireBinterpolation(initLink("interpolation","Path to the Interpolation component on scene"))
    , visualization(initData(&visualization, false, "visualization", "visualization of the implicit surface potential"))
    , PosMin(initData(&PosMin, defaulttype::Vec3d(-1.0,-1.0,-1.0), "PosMin", "position min for the visualization grid"))
    , PosMax(initData(&PosMax, defaulttype::Vec3d( 1.0, 1.0, 1.0), "PosMax", "position max for the visualization grid"))
    , init_domain(initData(&init_domain, 0, "domainId", "optional3: give an initial id for the domain"))
    , alarmDistance(initData(&alarmDistance, (Real) 0.0, "alarmDistance", "the point for which the potential is more than the threshold are desactivated"))
    , radius(initData(&radius, (Real)1.0, "radius", "radius: this value is used to include the thickness in the collision response"))
    , frictionCoef(initData(&frictionCoef, (Real)0.0, "frictionCoef", "coefficient of friction (Coulomb's law)"))
    , listBeams(initData(&listBeams, "listBeams", "list of beams used by Interpolation that are activated for contact (all by default)"))
    {
        _isoValue=0.0;
    }

     ~ImplicitSurfaceAdaptiveConstraint()
    {
    }

    /* deprecated */
    MechanicalState* getObject1() { return this->mstate1; }
    MechanicalState* getObject2() { return this->mstate2; }
    core::behavior::BaseMechanicalState* getMechModel1() { return this->mstate1; }
    core::behavior::BaseMechanicalState* getMechModel2() { return this->mstate2; }




    void init();

    void clear();

    void internalInit();

    void bwdInit();

    void reset();

    void reinit(){internalInit();}






    /////////////// standart interface of constraintset //////////////

    void buildConstraintMatrix(const sofa::core::ConstraintParams* cParams, DataMatrixDeriv &c1, DataMatrixDeriv &c2, unsigned int &cIndex
            , const DataVecCoord &, const DataVecCoord &x2) ;

    void getConstraintViolation(const sofa::core::ConstraintParams* cParams, defaulttype::BaseVector *v, const DataVecCoord &x1, const DataVecCoord &x2
            , const DataVecDeriv &v1, const DataVecDeriv &v2) ;

    void getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset);



    void draw(const core::visual::VisualParams* vparams);


private:
    sofa::helper::vector<Vec3> m_posSample;
    sofa::helper::vector<int> m_domainSample;

    struct potentialContact{
        unsigned int beamId;
        unsigned int posSampleId;
        Vec3 baryCoord;
        Vec3 n,t,s;
        Real d;
    };

    sofa::helper::vector<potentialContact> m_VecPotentialContact;



    // internal functions

    // get vec1 and vec2 orthogonal vectors from a given dir value
    void getOrthogonalVectors(const Vec3& dir, Vec3& vec1, Vec3& vec2);

    void detectPotentialContactOnImplicitSurface(const sofa::core::ConstVecCoordId &vecXId, sofa::helper::vector<int>& listBeam);

    void computeTangentialViolation(const Vec3 &Pos, const Vec3 &freePos, const Vec3 &t, const Vec3 &s,
                                    const Real& d, const Real& dfree, Real &dfree_t, Real &dfree_s );


};





} // namespace constraint

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINT_BILATERALINTERACTIONCONSTRAINT_H
