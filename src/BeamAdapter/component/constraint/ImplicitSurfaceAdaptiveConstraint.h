#ifndef SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_H
#define SOFA_COMPONENT_CONSTRAINT_IMPLICITSURFACEADAPTIVECONSTRAINT_H

#include <sofa/core/behavior/BaseInteractionConstraint.h>
#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/visual/VisualModel.h>
#include <sofa/helper/gl/template.h>
#include <iostream>

#include <sofa/defaulttype/Mat.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/helper/vector.h>
#include <SofaBaseTopology/TopologyData.h>

#include <SofaImplicitField/deprecated/ImplicitSurfaceContainer.h>

#include "../WireBeamInterpolation.h"

//#define SOFAEVE  //temporary fix for compilation, this component does not work without SOFAEVE (which is however not accessible easily)

#ifdef SOFAEVE
#include <../applications/plugins/SofaEVE/Implicit/Isosurface.h>
#endif


namespace sofa
{

namespace component
{

namespace constraint
{

namespace _implicitsurfaceadaptiveconstraint_
{

using namespace sofa::defaulttype;
using sofa::core::ConstVecCoordId;
using core::behavior::ConstraintResolution;
using sofa::component::fem::WireBeamInterpolation;
using sofa::helper::vector;
using sofa::core::behavior::MechanicalState;
using defaulttype::Vec3d;
using sofa::component::container::ImplicitSurfaceContainer;
using sofa::core::behavior::PairInteractionConstraint;
using defaulttype::Vec3d;
using core::behavior::BaseMechanicalState;
using core::visual::VisualParams;
using sofa::core::ConstraintParams;


/*!
 * \class ImplicitSurfaceAdaptiveConstraintResolution
 * \brief ImplicitSurfaceAdaptiveConstraintResolution Class
 */
template<class DataTypes>
class ImplicitSurfaceAdaptiveConstraintResolution : public ConstraintResolution
{
public:
    ImplicitSurfaceAdaptiveConstraintResolution(double frictionCoef, int line, WireBeamInterpolation<DataTypes>* wireInterpol)
        : ConstraintResolution(3)
        , m_wireInterpolation(wireInterpol)
        , m_mu(frictionCoef)
        , m_line(line)
    {
    }

    virtual void resolution(int line, double** w, double* d, double* force);

private:
    WireBeamInterpolation<DataTypes>* m_wireInterpolation;
    double m_mu;
    int m_line;
};


/*!
 * \class ImplicitSurfaceAdaptiveConstraint
 * \brief ImplicitSurfaceAdaptiveConstraint Class
 */
template<class DataTypes>
class ImplicitSurfaceAdaptiveConstraint : public PairInteractionConstraint<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ImplicitSurfaceAdaptiveConstraint,DataTypes),SOFA_TEMPLATE(PairInteractionConstraint,DataTypes));

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
    typedef Vec<3,Real> Vec3;
    typedef Vec<3,double> Vec3d;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef WireBeamInterpolation<DataTypes> WBInterpolation;

    typedef Data<VecCoord>		 DataVecCoord;
    typedef Data<VecDeriv>		 DataVecDeriv;
    typedef Data<MatrixDeriv>    DataMatrixDeriv;

    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;
    typedef typename vector<int>::iterator VectorIntIterator;

protected :
    SingleLink<ImplicitSurfaceAdaptiveConstraint<DataTypes>, WBInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_wireBinterpolation;

    Data<bool>  d_visualization;
    Data<Vec3d> d_posMin;
    Data<Vec3d> d_posMax;
    Data<int>   d_initDomain;
    Data<Real>  d_alarmDistance;
    Data<Real>  d_radius;
    Data<Real>  d_frictionCoef;
    Data<vector<int> > d_listBeams;

#ifdef SOFAEVE
    sofaeve::implicit::MarchingCube * mc;
#endif

    unsigned int m_cid;
    double m_cData[100000000];
    Real m_isoValue;
    vector<bool> m_activeList;
    bool m_friction;
    bool m_allActivated; /// list of beams to be considered for collision
    ImplicitSurfaceContainer* m_contactSurface;
    int m_nbConstraints;
    bool isHolonomic() {return false;} /// this constraint is NOT holonomic

public:

    ImplicitSurfaceAdaptiveConstraint(MechanicalState* object1, MechanicalState* object2);
    ImplicitSurfaceAdaptiveConstraint(MechanicalState* object);
    ImplicitSurfaceAdaptiveConstraint();

    ~ImplicitSurfaceAdaptiveConstraint(){}

    BaseMechanicalState* getMechModel1() { return this->mstate1; }
    BaseMechanicalState* getMechModel2() { return this->mstate2; }

    void init();
    void clear();
    void internalInit();
    void bwdInit();
    void reset();
    void reinit(){internalInit();}

    void buildConstraintMatrix(const ConstraintParams* cParams,
                               DataMatrixDeriv &c1, DataMatrixDeriv &c2, unsigned int &cIndex,
                               const DataVecCoord &, const DataVecCoord &x2) ;

    void getConstraintViolation(const ConstraintParams* cParams, BaseVector *v,
                                const DataVecCoord &x1, const DataVecCoord &x2,
                                const DataVecDeriv &v1, const DataVecDeriv &v2) ;

    void getConstraintResolution(std::vector<ConstraintResolution*>& resTab, unsigned int& offset);

    void draw(const VisualParams* vparams);


private:
    vector<Vec3> m_posSample;
    vector<int> m_domainSample;

    struct potentialContact{
        unsigned int beamId;
        unsigned int posSampleId;
        Vec3 baryCoord;
        Vec3 n,t,s;
        Real d;
    };

    vector<potentialContact> m_vecPotentialContact;

    void getOrthogonalVectors(const Vec3& dir, Vec3& vec1, Vec3& vec2);
    void detectPotentialContactOnImplicitSurface(const ConstVecCoordId &vecXId, vector<int>& listBeam);
    void computeTangentialViolation(const Vec3 &Pos, const Vec3 &freePos, const Vec3 &t, const Vec3 &s,
                                    const Real& d, const Real& dfree, Real &dfree_t, Real &dfree_s );


    using PairInteractionConstraint<DataTypes>::mstate1;
    using PairInteractionConstraint<DataTypes>::mstate2;

};

}

using _implicitsurfaceadaptiveconstraint_::ImplicitSurfaceAdaptiveConstraint;
using _implicitsurfaceadaptiveconstraint_::ImplicitSurfaceAdaptiveConstraintResolution;

} // namespace constraint

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_CONSTRAINT_BILATERALINTERACTIONCONSTRAINT_H
