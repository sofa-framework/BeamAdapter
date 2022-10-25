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
#include <sofa/core/behavior/BaseInteractionConstraint.h>
#include <sofa/core/behavior/PairInteractionConstraint.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <sofa/type/Mat.h>
#include <sofa/type/Vec.h>
#include <sofa/type/vector.h>
#include <sofa/defaulttype/VecTypes.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/behavior/ConstraintResolution.h>

#include <SofaImplicitField/deprecated/ImplicitSurfaceContainer.h>

#include <BeamAdapter/component/WireBeamInterpolation.h>

//#define SOFAEVE  //temporary fix for compilation, this component does not work without SOFAEVE (which is however not accessible easily)

#ifdef SOFAEVE
#include <../applications/plugins/SofaEVE/Implicit/Isosurface.h>
#endif


namespace sofa::component::constraint
{

namespace _implicitsurfaceadaptiveconstraint_
{

/*!
 * \class ImplicitSurfaceAdaptiveConstraintResolution
 * \brief ImplicitSurfaceAdaptiveConstraintResolution Class
 */
template<class DataTypes>
class ImplicitSurfaceAdaptiveConstraintResolution : public sofa::core::behavior::ConstraintResolution
{
public:
    ImplicitSurfaceAdaptiveConstraintResolution(double frictionCoef, int line, sofa::component::fem::WireBeamInterpolation<DataTypes>* wireInterpol)
        : ConstraintResolution(3)
        , m_wireInterpolation(wireInterpol)
        , m_mu(frictionCoef)
        , m_line(line)
    {
    }

    virtual void resolution(int line, double** w, double* d, double* force);

private:
    sofa::component::fem::WireBeamInterpolation<DataTypes>* m_wireInterpolation;
    double m_mu;
    int m_line;
};


/*!
 * \class ImplicitSurfaceAdaptiveConstraint
 * \brief ImplicitSurfaceAdaptiveConstraint Class
 */
template<class DataTypes>
class ImplicitSurfaceAdaptiveConstraint : public sofa::core::behavior::PairInteractionConstraint<DataTypes>
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
    typedef type::Vec<3,Real> Vec3;
    typedef type::Vec<3,double> Vec3d;
    typedef typename core::behavior::MechanicalState<DataTypes> MechanicalState;
    typedef typename component::fem::WireBeamInterpolation<DataTypes> WBInterpolation;

    typedef Data<VecCoord>		 DataVecCoord;
    typedef Data<VecDeriv>		 DataVecDeriv;
    typedef Data<MatrixDeriv>    DataMatrixDeriv;

    typedef typename defaulttype::SolidTypes<Real>::Transform Transform;
    typedef typename defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;
    typedef typename type::vector<int>::iterator VectorIntIterator;

protected :
    SingleLink<ImplicitSurfaceAdaptiveConstraint<DataTypes>, WBInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_wireBinterpolation;

    Data<bool>  d_visualization;
    Data<Vec3d> d_posMin;
    Data<Vec3d> d_posMax;
    Data<int>   d_initDomain;
    Data<Real>  d_alarmDistance;
    Data<Real>  d_radius;
    Data<Real>  d_frictionCoef;
    Data<type::vector<int> > d_listBeams;

#ifdef SOFAEVE
    sofaeve::implicit::MarchingCube * mc;
#endif

    unsigned int m_cid;
    double m_cData[100000000];
    Real m_isoValue;
    type::vector<bool> m_activeList;
    bool m_friction;
    bool m_allActivated; /// list of beams to be considered for collision
    component::container::ImplicitSurfaceContainer* m_contactSurface;
    int m_nbConstraints;
    bool isHolonomic() {return false;} /// this constraint is NOT holonomic

public:

    ImplicitSurfaceAdaptiveConstraint(MechanicalState* object1, MechanicalState* object2);
    ImplicitSurfaceAdaptiveConstraint(MechanicalState* object);
    ImplicitSurfaceAdaptiveConstraint();

    virtual ~ImplicitSurfaceAdaptiveConstraint() override = default;

    core::behavior::BaseMechanicalState* getMechModel1() { return this->mstate1; }
    core::behavior::BaseMechanicalState* getMechModel2() { return this->mstate2; }

    void init();
    void clear();
    void internalInit();
    void bwdInit();
    void reset();
    void reinit(){internalInit();}

    void buildConstraintMatrix(const core::ConstraintParams* cParams,
                               DataMatrixDeriv &c1, DataMatrixDeriv &c2, unsigned int &cIndex,
                               const DataVecCoord &, const DataVecCoord &x2) ;

    void getConstraintViolation(const core::ConstraintParams* cParams, linearalgebra::BaseVector *v,
                                const DataVecCoord &x1, const DataVecCoord &x2,
                                const DataVecDeriv &v1, const DataVecDeriv &v2) ;

    void getConstraintResolution(std::vector<core::behavior::ConstraintResolution*>& resTab, unsigned int& offset);

    void draw(const core::visual::VisualParams* vparams);


private:
    type::vector<Vec3> m_posSample;
    type::vector<int> m_domainSample;

    struct potentialContact{
        unsigned int beamId;
        unsigned int posSampleId;
        Vec3 baryCoord;
        Vec3 n,t,s;
        Real d;
    };

    type::vector<potentialContact> m_vecPotentialContact;

    void getOrthogonalVectors(const Vec3& dir, Vec3& vec1, Vec3& vec2);
    void detectPotentialContactOnImplicitSurface(const core::ConstVecCoordId &vecXId, type::vector<int>&listBeam);
    void computeTangentialViolation(const Vec3 &Pos, const Vec3 &freePos, const Vec3 &t, const Vec3 &s,
                                    const Real& d, const Real& dfree, Real &dfree_t, Real &dfree_s );


    using core::behavior::PairInteractionConstraint<DataTypes>::mstate1;
    using core::behavior::PairInteractionConstraint<DataTypes>::mstate2;

};

} // namespace _implicitsurfaceadaptiveconstraint_

using _implicitsurfaceadaptiveconstraint_::ImplicitSurfaceAdaptiveConstraint;
using _implicitsurfaceadaptiveconstraint_::ImplicitSurfaceAdaptiveConstraintResolution;

} // namespace sofa::component::constraint
