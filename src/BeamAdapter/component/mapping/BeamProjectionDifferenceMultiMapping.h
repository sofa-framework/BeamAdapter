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
//
// Author: Eulalie Coevoet
//
// Copyright: See COPYING file that comes with this distribution

#pragma once

#include <sofa/core/BaseMapping.h>
#include <sofa/core/config.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/Multi2Mapping.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/linearalgebra/EigenSparseMatrix.h>
#include <BeamAdapter/component/BeamInterpolation.h>

#include <BeamAdapter/config.h>


namespace beamadapter::mapping
{
using sofa::defaulttype::SolidTypes ;
using sofa::type::Matrix3;
using sofa::type::Matrix4;
using sofa::type::Vec3;
using sofa::type::Vec6;
using std::get;
using sofa::type::vector;

/**
* \class BeamProjectionDifferenceMultiMapping
* @brief Computes the difference between a model's points and their projection on a beam
*/


template <class TIn1, class TIn2, class TOut>
class SOFA_BEAMADAPTER_API BeamProjectionDifferenceMultiMapping : public sofa::core::Multi2Mapping<TIn1, TIn2, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE3(BeamProjectionDifferenceMultiMapping, TIn1,TIn2, TOut), SOFA_TEMPLATE3(sofa::core::Multi2Mapping, TIn1, TIn2, TOut) );
    typedef sofa::core::Multi2Mapping<TIn1, TIn2, TOut> Inherit;

    /// Input Model Type
    typedef TIn1 In1;
    typedef TIn2 In2;

    /// Output Model Type
    typedef TOut Out;

    typedef typename In1::Coord In1Coord;
    typedef typename In1::Deriv In1Deriv;
    typedef typename In1::VecCoord In1VecCoord;
    typedef typename In1::VecDeriv In1VecDeriv;
    typedef typename In1::MatrixDeriv In1MatrixDeriv;
    typedef sofa::Data<In1VecCoord> In1DataVecCoord;
    typedef sofa::Data<In1VecDeriv> In1DataVecDeriv;
    typedef sofa::Data<In1MatrixDeriv> In1DataMatrixDeriv;

    typedef sofa::defaulttype::Rigid3dTypes::Coord Rigid;
    
    typedef typename In2::Coord::value_type Real;
    typedef typename In2::Coord             In2Coord;
    typedef typename In2::Deriv             In2Deriv;
    typedef typename In2::VecCoord In2VecCoord;
    typedef typename In2::VecDeriv In2VecDeriv;
    typedef typename In2::MatrixDeriv In2MatrixDeriv;
    typedef sofa::Data<In2VecCoord> In2DataVecCoord;
    typedef sofa::Data<In2VecDeriv> In2DataVecDeriv;
    typedef sofa::Data<In2MatrixDeriv> In2DataMatrixDeriv;
    typedef sofa::type::Mat<4,4,Real> Mat4x4;

    typedef typename Out::VecCoord OutVecCoord;
    typedef typename Out::Coord OutCoord;
    typedef typename Out::Deriv OutDeriv;
    typedef typename Out::VecDeriv OutVecDeriv;
    typedef typename Out::MatrixDeriv OutMatrixDeriv;
    typedef sofa::Data<OutVecCoord> OutDataVecCoord;
    typedef sofa::Data<OutVecDeriv> OutDataVecDeriv;
    typedef sofa::Data<OutMatrixDeriv> OutDataMatrixDeriv;

    typedef typename SolidTypes<Real>::Transform Transform;

    typedef sofa::type::vector<sofa::topology::Edge> SeqEdges;

    enum
    {
        N = Out::spatial_dimensions
    };
    enum
    {
        NIn1 = sofa::defaulttype::DataTypeInfo<In1Deriv>::Size
    };
    enum
    {
        NIn2 = sofa::defaulttype::DataTypeInfo<In2Deriv>::Size
    };
    enum
    {
        NOut = sofa::defaulttype::DataTypeInfo<OutDeriv>::Size
    };

protected:
    BeamProjectionDifferenceMultiMapping() ;
    ~BeamProjectionDifferenceMultiMapping()  override {}

public:
    void init() override;
    void draw(const sofa::core::visual::VisualParams* vparams) override;

    void apply(
            const sofa::core::MechanicalParams* mparams, const sofa::type::vector<OutDataVecCoord*>& dataVecOutPos,
            const sofa::type::vector<const In1DataVecCoord*>& dataVecIn1Pos ,
            const sofa::type::vector<const In2DataVecCoord*>& dataVecIn2Pos) override;

    void applyJ(
            const sofa::core::MechanicalParams* mparams, const sofa::type::vector< OutDataVecDeriv*>& dataVecOutVel,
            const sofa::type::vector<const In1DataVecDeriv*>& dataVecIn1Vel,
            const sofa::type::vector<const In2DataVecDeriv*>& dataVecIn2Vel) override;

    void applyJT(
            const sofa::core::MechanicalParams* mparams, const sofa::type::vector< In1DataVecDeriv*>& dataVecOut1Force,
            const sofa::type::vector< In2DataVecDeriv*>& dataVecOut2RootForce,
            const sofa::type::vector<const OutDataVecDeriv*>& dataVecInForce) override;

    void applyDJT(const sofa::core::MechanicalParams* mparams,
                  sofa::core::MultiVecDerivId inForce,
                  sofa::core::ConstMultiVecDerivId outForce) override;

    virtual void applyJT(
            const sofa::core::ConstraintParams*  cparams , const sofa::type::vector< In1DataMatrixDeriv*>& dataMatOut1Const  ,
            const sofa::type::vector< In2DataMatrixDeriv*>&  dataMatOut2Const ,
            const sofa::type::vector<const OutDataMatrixDeriv*>&  dataMatInConst) override;

    virtual const sofa::type::vector<sofa::linearalgebra::BaseMatrix*>* getJs() override;

    void computeProjection(const In1VecCoord &xFrom, const In2VecCoord &xTo, const bool &updateOrientation);

public:
    sofa::Data<vector<unsigned int>> d_indices;
    sofa::Data<sofa::type::vector<bool>> d_directions;
    sofa::Data<bool> d_updateProjectionPosition;
    sofa::Data<bool> d_updateProjectionOrientation;
    sofa::Data<bool> d_draw;
    sofa::Data<Real> d_drawSize;

    sofa::SingleLink<BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>, sofa::core::topology::BaseMeshTopology, sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_in2Topology;
    sofa::SingleLink<BeamProjectionDifferenceMultiMapping<TIn1, TIn2, TOut>, sofa::component::fem::BeamInterpolation<TIn2>, sofa::BaseLink::FLAG_STOREPATH | sofa::BaseLink::FLAG_STRONGLINK> l_interpolation;

    using sofa::core::Multi2Mapping<TIn1, TIn2, TOut>::d_componentState ;

protected:
    sofa::core::State<In1>* m_fromModel1;
    sofa::core::State<In2>* m_fromModel2;
    sofa::core::State<Out>* m_toModel;

    typedef sofa::linearalgebra::EigenSparseMatrix<In1,Out> SparseMatrixEigen1;
    SparseMatrixEigen1 m_eigenJacobian1;
    typedef sofa::linearalgebra::EigenSparseMatrix<In1,Out> SparseMatrixEigen2;
    SparseMatrixEigen2 m_eigenJacobian2;                      ///< Jacobian of the mapping used by getJs
    sofa::type::vector<sofa::linearalgebra::BaseMatrix*> m_eigenJacobians; /// used by getJs

    bool m_updateJ;

private:

    typedef struct {
        int edgeIndex;
        int pi1, pi2;  // edge's points index
        Transform interpolatedTransform;
        Real alpha;
        bool onBeam;
    } MappedPoint;

    sofa::type::vector<MappedPoint> m_mappedPoints;

};

// Declares template as extern to avoid the code generation of the template for
// each compilation unit. see: http://www.stroustrup.com/C++11FAQ.html#extern-templates
#if !defined(BEAMADAPTER_MAPPING_BEAMPROJECTIONDIFFERENCEMULTIMAPPING_CPP)
//extern template class BeamProjectionDifferenceMultiMapping< sofa::defaulttype::Vec3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Vec3Types >;
extern template class SOFA_BEAMADAPTER_API BeamProjectionDifferenceMultiMapping< sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types, sofa::defaulttype::Rigid3Types >;
#endif

} // namespace
