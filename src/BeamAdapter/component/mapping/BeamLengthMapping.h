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
// C++ Implementation : AdaptiveBeamMapping
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_MAPPING_BEAMLENGTHMAPPING_H
#define SOFA_COMPONENT_MAPPING_BEAMLENGTHMAPPING_H

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/type/vector.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/type/Vec.h>
#include <sofa/type/Mat.h>
#include <sofa/core/Mapping.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/visual/VisualParams.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyModifier.h>

#include <BeamAdapter/config.h>
#include <BeamAdapter/component/BeamInterpolation.h>
#include <BeamAdapter/component/controller/AdaptiveBeamController.h>

#include <sofa/linearalgebra/EigenSparseMatrix.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Forward declarations, see https://en.wikipedia.org/wiki/Forward_declaration
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////


namespace sofa::component::mapping
{

/////////////////////////////////// private namespace pattern //////////////////////////////////////
/// To avoid the lacking of names imported with with 'using' in the other's component namespace
/// you should use a private namespace and "export" only this one in the public namespace.
/// This is done at the end of this file, have a look if you are not used to this pattern.
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace _beamlengthmapping_
{

using namespace sofa::component::fem;
using namespace sofa::core::objectmodel;

using sofa::core::State ;
using core::Mapping;
using sofa::type::Vec;
using sofa::type::Mat;
using sofa::core::topology::BaseMeshTopology;
using defaulttype::SolidTypes;
using std::pair;
using sofa::component::fem::BeamInterpolation;
using sofa::type::vector;
using std::string;
using core::MechanicalParams;
using core::ConstraintParams;
using core::visual::VisualParams;
using sofa::core::topology::TopologyContainer ;

/*!
 * \class BeamLengthMapping
 * @brief Computes and map the length of the beams
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */
template <class TIn, class TOut>
class BeamLengthMapping : public Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(BeamLengthMapping,TIn,TOut),
               SOFA_TEMPLATE2(Mapping,TIn,TOut));

    typedef Mapping<TIn, TOut> Inherit;
    typedef TIn In;
    typedef TOut Out;

    typedef typename Out::Coord::value_type Real          ;
    typedef typename Out::Coord             Coord         ;
    typedef typename Out::Deriv             Deriv         ;
    typedef typename Out::VecCoord          VecCoord      ;
    typedef typename Out::VecDeriv          VecDeriv      ;
    typedef typename Out::MatrixDeriv       OutMatrixDeriv;

    typedef typename In::Coord::value_type  InReal       ;
    typedef typename In::Deriv              InDeriv      ;
    typedef typename In::VecCoord           InVecCoord   ;
    typedef typename In::VecDeriv           InVecDeriv   ;
    typedef typename In::MatrixDeriv        InMatrixDeriv;
    enum {Nin = In::deriv_total_size, Nout = Out::deriv_total_size };

    typedef vector<unsigned int>             VecIndex;
    typedef BaseMeshTopology::EdgeID         ElementID;
    typedef vector<BaseMeshTopology::Edge>   VecEdges    ;
    typedef vector<BaseMeshTopology::EdgeID> VecElementID;

    typedef typename SolidTypes<InReal>::Transform      Transform       ;
    typedef pair<int, Transform>                        IndexedTransform;
    typedef typename SolidTypes< InReal>::SpatialVector SpatialVector   ;
    typedef linearalgebra::EigenSparseMatrix<TIn,TOut>   SparseMatrixEigen;
    typedef linearalgebra::EigenSparseMatrix<TIn,TIn>    SparseKMatrixEigen;

    typedef Vec<3, Real>   Vec3;
    typedef Vec<6, Real>   Vec6;
    typedef Mat<3,3,Real>  Mat3;
    typedef Mat<3,6,Real>  Mat3x6;
    typedef Mat<3,12,Real> Mat3x12;
    typedef Mat<12,3,Real> Mat12x3;
    typedef Mat<6,12,Real> Mat6x12;
    typedef Mat<12,6,Real> Mat12x6;
    typedef BeamInterpolation<TIn> BInterpolation;



public:


    SingleLink<BeamLengthMapping<TIn, TOut>,
               BInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_adaptativebeamInterpolation;

    Data< unsigned >       d_geometricStiffness; ///< how to compute geometric stiffness (0->no GS, 1->exact GS, 2->stabilized GS)


    BeamLengthMapping(State< In >* from=NULL,
                        State< Out >* to=NULL,
                        BeamInterpolation< TIn >* interpolation=NULL) ;


    void setInterpolation(BeamInterpolation< TIn >* interpolation)
    {
        l_adaptativebeamInterpolation.set(interpolation);
    }

    virtual ~BeamLengthMapping(){}


    virtual void init() override;     // get the interpolation
    virtual void bwdInit() override;  // get the points
    virtual void reset() override;
    virtual void reinit() override;
    virtual void draw(const VisualParams*) override;

    // interface of mapping.h
    virtual void apply(const MechanicalParams *mparams, Data<VecCoord>& out, const Data<InVecCoord>& in) override;
    virtual void applyJ(const MechanicalParams *mparams, Data<VecDeriv>& out, const Data<InVecDeriv>& in) override;
    virtual void applyJT(const MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<VecDeriv>& in) override;
    virtual void applyJT(const ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in) override;
    virtual void applyDJT(const MechanicalParams* mparams, core::MultiVecDerivId parentDfId, core::ConstMultiVecDerivId childDfId) override;

    // interface of baseMapping.h
    virtual void updateK( const MechanicalParams* /*mparams*/, core::ConstMultiVecDerivId /*outForce*/ ) override;
    const linearalgebra::BaseMatrix* getK() override;





    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using Mapping<TIn, TOut>::toModel ;
    using Mapping<TIn, TOut>::fromModel ;
    ////////////////////////////////////////////////////////////////////////////

 protected:
   /* Vec3 F0_buf, F1_buf, F2_buf, F3_buf; // Used for debug */
    SparseKMatrixEigen K_geom;

    // used for applyJ on one beam
    void computeJSpline(Real &dlength, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3& P3,
                                       const Vec3& dP0, const Vec3& dP1, const Vec3& dP2, const Vec3& dP3);

    // used for applyJt on one beam
    void computeJtSpline(const Real &f_input, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3& P3,
                                                     Vec3& F0,  Vec3& F1,  Vec3& F2, Vec3& F3);

    // compute stiffness of forces for applyJt on one beam
    void computeDJtSpline(const Real &f_input, const Vec3& P0, const Vec3& P1, const Vec3& P2, const Vec3& P3,
                                                        Mat<4,4,Mat3> &Mat);


    // useful: create a cross matrix (to produce the cross product)
    void createCrossMatrix(const Vec3& v, Mat3& result){ ; result.clear();
                                     result[0][1]=-v[2]; result[0][2]=v[1];
                                     result[1][0]=v[2]; result[1][2]=-v[0];
                                     result[2][0]=-v[1]; result[2][1]=v[0];}



};

#ifndef BEAMADAPTER_BEAMLENGTHMAPPING_CPP
extern template class SOFA_BEAMADAPTER_API BeamLengthMapping<defaulttype::Rigid3dTypes, defaulttype::Vec1dTypes   >;
#endif

} /// _beamlengthmapping_

using _beamlengthmapping_::BeamLengthMapping ;



} // namespace sofa::component::mapping
#endif  /* SOFA_COMPONENT_MAPPING_BEAMLENGTHMAPPING_H */
