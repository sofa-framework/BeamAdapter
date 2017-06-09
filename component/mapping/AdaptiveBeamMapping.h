/******************************************************************************
*       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
*                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
*                                                                             *
* This library is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This library is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this library; if not, write to the Free Software Foundation,     *
* Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA.          *
*******************************************************************************
*                               SOFA :: Modules                               *
*                                                                             *
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
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
#ifndef SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_H
#define SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_H

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <sofa/core/Mapping.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/visual/VisualParams.h>
#include <SofaBaseTopology/EdgeSetTopologyModifier.h>

#include "../initBeamAdapter.h"
#include "../BeamInterpolation.h"
#include "../controller/AdaptiveBeamController.h"


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Forward declarations, see https://en.wikipedia.org/wiki/Forward_declaration
////////////////////////////////////////////////////////////////////////////////////////////////////



////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa
{

namespace component
{

namespace mapping
{

/////////////////////////////////// private namespace pattern //////////////////////////////////////
/// To avoid the lacking of names imported with with 'using' in the other's component namespace
/// you should use a private namespace and "export" only this one in the public namespace.
/// This is done at the end of this file, have a look if you are not used to this pattern.
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace _adaptivebeammapping_
{

using namespace sofa::component::fem;
using namespace sofa::core::objectmodel;

using sofa::core::State ;
using core::Mapping;
using sofa::defaulttype::Vec;
using sofa::defaulttype::Mat;
using sofa::core::topology::BaseMeshTopology;
using defaulttype::SolidTypes;
using std::pair;
using sofa::component::fem::BeamInterpolation;
using sofa::helper::vector;
using std::string;
using core::MechanicalParams;
using core::ConstraintParams;
using core::visual::VisualParams;

/*!
 * \class AdaptiveBeamMapping
 * @brief Set the positions and velocities of points attached to a beam using linear interpolation between DOFs
 *
 * This is a component:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 */
template <class TIn, class TOut>
class AdaptiveBeamMapping : public Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(AdaptiveBeamMapping,TIn,TOut),
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

    typedef vector<unsigned int>             VecIndex;
    typedef BaseMeshTopology::EdgeID         ElementID;
    typedef vector<BaseMeshTopology::Edge>   VecEdges    ;
    typedef vector<BaseMeshTopology::EdgeID> VecElementID;

    typedef typename SolidTypes<InReal>::Transform      Transform       ;
    typedef pair<int, Transform>                        IndexedTransform;
    typedef typename SolidTypes< InReal>::SpatialVector SpatialVector   ;

    typedef Vec<3, Real>   Vec3;
    typedef Vec<6, Real>   Vec6;
    typedef Mat<3,12,Real> Mat3x12;
    typedef Mat<12,3,Real> Mat12x3;
    typedef Mat<6,12,Real> Mat6x12;
    typedef Mat<12,6,Real> Mat12x6;
    typedef BeamInterpolation<TIn> BInterpolation;

    typedef pair<unsigned int, Vec3> BeamIdAndBaryCoord;
    typedef struct
    {
       unsigned int beamId;
       /// A bary point has 3 components
       /// -The first denote the curvilinear coordinate
       /// -The two followings denote the planar coordinate on the perpendicular cross section on the curve
       Vec3 baryPoint;
    } PosPointDefinition;

public:

    Data<bool> d_useCurvAbs;			  /*!< true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same */
    Data<vector<Vec3>> d_points;	      /*!< defines the mapped points along the beam axis (in beam frame local coordinates) */
    Data<double> d_proximity;			  /*!< if positive, the mapping is modified for the constraints to take into account the lever created by the proximity */
    Data<bool> d_contactDuplicate;		  /*!< if true, this mapping is a copy of an input mapping and is used to gather contact points (ContinuousFrictionContact Response) */
    Data<string> d_inputMapName;		  /*!< if contactDuplicate==true, it provides the name of the input mapping */
    Data<double> d_nbPointsPerBeam;		  /*!< if non zero, we will adapt the points depending on the discretization, with this num of points per beam (compatible with useCurvAbs)*/
    Data<vector<Real>> d_segmentsCurvAbs; /*!< (output) the abscissa of each created point on the collision model */

    SingleLink<AdaptiveBeamMapping<TIn, TOut>,
               BInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_adaptativebeamInterpolation;


    AdaptiveBeamMapping(State< In >* from=NULL,
                        State< Out >* to=NULL,
                        BeamInterpolation< TIn >* interpolation=NULL,
                        bool isSubMapping=false) ;

    virtual ~AdaptiveBeamMapping(){}


    virtual void init() override;     // get the interpolation
    virtual void bwdInit() override;  // get the points
    virtual void reset() override;
    virtual void reinit() override;
    virtual void draw(const VisualParams*) override;


    virtual void apply(const MechanicalParams *mparams, Data<VecCoord>& out, const Data<InVecCoord>& in) override;
    virtual void applyJ(const MechanicalParams *mparams, Data<VecDeriv>& out, const Data<InVecDeriv>& in) override;
    virtual void applyJT(const MechanicalParams *mparams, Data<InVecDeriv>& out, const Data<VecDeriv>& in) override;
    virtual void applyJT(const ConstraintParams *cparams, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in) override;


    void printIstrumentInfo()const;
    int addPoint(const Coord& c, int indexFrom);
    int addContactPoint(const Vec3& bary);
    void setBarycentricMapping();
    int addBaryPoint(const int& beamId, const Vec3& baryCoord, bool straightlineSplineOption) ;

    //clear the mapping in functions of size given
    void clear(int size) ;
    void computeIdxAndBaryCoordsForAbs(unsigned int &b, Real &xBary, const Real &xAbs);

    void beginAddContactPoint();
    void clearIdPointSubMap() {m_idPointSubMap.clear();}
    void addIdPointSubMap(unsigned int id) {m_idPointSubMap.push_back(id);}
    void setUseCurvAbs(bool value) {d_useCurvAbs.setValue(value);}

    const vector<PosPointDefinition>& getPointBeamDistribution() const {return m_pointBeamDistribution;}

protected:

    ////////////////////////// Inherited attributes ////////////////////////////
    /// https://gcc.gnu.org/onlinedocs/gcc/Name-lookup.html
    /// Bring inherited attributes and function in the current lookup context.
    /// otherwise any access to the base::attribute would require
    /// the "this->" approach.
    using Mapping<TIn, TOut>::toModel ;
    using Mapping<TIn, TOut>::fromModel ;
    ////////////////////////////////////////////////////////////////////////////

    TopologyContainer* m_topology;

    bool m_isXBufferUsed;
    typename In::VecCoord m_xBuffer;

    vector< PosPointDefinition > m_pointBeamDistribution;

    /// for continuous_friction_contact:
    AdaptiveBeamMapping<TIn, TOut>* m_inputMapping;
    vector<unsigned int> m_idPointSubMap;
    bool m_isSubMapping ;
    bool m_isBarycentricMapping;


    void applyJonPoint(unsigned int i, SpatialVector& VNode0input, SpatialVector& VNode1input, Deriv& vOutput, const InVecCoord& x);
    void applyJTonPoint(unsigned int i, const Deriv& finput, SpatialVector& FNode0output, SpatialVector& FNode1output, const InVecCoord& x);
    void computeJacobianOnPoint(unsigned int i, const typename In::VecCoord& x);
    void computeDistribution();
};


} /// _adaptivebeammappin_

using _adaptivebeammapping_::AdaptiveBeamMapping ;

} /// namespace mapping

} /// namespace component

} /// namespace sofa

#endif  /* SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_H */
