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
// C++ Implementation : SutureController
//
// Description:
// This controller allows to drive the discretization of a suture model computed with adaptive beam elements
// + The discretization is driven by: the curvature, some particular points on the wire (the node on the boundary between needle and thread), the suture constraints
// + It allows to detect and "rigidify" the nodes realized during suture
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
//
#pragma once

#include <sofa/helper/set.h>

#include <sofa/component/controller/MechanicalStateController.h>
#include <sofa/component/collision/geometry/PointModel.h>
#include <sofa/component/collision/geometry/LineModel.h>
#include <sofa/defaulttype/RigidTypes.h>

#include <BeamAdapter/component/WireBeamInterpolation.h>
#include <BeamAdapter/config.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalPropagateOnlyPositionAndVelocityVisitor.h>
#include <sofa/simulation/mechanicalvisitor/MechanicalProjectPositionAndVelocityVisitor.h>
#include <sofa/component/topology/container/dynamic/EdgeSetGeometryAlgorithms.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyModifier.h>
#include <sofa/defaulttype/RigidTypes.h>

////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////DECLARATIONS /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa::component::controller
{

namespace _suturecontroller_
{

using std::set ;
using sofa::type::vector ;
using sofa::type::Vec ;
using sofa::defaulttype::Rigid3dTypes ;
using sofa::defaulttype::Rigid3fTypes ;
using sofa::defaulttype::SolidTypes ;
using sofa::core::topology::TopologyContainer ;
using sofa::core::CollisionModel ;
using sofa::core::topology::BaseMeshTopology ;
using sofa::component::fem::WireBeamInterpolation;
using sofa::simulation::mechanicalvisitor::MechanicalProjectPositionAndVelocityVisitor;
using sofa::simulation::mechanicalvisitor::MechanicalPropagateOnlyPositionAndVelocityVisitor;

/**
 * \class SutureController
 * @brief SutureController Class
 *
 * Provides a Mouse & Keyboard user control on an EdgeSet Topology.
 */
template<class DataTypes>
class SutureController : public MechanicalStateController<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(SutureController, DataTypes),
               SOFA_TEMPLATE(MechanicalStateController, DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;
    typedef Vec<3, Real> Vec3;
    typedef Vec<2, Real> Vec2;

    typedef BaseMeshTopology::EdgeID ElementID;
    typedef type::vector< ElementID > VecElementID;

    typedef typename std::list< Real >::iterator ListRealIterator;

    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;
    typedef WireBeamInterpolation<DataTypes> WInterpolation;

    typedef typename set<Real>::const_iterator RealConstIterator;

    using Inherit1::getContext;
    using Inherit1::f_listening;
    using Inherit1::getMechanicalState;

public:
    SutureController(WInterpolation* _adaptiveinterpolation = nullptr) ;
    virtual ~SutureController() override = default;

    /////////////////// Inherited from Base //////////////////////////////////////////////////
    virtual void init() override;
    virtual void bwdInit() override {applyController();}
    virtual void reinit() override;
    virtual void reset() override { init(); applyController();}
    virtual void draw(const core::visual::VisualParams*) override;

    /////////////////// Inherited from Controller //////////////////////////////////////////////////
    virtual void onMouseEvent(core::objectmodel::MouseEvent *)override {}
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *) override {}
    virtual void onBeginAnimationStep(const double dt) override ;
    virtual void onEndAnimationStep(const double dt) override ;

private:
    void applyController() ;
    void insertActualNoticeablePoint(Real _absc);


    /// addNodeOnXcurv (const Real& x_curv) will add a node at abs curv "x_curv" in the list of the imposed node
    /// this list is the used in the controller for the sampling of nodes of the suture model
    void addNodeOnXcurv(const Real& x_curv){m_listOfImposedNodesOnXcurv.push_back(x_curv);}
    void clearNodesOnXcurv(){m_listOfImposedNodesOnXcurv.clear();}


    /// Checks if the controlled MechanicalState and Topology have already been initialized.
    bool wireIsAlreadyInitialized();


    /// "Wire" initialization from the starting position, rest shape, propozed discretization...
    void initWireModel();


    void recreateTopology();
    void addNodesAndEdge(unsigned int num, Real &xend);
    void removeNodesAndEdge(unsigned int num);


    /// Computes the bending angle between curv abs xmin and xmax
    /// this function calls ComputeTotalBendingRotationAngle on the beams between xmin and xmax
    Real computeBendingAngle(const Real& xmin, const Real& xmax, const Real& dx_comput, const VecCoord& Pos);


    /// Computes the tangent value on a series of discrete points (store also the curv_abs of these discrete points)
    void computeTangentOnDiscretePoints(type::vector<Vec3> TangTable, type::vector<Real> xTable,  unsigned int numDiscretePoints, const VecCoord& Pos);


    /// fill the list rigidBeamList
    void detectRigidBeams(const type::vector<Real> &newCurvAbs);


    /// when the sampling is computed, this function allows for reinterpolate the position and the velocity
    void applyNewSampling(const type::vector<Real> &newCurvAbs, const type::vector<Real> &oldCurvAbs, VecCoord &x, VecDeriv &v);


    /// Add the curv abs of the nodes at the extremity of the rigid segment
    /// if a node already exists or is very close (< tol), do not add any point
    void addRigidCurvAbs(type::vector<Real> &newCurvAbs, const Real &tol);


    /// Add the nodes that are imposed at a given curv abs
    /// if a node already exists or is very close (< tol), do not add any point
    void addImposedCurvAbs(type::vector<Real> &newCurvAbs, const Real &tol);


    /// Computes sampling
    void computeSampling(type::vector<Real> &newCurvAbs, VecCoord& x);

    /// little function to verify that the rigidCurveSegments are sorted
    bool verifyRigidCurveSegmentSort();

    /// make sure the sampling does not change for the rigid segments
    void verifyRigidSegmentsSampling(type::vector<Real> &newCurvAbs);

    void storeRigidSegmentsTransformations();
    void verifyRigidSegmentsTransformations();

    void updateControlPointsPositions();

private:
    /// for point and line activer
    type::vector<Real> m_XAbsCollisionPointsBuffer;

    /// Data:
    Data< Coord >     d_startingPos;
    Data< Real >      d_threshold;
    Data< Real >      d_maxBendingAngle;
    Data< bool >      d_useDummyController;
    Data< bool >      d_fixRigidTransforms;

    /// Values in the set are interpreted as pairs (start - end) curve absciss.
    Data< set<Real> > d_rigidCurvAbs;
    SingleLink<SutureController<DataTypes>, WInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_adaptiveInterpolation;

    /// for rigidity control
    type::vector< std::pair<Real, Real> > m_rigidCurveSegments;
    type::vector< std::pair<Real, Real> > m_prevRigidCurvSegments;
    type::vector< bool >            m_rigidBeamList;
    type::vector<Transform>         m_vecGlobalHGravityCenter;
    std::map<Real, Transform> m_prevRigidTransforms;

    /// for re-interpolation
    type::vector<Transform>    m_vecGlobalHNode;
    type::vector<Deriv>        m_vecGlobalVelNode;

    /// for imposing nodes along the spline
    std::list< Real >    m_listOfImposedNodesOnXcurv;

    Data< type::vector<Real> > d_nodeCurvAbs;
    Data< type::vector<Vec2> > d_curvatureList;
    Data< VecCoord >     d_controlPoints;

    void dummyController(type::vector<Real> &newCurvAbs);

    /// If true update interpolation and subgraph on beginAnimationStep
    Data< bool >        d_updateOnBeginAnimationStep;
    Data< bool>         d_applyOrientationFirstInCreateNeedle;
    Data< bool >        d_reinitilizeWireOnInit;
    Data<type::vector<Real> > d_actualStepNoticeablePoints;

    /// Interface for topology changes
    TopologyContainer*  m_topology {nullptr} ;

    /// internal data necessary for computeSampling
    /// bool addingBeams is true when we are currently "adding" new Beams (for more precise computation)
    bool m_addingBeams;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This pattern is used to minize the amount of code inclusion into .h file and thus minimize the
/// overall compilation time of SOFA. In the .h these instances are declared as 'extern' meaning
/// they will not be instanciated. The actual instanciation is done in the corresponding .cpp file.
////////////////////////////////////////////////////////////////////////////////////////////////////
#if !defined(SOFA_PLUGIN_BEAMADAPTER_SUTURECONTROLLER_CPP)
extern template class SOFA_BEAMADAPTER_API SutureController<sofa::defaulttype::Rigid3Types>;
#endif
} /// namespace _suturecontroller_

using _suturecontroller_::SutureController ;

} /// namespace sofa::component::controller
