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

#ifndef SOFA_COMPONENT_CONTROLLER_SUTURECONTROLLER_H
#define SOFA_COMPONENT_CONTROLLER_SUTURECONTROLLER_H
#include <sofa/helper/set.h>

#include <SofaUserInteraction/MechanicalStateController.h>
#include <SofaMeshCollision/PointModel.h>
#include <SofaMeshCollision/LineModel.h>

#include "../WireBeamInterpolation.h"
#include "../initBeamAdapter.h"

////////////////////////////////////////////////////////////////////////////////////////////////////
/// Forward declarations, see https://en.wikipedia.org/wiki/Forward_declaration
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa {
namespace component {
namespace topology {
    template <class T> class EdgeSetGeometryAlgorithms;
    class EdgeSetTopologyModifier;
    }
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////DECLARATIONS /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa
{
namespace component
{
namespace controller
{
namespace _suturecontroller_
{

using std::set ;
using sofa::helper::vector ;
using sofa::defaulttype::Vec ;
using sofa::defaulttype::Rigid3dTypes ;
using sofa::defaulttype::Rigid3fTypes ;
using sofa::defaulttype::SolidTypes ;
using sofa::core::topology::TopologyContainer ;
using sofa::core::CollisionModel ;
using sofa::component::collision::LineActiver;
using sofa::component::collision::PointActiver;
using sofa::core::topology::BaseMeshTopology ;
using sofa::component::fem::WireBeamInterpolation;

/**
 * \class SutureController
 * @brief SutureController Class
 *
 * Provides a Mouse & Keyboard user control on an EdgeSet Topology.
 */
template<class DataTypes>
class SutureController : public MechanicalStateController<DataTypes>,
                         public PointActiver,
                         public LineActiver
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
    typedef vector< ElementID > VecElementID;

    typedef typename std::list< Real >::iterator ListRealIterator;

    typedef MechanicalStateController<DataTypes> Inherit;

    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;
    typedef WireBeamInterpolation<DataTypes> WInterpolation;

    typedef typename set<Real>::const_iterator RealConstIterator;

public:
    SutureController(WInterpolation* _adaptiveinterpolation = NULL) ;
    virtual ~SutureController(){}

    /////////////////// Inherited from Base //////////////////////////////////////////////////
    virtual void init();
    virtual void bwdInit(){applyController();}
    virtual void reinit();
    virtual void reset(){ init(); applyController();}
    virtual void draw(const core::visual::VisualParams*);

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const SutureController<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }

    /////////////////// Inherited from Controller //////////////////////////////////////////////////
    virtual void onMouseEvent(core::objectmodel::MouseEvent *)override {}
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *) override {}
    virtual void onBeginAnimationStep(const double dt) override ;
    virtual void onEndAnimationStep(const double dt) override ;

    /////////////////// Inherited from PointActiver //////////////////////////////////////////////////
    virtual bool activePoint(int index, CollisionModel * /*cm*/ = 0) override ;

    /////////////////// Inherited from LineActiver //////////////////////////////////////////////////
    virtual bool activeLine(int index, CollisionModel * /*cm*/ = 0) override ;

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
    Real computeBendingAngle(const Real& xmin, const Real& xmax, const Real& dx_comput, const VecCoord& Pos);


    /// Computes the tangent value on a series of discrete points (store also the curv_abs of these discrete points)
    void computeTangentOnDiscretePoints(vector<Vec3> TangTable, vector<Real> xTable,  unsigned int numDiscretePoints, const VecCoord& Pos);


    /// fill the list rigidBeamList
    void detectRigidBeams(const vector<Real> &newCurvAbs);


    /// when the sampling is computed, this function allows for reinterpolate the position and the velocity
    // TODO !! ADD 1/ RIGIDIFICATIONS + 2/ TOPOLOGY CHANGES + 3/ADD_BEAM (on the adaptive interpolation)
    void applyNewSampling(const vector<Real> &newCurvAbs, const vector<Real> &oldCurvAbs, VecCoord &x, VecDeriv &v);


    /// Add the curv abs of the nodes at the extremity of the rigid segment
    /// if a node already exists or is very close (< tol), do not add any point
    void addRigidCurvAbs(vector<Real> &newCurvAbs, const Real &tol);


    /// Add the nodes that are imposed at a given curv abs
    /// if a node already exists or is very close (< tol), do not add any point
    void addImposedCurvAbs(vector<Real> &newCurvAbs, const Real &tol);


    /// Computes sampling
    void computeSampling(vector<Real> &newCurvAbs, VecCoord& x);

    /// little function to verify that the rigidCurveSegments are sorted
    bool verifyRigidCurveSegmentSort();

    /// make sure the sampling does not change for the rigid segments
    void verifyRigidSegmentsSampling(vector<Real> &newCurvAbs);

    void storeRigidSegmentsTransformations();
    void verifyRigidSegmentsTransformations();

    void updateControlPointsPositions();

private:
    ////// for point and line activer
    vector<Real> m_XAbsCollisionPointsBuffer;

    ///// Data:
    Data< Coord >     d_startingPos;
    Data< Real >      d_threshold;
    Data< Real >      d_maxBendingAngle;
    Data< bool >      d_useDummyController;
    Data< bool >      d_fixRigidTransforms;

    /// Values in the set are interpreted as pairs (start - end) curve absciss.
    Data< set<Real> > d_rigidCurvAbs;
    SingleLink<SutureController<DataTypes>, WInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_adaptiveinterpolation;

    /////// for rigidity control
    vector< std::pair<Real, Real> > m_rigidCurveSegments;
    vector< std::pair<Real, Real> > m_prevRigidCurvSegments;
    vector< bool >            m_rigidBeamList;
    vector<Transform>         m_VecGlobalHGravityCenter;
    std::map<Real, Transform> m_prevRigidTransforms;

    /////// for re-interpolation
    vector<Transform>    m_VecGlobalHNode;
    vector<Deriv>        m_vecGlobalVelNode;

    /////// for imposing nodes along the spline
    std::list< Real >    m_listOfImposedNodesOnXcurv;

    Data< vector<Real> > d_nodeCurvAbs;
    Data< vector<Vec2> > d_curvatureList;
    Data< VecCoord >     d_controlPoints;

    /////////// Interface for topology changes
    TopologyContainer*  m_topology {nullptr} ;

    void dummyController(vector<Real> &newCurvAbs);

    //// If true update interpolation and subgraph on beginAnimationStep
    Data< bool >        d_updateOnBeginAnimationStep;
    Data< bool>         d_applyOrientationFirstInCreateNeedle;
    Data< bool >        d_reinitilizeWireOnInit;
    Data<vector<Real> > d_actualStepNoticeablePoints;

    /// internal data necessary for computeSampling
    /// bool addingBeams is true when we are currently "adding" new Beams (for more precise computation)
    bool addingBeams;
};

////////////////////////////////////////////////////////////////////////////////////////////////////
/// This pattern is used to minize the amount of code inclusion into .h file and thus minimize the
/// overall compilation time of SOFA. In the .h these instances are declared as 'extern' meaning
/// they will not be instanciated. The actual instanciation is done in the corresponding .cpp file.
////////////////////////////////////////////////////////////////////////////////////////////////////
#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_SUTURECONTROLLER_CPP)
#ifdef SOFA_WITH_DOUBLE
extern template class SOFA_BEAMADAPTER_API SutureController<Rigid3dTypes>;
#endif
#ifdef SOFA_WITH_FLOAT
extern template class SOFA_BEAMADAPTER_API SutureController<Rigid3fTypes>;
#endif
#endif
} /// namespace _suturecontroller_

using _suturecontroller_::SutureController ;

} /// namespace controller

} /// namespace component

} /// namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_SUTURECONTROLLER_H */
