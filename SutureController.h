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
// This controller allows to drive the disretization of a suture model computed with adaptive beam elements
// + The disretization is driven by: the curvature, some particular points on the wire (the node on the boundary between needle and thread), the suture constraints
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

#include "WireRestShape.h" // needed ??
#include "WireBeamInterpolation.h"
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/controller/MechanicalStateController.h>
#include <sofa/component/topology/EdgeSetTopologyModifier.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/core/DataEngine.h>
#include <sofa/component/collision/PointModel.h>
#include <sofa/component/collision/LineModel.h>


using namespace sofa::component::fem;

namespace sofa
{

namespace component
{
	namespace topology
	{
		template <class T>
		class EdgeSetGeometryAlgorithms;

		class EdgeSetTopologyModifier;
	}


namespace controller
{


/**
 * @brief SutureController Class
 *
 * Provides a Mouse & Keyboard user control on an EdgeSet Topology.
 */
template<class DataTypes>
    class SOFA_BEAMADAPTER_API SutureController : public MechanicalStateController<DataTypes> ,
                                                        public sofa::component::collision::PointActiver,
                                                        public sofa::component::collision::LineActiver
{
public:
  SOFA_CLASS(SOFA_TEMPLATE(SutureController,DataTypes),SOFA_TEMPLATE(MechanicalStateController,DataTypes));
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::Coord    Coord   ;
	typedef typename DataTypes::Deriv    Deriv   ;
	typedef typename Coord::value_type   Real    ;
     typedef Vec<3, Real> Vec3;

    typedef sofa::core::topology::BaseMeshTopology::EdgeID ElementID;
    typedef sofa::helper::vector<sofa::core::topology::BaseMeshTopology::EdgeID> VecElementID;

    typedef typename std::list< Real >::iterator ListRealIterator;

    typedef MechanicalStateController<DataTypes> Inherit;

    typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
    typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;
    typedef sofa::component::fem::WireBeamInterpolation<DataTypes> WInterpolation;


protected :
    WInterpolation* m_adaptiveinterpolation;
    //sofa::core::objectmodel::DataObjectRef m_interpolationPath;

public :
//    void setPathToInterpolation(const std::string &o){
//        m_interpolationPath.getValue()[0] = o;
//    }

    /////////////// Point & Line Activer interface
    bool activePoint(int index, core::CollisionModel * /*cm*/ = 0)
    {

        if (index >= (int)xAbs_collisionPoints_buf.size() || index<0)
        return false;

        if(xAbs_collisionPoints_buf[index]>10.0)
            return true;

        return false;

    }

    bool activeLine(int index, core::CollisionModel * /*cm*/ = 0)
    {

        if ((index+1) >= (int)xAbs_collisionPoints_buf.size() || (index+1)<0)
        return false;

        if(xAbs_collisionPoints_buf[index+1]>10.0)
            return true;

        return false;
    }


    /////////////////////



	/**
	 * @brief Default Constructor.
	 */
    SutureController(WireBeamInterpolation<DataTypes>* _adaptiveinterpolation);
    SutureController();

	/**
	 * @brief Default Destructor.
	 */
    virtual ~SutureController(){};

	/**
	 * @brief SceneGraph callback initialization method.
	 */
	virtual void init();

    virtual void bwdInit(){applyController();}

    virtual void reinit();

    virtual void reset(){ init(); applyController();}

	/**
	 * @name Controller Interface
	 */
	//@{

	/**
	 * @brief Mouse event callback.
	 */
    virtual void onMouseEvent(core::objectmodel::MouseEvent *){}

	/**
	 * @brief Keyboard key pressed event callback.
	 */
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *){}


	/**
	 * @brief Begin Animation event callback.
	 */
        virtual void onEndAnimationStep(const double dt);

	//@}

	/**
	 * @name Accessors
	 */
	//@{

	virtual std::string getTemplateName() const
    {
      return templateName(this);
    }

    static std::string templateName(const SutureController<DataTypes>* = NULL)
    {
      return DataTypes::Name();
    }

	//@}


	/**
	 * @brief Apply the controller modifications to the controlled MechanicalState.
	 */
	virtual void applyController(void);

	/**
	 * @brief
	 */
        virtual bool modifyTopology(void){ return false;}
 
	/**
	 * @brief
	 */
        virtual void draw(const core::visual::VisualParams*);




    /**
      * @brief addNodeOnXcurv (const Real& x_curv) will add a node at abs curv "x_curv" in the list of the imposed node
      *        this list is the used in the controller for the sampling of nodes of the suture model
      */
    void addNodeOnXcurv(const Real& x_curv){   listOfImposedNodesOnXcurv.push_back(x_curv);}
    void clearNodesOnXcurv(){listOfImposedNodesOnXcurv.clear();}



private:
    void recreateTopology();
    void addNodesAndEdge(unsigned int num, Real &xend);
    void removeNodesAndEdge(unsigned int num);

    // internal function: compute the bending angle between curv abs xmin and xmax
    void computeBendingAngle(Real& angle, const Real& xmin, const Real& xmax, const Real& dx_comput, const VecCoord& Pos);


    // this function computes the tangent value on a series of discrete points (store also the curv_abs of these discrete points)
    void computeTangentOnDiscretePoints(sofa::helper::vector<Vec3> TangTable, sofa::helper::vector<Real> xTable,  unsigned int numDiscretePoints, const VecCoord& Pos);

    // fill the list rigidBeamList
    void detectRigidBeams(sofa::helper::vector<Real> &newCurvAbs);

    // when the sampling is computed, this function allows for reinterpolate the position and the velocity
    // TODO !! ADD 1/ RIGIDIFICATIONS + 2/ TOPOLOGY CHANGES + 3/ADD_BEAM (on the adaptive interpolation)
    void applyNewSampling(sofa::helper::vector<Real> &newCurvAbs, sofa::helper::vector<Real> &oldCurvAbs, VecCoord &x, VecDeriv &v);


    // add the curv abs of the nodes at the extremity of the rigid segment
    // if a node already exists or is very close (< tol), do not add any point
    void addRigidCurvAbs(sofa::helper::vector<Real> &newCurvAbs, const Real &tol);


    // add the nodes that are imposed at a given curv abs
    // if a node already exists or is very close (< tol), do not add any point
    void addImposedCurvAbs(sofa::helper::vector<Real> &newCurvAbs, const Real &tol);


    // compute sampling
    void computeSampling(sofa::helper::vector<Real> &newCurvAbs, VecCoord& x);
    // internal data necessary for computeSampling
    // bool addingBeams is true when we are currently "adding" new Beams (for more precise computation)
    bool addingBeams;


    // little function to verify that the rigidCurveSegments are sorted
    bool verifyRigidCurveSegmentSort();




protected:

    ////// for point and line activer
    sofa::helper::vector<Real> xAbs_collisionPoints_buf;

    ///// Data:
    Data<Coord> startingPos;
    Data<Real> threshold;
    Data< Real > maxBendingAngle;
    Data< helper::vector< std::string > > m_interpolationPath;
    Data<bool> useDummyController;

    /////// for rigidity control
    sofa::helper::vector< std::pair<Real, Real> > rigidCurveSegments;
    sofa::helper::vector< bool > rigidBeamList;
    sofa::helper::vector<Transform> vec_global_H_gravityCenter;

    /////// for re-interpolation
    sofa::helper::vector<Transform> vec_global_H_node;
    sofa::helper::vector<Deriv> vec_global_Vel_node;

    /////// for imposing nodes along the spline
    std::list< Real > listOfImposedNodesOnXcurv;





    sofa::helper::vector<Real> cutCurvAbs; // store the curv abs where the thread is cut
    sofa::helper::vector<Real> nodeCurvAbs;





    /////////// Interface for topology changes

    sofa::core::topology::TopologyContainer* _topology;
    sofa::component::topology::EdgeSetGeometryAlgorithms<DataTypes>* edgeGeo;
    sofa::component::topology::EdgeSetTopologyModifier* edgeMod;

    void dummyController(sofa::helper::vector<Real> &newCurvAbs);


};

} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_SUTURECONTROLLER_H */
