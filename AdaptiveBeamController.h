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
// C++ Implementation : AdaptiveBeamController
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
//

#ifndef SOFA_COMPONENT_CONTROLLER_ADAPTIVEBEAMCONTROLLER_H
#define SOFA_COMPONENT_CONTROLLER_ADAPTIVEBEAMCONTROLLER_H

#include "BeamInterpolation.h"
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/component/controller/MechanicalStateController.h>
#include <sofa/component/topology/EdgeSetTopologyModifier.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/component/projectiveconstraintset/FixedConstraint.h>
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
 * @brief AdaptiveBeamController Class
 *
 * Provides a Mouse & Keyboard user control on an EdgeSet Topology.
 */
template<class DataTypes>
    class SOFA_BEAMADAPTER_API AdaptiveBeamController : public MechanicalStateController<DataTypes> ,
                                                        public sofa::component::collision::PointActiver,
                                                        public sofa::component::collision::LineActiver
{
public:
  SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamController,DataTypes),SOFA_TEMPLATE(MechanicalStateController,DataTypes));
	typedef typename DataTypes::VecCoord VecCoord;
	typedef typename DataTypes::VecDeriv VecDeriv;
	typedef typename DataTypes::Coord    Coord   ;
	typedef typename DataTypes::Deriv    Deriv   ;
	typedef typename Coord::value_type   Real    ;

    typedef sofa::core::topology::BaseMeshTopology::EdgeID ElementID;
    typedef sofa::helper::vector<sofa::core::topology::BaseMeshTopology::EdgeID> VecElementID;

	typedef MechanicalStateController<DataTypes> Inherit;

    typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
    typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;
      typedef  sofa::component::fem::BeamInterpolation<DataTypes> BInterpolation;

protected :
    BInterpolation* m_adaptiveinterpolation;
     Data< helper::vector< std::string > > m_interpolationPath;

public :

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

    typedef Vec<3, Real> Vec3;

	/**
	 * @brief Default Constructor.
	 */
    AdaptiveBeamController();

	/**
	 * @brief Default Destructor.
	 */
    virtual ~AdaptiveBeamController(){};

	/**
	 * @brief SceneGraph callback initialization method.
	 */
	virtual void init();

    virtual void reinit();

	/**
	 * @name Controller Interface
	 */
	//@{

	/**
	 * @brief Mouse event callback.
	 */
    virtual void onMouseEvent(core::objectmodel::MouseEvent *);

	/**
	 * @brief Keyboard key pressed event callback.
	 */
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *);


	/**
	 * @brief Begin Animation event callback.
	 */
	virtual void onBeginAnimationStep(const double dt);

	//@}

	/**
	 * @name Accessors
	 */
	//@{

	virtual std::string getTemplateName() const
    {
      return templateName(this);
    }

    static std::string templateName(const AdaptiveBeamController<DataTypes>* = NULL)
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
    virtual void draw(const core::visual::VisualParams*){}


protected:

    ////// for point and line activer
    sofa::helper::vector<Real> xAbs_collisionPoints_buf;

    Data<int> controlledInstrument;
    Data< sofa::helper::vector<Real> > xtip;
    Data< sofa::helper::vector<Real> > rotationInstrument;
    Data<Real> step;
    Data<Real> angularStep;
    Data<Real> speed;
    Data<Coord> startingPos;
    Data<Real> threshold;



    bool FF, RW;

    sofa::component::projectiveconstraintset::FixedConstraint<DataTypes> *_fixedConstraint;
    sofa::helper::vector<int> droppedInstruments;


    sofa::helper::vector<Real> nodeCurvAbs;
    sofa::helper::vector< sofa::helper::vector<int> > id_instrument_curvAbs_table;
    unsigned int numControlledNodes; // excluding the nodes that are "dropped"

    bool dropCall;


    /////////// Interface for other Adaptive Control

	sofa::core::topology::BaseMeshTopology* _topology;
	sofa::component::topology::EdgeSetGeometryAlgorithms<DataTypes>* edgeGeo;
	sofa::component::topology::EdgeSetTopologyModifier* edgeMod;
	Coord refPos;
    helper::vector<Real> vertexT; //=> replace by curvilinearAbs;
	Real edgeTLength;




};

} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_ADAPTIVEBEAMCONTROLLER_H */
