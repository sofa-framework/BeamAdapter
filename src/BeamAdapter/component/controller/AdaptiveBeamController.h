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

//////////////////////// Inclusion of headers...from wider to narrower/closer //////////////////////
#include <SofaUserInteraction/MechanicalStateController.h>
#include <SofaBoundaryCondition/FixedConstraint.h>
#include <SofaMeshCollision/PointModel.h>
#include <SofaMeshCollision/LineModel.h>

#include <BeamAdapter/component/BeamInterpolation.h>
#include <BeamAdapter/component/BeamInterpolation.h>


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Forward declarations, see https://en.wikipedia.org/wiki/Forward_declaration
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa {
    namespace component {
        namespace topology {
            template <class T>
            class EdgeSetGeometryAlgorithms;
            class EdgeSetTopologyModifier;
        }
    }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// Declarations
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace sofa
{
namespace component
{
namespace controller
{

/////////////////////////////////// private namespace pattern //////////////////////////////////////
/// To avoid the lacking of names imported with with 'using' in the other's component namespace
/// you should use a private namespace and "export" only this one in the public namespace.
/// This is done at the end of this file, have a look if you are not used to this pattern.
////////////////////////////////////////////////////////////////////////////////////////////////////
namespace _adaptivebeamcontroller_
{

using sofa::component::projectiveconstraintset::FixedConstraint ;
using sofa::component::topology::EdgeSetTopologyModifier ;
using sofa::component::topology::EdgeSetGeometryAlgorithms ;
using sofa::component::fem::BeamInterpolation ;
using sofa::core::objectmodel::KeypressedEvent ;
using sofa::core::objectmodel::MouseEvent ;
using sofa::core::topology::BaseMeshTopology ;
using sofa::core::CollisionModel ;
using sofa::defaulttype::SolidTypes ;
using sofa::type::Vec ;
using sofa::type::vector ;
using std::string;

// TODO(dmarchal 2017-05-17) to eulalie & christian is the following still valid ?
/**
 * \class AdaptiveBeamController
 * @brief AdaptiveBeamController Mouse & Keyboard controller for EdgeSetTopology
 *
 * This component provides an interaction technique based on Mouse & Keyboard that allow user to
 * control an EdgeSet Topology.
 *
 * More informations about SOFA components:
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/create-your-component/
 * https://www.sofa-framework.org/community/doc/programming-with-sofa/components-api/components-and-datas/
 */
template<class DataTypes>
class AdaptiveBeamController : public MechanicalStateController<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(AdaptiveBeamController,DataTypes),
               SOFA_TEMPLATE(MechanicalStateController,DataTypes));

    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    typedef BaseMeshTopology::EdgeID ElementID;
    typedef vector<BaseMeshTopology::EdgeID> VecElementID;

    typedef MechanicalStateController<DataTypes> Inherit;

    typedef typename SolidTypes<Real>::Transform Transform;
    typedef typename SolidTypes<Real>::SpatialVector SpatialVector;
    typedef BeamInterpolation<DataTypes> BInterpolation;

    typedef Vec<3, Real> Vec3;

public :
    AdaptiveBeamController();
    virtual ~AdaptiveBeamController(){}

    /////////////// Inherited from BaseObject  /////////////////////////////////////////////////////
    virtual void init() override ;
    virtual void reinit() override ;

    /////////////// Inherited from MechanicalStateController  //////////////////////////////////////
    virtual void onMouseEvent(MouseEvent *) override ;
    virtual void onKeyPressedEvent(KeypressedEvent *) override ;
    virtual void onBeginAnimationStep(const double dt) override ;

    //TODO(dmarchal 2017-05-17) Check that these two are really needed (remove 1 one year if not done)
    virtual string getTemplateName() const
    {
      return templateName(this);
    }

    static string templateName(const AdaptiveBeamController<DataTypes>* = NULL)
    {
      return DataTypes::Name();
    }

protected:
    void applyController(void) ;

    Data<vector<string>>   d_interpolationPath;
    Data<int>                   d_controlledInstrument;
    Data<vector<Real>>          d_xtip;
    Data<vector<Real>>          d_rotationInstrument;
    Data<Real>                  d_step;
    Data<Real>                  d_angularStep;
    Data<Real>                  d_speed;

    BInterpolation*             m_adaptiveinterpolation {nullptr};

    ////// for point and line activer
    vector<Real>                m_xAbsCollisionPointsBuffer;

    bool                        FF {false} ;
    bool                        RW {false} ;
};

} /// namespace _adaptivebeamcontroller_

/// 'Export' the objects defined in the private namespace into the 'public' one.
using _adaptivebeamcontroller_::AdaptiveBeamController ;

} /// namespace controller

} /// namespace component

} /// namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_ADAPTIVEBEAMCONTROLLER_H */
