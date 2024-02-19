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
// C++ Implementation : InterventionalRadiologyController
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
#pragma once

#include <sofa/component/controller/MechanicalStateController.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyModifier.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <sofa/component/constraint/projective/FixedProjectiveConstraint.h>
#include <sofa/core/DataEngine.h>
#include <sofa/component/collision/geometry/PointModel.h>
#include <sofa/component/collision/geometry/LineModel.h>

#include <BeamAdapter/utils/BeamActions.h>
#include <BeamAdapter/component/WireBeamInterpolation.h>
#include <sofa/component/topology/container/dynamic/EdgeSetGeometryAlgorithms.h>
#include <sofa/component/topology/container/dynamic/EdgeSetTopologyModifier.h>

namespace sofa::component::controller
{

namespace _interventionalradiologycontroller_
{

using sofa::type::Vec;
using sofa::type::Vec3d;
using sofa::core::topology::BaseMeshTopology;
using sofa::component::constraint::projective::FixedProjectiveConstraint;

/*!
 * \class InterventionalRadiologyController
 * @brief InterventionalRadiologyController Class
 *
 * Provides a Mouse & Keyboard user control on an EdgeSet Topology.
 */
template<class DataTypes>
class InterventionalRadiologyController : public MechanicalStateController<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(InterventionalRadiologyController,DataTypes),SOFA_TEMPLATE(MechanicalStateController,DataTypes));
    typedef typename DataTypes::VecCoord                            VecCoord;
    typedef typename DataTypes::VecDeriv                            VecDeriv;
    typedef typename DataTypes::Coord                               Coord   ;
    typedef typename DataTypes::Deriv                               Deriv   ;
    typedef typename Coord::value_type                              Real    ;
    typedef typename defaulttype::SolidTypes<Real>::Transform       Transform;
    typedef typename defaulttype::SolidTypes<Real>::SpatialVector   SpatialVector;
    typedef typename std::vector<Real>::const_iterator              RealConstIterator;

    typedef Vec<3, Real>                            Vec3;
    typedef BaseMeshTopology::EdgeID                ElementID;
    typedef type::vector<BaseMeshTopology::EdgeID>        VecElementID;
    typedef MechanicalStateController<DataTypes>    Inherit;
    typedef fem::WireBeamInterpolation<DataTypes>   WBeamInterpolation;

public:
    InterventionalRadiologyController();
    virtual ~InterventionalRadiologyController() = default;

    ////////////////////// Inherited from BaseObject ///////////////////////////////////////////////
    virtual void init() override ;
    virtual void bwdInit() override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////// Inherited from Controller ///////////////////////////////////////////////
    virtual void onMouseEvent(core::objectmodel::MouseEvent *) override ;
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *) override ;
    virtual void onBeginAnimationStep(const double dt) override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////

    virtual bool modifyTopology(void);
    void interventionalRadiologyCollisionControls(type::vector<Real> &x_point_list,
                                                  type::vector<int> &id_instrument_list,
                                                  type::vector<int> &removeEdge);
    void getInstrumentList(type::vector<sofa::component::fem::WireBeamInterpolation<DataTypes>*>& list);
    const type::vector< type::vector<int> >& get_id_instrument_curvAbs_table()const;
    int getTotalNbEdges()const;

    void applyAction(sofa::beamadapter::BeamAdapterAction action);
    /// Method to warn this controller that a BeamActionController is controlling the scene. Will bypass the event handling in this component.
    void useBeamAction(bool value) { m_useBeamActions = value; }

    /// Getter to the tools curviline abscisses sorted @sa m_nodeCurvAbs at the current timestep.
    [[nodiscard]] const type::vector<Real>& getCurrentCurvAbscisses() const { return m_nodeCurvAbs; }


    /////////////////////////// Deprecated Methods  ////////////////////////////////////////// 
    [[deprecated("Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06")]]
    void processDrop(unsigned int& previousNumControlledNodes, unsigned int& seg_remove)
    {
        SOFA_UNUSED(previousNumControlledNodes);
        SOFA_UNUSED(seg_remove);
        msg_warning() << "Releasing catheter or brokenIn2 mode is not anymore supported. Feature has been removed after release v23.06";
    }


public:

    using Inherit1::f_printLog;
    using Inherit1::f_listening;
    using Inherit1::getContext;
    using Inherit1::getMechanicalState;

    /// Conditional elements for construction of InterventionalRadiologyController
    Data< type::vector< std::string > >  d_instrumentsPath;
    type::vector< WBeamInterpolation * > m_instrumentsList;

    /// For point and line activer
    type::vector<bool> m_activatedPointsBuf;

    /// Interface for interventionalRadiology instruments:
    virtual void applyInterventionalRadiologyController(void);
    
    
private:
    /** Compute the sambling curv abscisses using each instrument sampling and key points parameters
    * Will call @sa sortCurvAbs to sort the curv abs and remove doubloon
    * Need each tool starting position to sample only activated nodes and tool total length (combined deployed tool lengths)
    **/
    void computeInstrumentsCurvAbs(type::vector<Real>& newCurvAbs, const type::vector<Real>& tools_xBegin, const Real& totalLength);

    /// Method to sort the curv Abs in the ascending order and remove doubloon that are closer than d_threshold
    void sortCurvAbs(type::vector<Real>& curvAbs);

    /// Method to fill the id_instrument_table based on curvAbs and each tool begin and end. The table as the same size as the curvAbs buffer and store for each curvAbs[i] a vector with the ids of the instruments that are present at this position.
    void fillInstrumentCurvAbsTable(const type::vector<Real>& curvAbs, const type::vector<Real>& tools_xBegin, const type::vector<Real>& tools_xEnd, type::vector< type::vector<int> >& id_instrument_table);

public:
    void totalLengthIsChanging(const type::vector<Real>& newNodeCurvAbs, type::vector<Real>& modifiedNodeCurvAbs, const type::vector< type::vector<int> >& newTable);
    void fixFirstNodesWithUntil(unsigned int first_simulated_Node);
    void activateBeamListForCollision( type::vector<Real> &curv_abs, type::vector< type::vector<int> > &id_instrument_table);
    void loadMotionData(std::string filename);

    Data<int>            d_controlledInstrument;
    Data<type::vector<Real>>   d_xTip;
    Data<type::vector<Real>>   d_rotationInstrument;
    Data<Real>           d_step;
    Data<Real>           d_angularStep;
    Data<Real>           d_speed;
    Data<Coord>          d_startingPos;
    Data<Real>           d_threshold;
    Data<type::vector<Real>>   d_rigidCurvAbs; // Pairs (start - end)
    Data<std::string>    d_motionFilename;
    Data<unsigned int>   d_indexFirstNode; // First Node simulated
    
    
    bool m_useBeamActions = false;
    bool m_FF, m_RW, m_sensored;
    FixedProjectiveConstraint<DataTypes> *    m_fixedConstraint;
    type::vector<Vec3d>                   m_sensorMotionData;
    unsigned int                    m_currentSensorData;
    type::vector<Real>                    m_nodeCurvAbs;
    type::vector< type::vector<int> >           m_idInstrumentCurvAbsTable;
    unsigned int                    m_numControlledNodes; // Excluding the nodes that are "dropped"
};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_INTERVENTIONALRADIOCONTROLLER_CPP)
extern template class SOFA_BEAMADAPTER_API InterventionalRadiologyController<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace _interventionalradiologycontroller_

using _interventionalradiologycontroller_::InterventionalRadiologyController;

} // namespace sofa::component::controller
