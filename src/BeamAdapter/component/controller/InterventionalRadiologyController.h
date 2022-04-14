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

#ifndef SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGYCONTROLLER_H
#define SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGYCONTROLLER_H

#include <SofaUserInteraction/MechanicalStateController.h>
#include <SofaBaseTopology/EdgeSetTopologyModifier.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <SofaBoundaryCondition/FixedConstraint.h>
#include <sofa/core/DataEngine.h>
#include <SofaMeshCollision/PointModel.h>
#include <SofaMeshCollision/LineModel.h>

#include <BeamAdapter/component/WireBeamInterpolation.h>


///FORWARD DECLARATION
namespace sofa::component::topology::container::dynamic
{
    template <class T> 
    class EdgeSetGeometryAlgorithms;
    class EdgeSetTopologyModifier;
}


namespace sofa
{
namespace component
{
namespace controller
{
namespace _interventionalradiologycontroller_
{

using sofa::type::Vec;
using sofa::type::Vec3d;
using namespace sofa::component::fem;
using namespace sofa::helper;
using sofa::core::topology::BaseMeshTopology;
using sofa::type::vector;
using sofa::component::projectiveconstraintset::FixedConstraint;


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
    typedef vector<BaseMeshTopology::EdgeID>        VecElementID;
    typedef MechanicalStateController<DataTypes>    Inherit;
    typedef fem::WireBeamInterpolation<DataTypes>   WBeamInterpolation;

public:
    InterventionalRadiologyController();
    virtual ~InterventionalRadiologyController(){}

    ////////////////////// Inherited from BaseObject ///////////////////////////////////////////////
    virtual void init() override ;
    virtual void bwdInit() override ;
    virtual void reinit() override;
    virtual void draw(const core::visual::VisualParams*) override {}
    virtual std::string getTemplateName() const override;
    static std::string templateName(const InterventionalRadiologyController<DataTypes>* = NULL);
    ////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////// Inherited from Controller ///////////////////////////////////////////////
    virtual void onMouseEvent(core::objectmodel::MouseEvent *) override ;
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *) override ;
    virtual void onBeginAnimationStep(const double dt) override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////////////////////////////

    virtual bool modifyTopology(void);
    void interventionalRadiologyCollisionControls(vector<Real> &x_point_list,
                                                  vector<int> &id_instrument_list,
                                                  vector<int> &removeEdge);
    void getInstrumentList(vector<sofa::component::fem::WireBeamInterpolation<DataTypes>*>& list);
    const vector< vector<int> >& get_id_instrument_curvAbs_table()const;
    int getTotalNbEdges()const;

public:

    using Inherit1::f_printLog;
    using Inherit1::f_listening;
    using Inherit1::getContext;
    using Inherit1::getMechanicalState;

    /// Conditional elements for construction of InterventionalRadiologyController
    Data< vector< std::string > >  d_instrumentsPath;
    vector< WBeamInterpolation * > m_instrumentsList;

    /// For point and line activer
    vector<bool> m_activatedPointsBuf;

    /// Interface for interventionalRadiology instruments:
    virtual void applyInterventionalRadiologyController(void);
    virtual void computeVertexT();

    void processDrop(unsigned int &previousNumControlledNodes,  unsigned int &seg_remove);
    void interventionalRadiologyComputeSampling(vector<Real> &newCurvAbs, vector< vector<int> > &id_instrument_table, const vector<Real> &xBegin, const Real& xEnd);
    /// Sort the curv Abs in the ascending order and avoid doubloon
    void sortCurvAbs(vector<Real> &CurvAbs,  vector< vector<int> >& id_instrument_table);
    void totalLengthIsChanging(const vector<Real>& newNodeCurvAbs, vector<Real>& modifiedNodeCurvAbs, const vector< vector<int> >& newTable);
    void fixFirstNodesWithUntil(unsigned int first_simulated_Node);
    void activateBeamListForCollision( vector<Real> &curv_abs, vector< vector<int> > &id_instrument_table);
    void loadMotionData(std::string filename);

    Data<int>            d_controlledInstrument;
    Data<vector<Real>>   d_xTip;
    Data<vector<Real>>   d_rotationInstrument;
    Data<Real>           d_step;
    Data<Real>           d_angularStep;
    Data<Real>           d_speed;
    Data<Coord>          d_startingPos;
    Data<Real>           d_threshold;
    Data<vector<Real>>   d_rigidCurvAbs; // Pairs (start - end)
    Data<std::string>    d_motionFilename;
    Data<unsigned int>   d_indexFirstNode; // First Node simulated
    Data<vector<Real>>   d_curvAbs;

    bool m_FF, m_RW, m_sensored;
    FixedConstraint<DataTypes> *    m_fixedConstraint;
    vector<int>                     m_droppedInstruments;
    vector<Vec3d>                   m_sensorMotionData;
    unsigned int                    m_currentSensorData;
    vector<Real>                    m_nodeCurvAbs;
    vector< vector<int> >           m_idInstrumentCurvAbsTable;
    unsigned int                    m_numControlledNodes; // Excluding the nodes that are "dropped"
    bool                            m_dropCall;
};

} // namespace _interventionalradiologycontroller_

using _interventionalradiologycontroller_::InterventionalRadiologyController;

} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGYCONTROLLER_H */
