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

#include "../WireBeamInterpolation.h"
#include <SofaUserInteraction/MechanicalStateController.h>
#include <SofaBaseTopology/EdgeSetTopologyModifier.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <SofaBoundaryCondition/FixedConstraint.h>
#include <sofa/core/DataEngine.h>
#include <SofaMeshCollision/PointModel.h>
#include <SofaMeshCollision/LineModel.h>


/////////////////////////// FORWARD DECLARATION ////////////////////////////////////////////////////
namespace sofa {
namespace component {
namespace topology {
    template <class T> class EdgeSetGeometryAlgorithms;
    class EdgeSetTopologyModifier;
    }
  }
}
////////////////////////////////////////////////////////////////////////////////////////////////////



/////////////////////////////////// DECLARATION ////////////////////////////////////////////////////
namespace sofa
{
namespace component
{
namespace controller
{
using sofa::defaulttype::Vec;
using sofa::defaulttype::Vec3d;
using namespace sofa::component::fem;
using namespace sofa::helper;


/*!
 * \class InterventionalRadiologyController
 * @brief InterventionalRadiologyController Class
 *
 * Provides a Mouse & Keyboard user control on an EdgeSet Topology.
 */
template<class DataTypes>
class InterventionalRadiologyController : public MechanicalStateController<DataTypes>,
                                          public collision::PointActiver,
                                          public collision::LineActiver
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(InterventionalRadiologyController,DataTypes),SOFA_TEMPLATE(MechanicalStateController,DataTypes));
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord    Coord   ;
    typedef typename DataTypes::Deriv    Deriv   ;
    typedef typename Coord::value_type   Real    ;

    typedef sofa::core::topology::BaseMeshTopology::EdgeID ElementID;
    typedef sofa::helper::vector<sofa::core::topology::BaseMeshTopology::EdgeID> VecElementID;

    typedef MechanicalStateController<DataTypes> Inherit;
    typedef sofa::component::fem::WireBeamInterpolation<DataTypes> WBeamInterpolation;

    typedef typename sofa::defaulttype::SolidTypes<Real>::Transform Transform;
    typedef typename sofa::defaulttype::SolidTypes<Real>::SpatialVector SpatialVector;

    typedef typename std::vector<Real>::const_iterator RealConstIterator;
    typedef Vec<3, Real> Vec3;

public:
    InterventionalRadiologyController();
    virtual ~InterventionalRadiologyController(){}

    ////////////////////// Inherited from BaseObject ///////////////////////////////////////////////
    virtual void init() override ;
    virtual void bwdInit() override ;
    virtual void reinit() override;
    virtual void draw(const core::visual::VisualParams*) override {}

    virtual std::string getTemplateName() const
    {
        return templateName(this);
    }

    static std::string templateName(const InterventionalRadiologyController<DataTypes>* = NULL)
    {
        return DataTypes::Name();
    }


    ////////////////////// Inherited from Controller ///////////////////////////////////////////////
    virtual void onMouseEvent(core::objectmodel::MouseEvent *) override ;
    virtual void onKeyPressedEvent(core::objectmodel::KeypressedEvent *) override ;
    virtual void onBeginAnimationStep(const double dt) override ;


    ////////////////////// For Point & Line Activer interface //////////////////////////////////////
    bool activePoint(int index, core::CollisionModel * /*cm*/ = 0)
    {
        return activated_Points_buf[index];
    }

    bool activeLine(int index, core::CollisionModel * /*cm*/ = 0)
    {
        return activated_Points_buf[index+1];
    }

    virtual bool modifyTopology(void){ return false;}

    void interventionalRadiologyCollisionControls(sofa::helper::vector<Real> &x_point_list,
                                                  sofa::helper::vector<int> &id_instrument_list,
                                                  sofa::helper::vector<int> &removeEdge);

    void getInstrumentList(sofa::helper::vector< sofa::component::fem::WireBeamInterpolation<DataTypes>  *>& list){
        list = m_instrumentsList;
    }

    const vector< vector<int> >& get_id_instrument_curvAbs_table()const {return id_instrument_curvAbs_table;}
    int getTotalNbEdges()const { return this->getContext()->getMeshTopology()->getNbEdges();}

protected:
    //conditional elements for construction of InterventionalRadiologyController
    Data< helper::vector< std::string > >  d_instrumentsPath;
    sofa::helper::vector< WBeamInterpolation * > m_instrumentsList;

    ///////////////////////////////// for point and line activer
    sofa::helper::vector<bool> activated_Points_buf;

    /////////////////////////// Interface for interventionalRadiology instruments:
    virtual void applyInterventionalRadiologyController(void);
    void processDrop(unsigned int &previousNumControlledNodes,  unsigned int &seg_remove);
    void interventionalRadiologyComputeSampling(sofa::helper::vector<Real> &newCurvAbs,
                                                sofa::helper::vector< sofa::helper::vector<int> > &id_instrument_table,
                                                const sofa::helper::vector<Real> &xBegin, const Real& xEnd);
    void sortCurvAbs(sofa::helper::vector<Real> &CurvAbs,  sofa::helper::vector< sofa::helper::vector<int> >& id_instrument_table);
    void totalLengthIsChanging(const sofa::helper::vector<Real>& newNodeCurvAbs, sofa::helper::vector<Real>& modifiedNodeCurvAbs,
                               const sofa::helper::vector< sofa::helper::vector<int> >& newTable);
    void fixFirstNodesWithUntil(unsigned int first_simulated_Node);

    void activateBeamListForCollision( sofa::helper::vector<Real> &curv_abs,
                                       sofa::helper::vector< sofa::helper::vector<int> > &id_instrument_table);

    void loadMotionData(std::string filename);

    Data<int> d_controlledInstrument;
    Data< sofa::helper::vector<Real> > d_xtip;
    Data< sofa::helper::vector<Real> > d_rotationInstrument;
    Data<Real> d_step;
    Data<Real> d_angularStep;
    Data<Real> d_speed;
    Data<Coord> d_startingPos;
    Data<Real> d_threshold;
    Data< helper::vector<Real> > d_rigidCurvAbs;	// Pairs (start - end)
    Data <std::string> d_motionFilename;
    Data<unsigned int> d_indexFirstNode; // First Node simulated
    Data< sofa::helper::vector<Real> > d_CurvAbs; //

    bool FF, RW, m_sensored;
    sofa::component::projectiveconstraintset::FixedConstraint<DataTypes> *m_fixedConstraint;
    sofa::helper::vector<int> droppedInstruments;

    sofa::helper::vector<Vec3d> sensorMotionData;
    unsigned int currentSensorData;

    sofa::helper::vector<Real> nodeCurvAbs;
    sofa::helper::vector< sofa::helper::vector<int> > id_instrument_curvAbs_table;
    unsigned int numControlledNodes; // excluding the nodes that are "dropped"

    bool dropCall;

    ///////////////////////// Interface for other Adaptive Control /////////////////////////////////
    sofa::core::topology::BaseMeshTopology* _topology;
    sofa::component::topology::EdgeSetGeometryAlgorithms<DataTypes>* edgeGeo;
    sofa::component::topology::EdgeSetTopologyModifier* edgeMod;
    Coord refPos;
    helper::vector<Real> vertexT; //=> replace by curvilinearAbs;

    virtual void computeVertexT();

    Real edgeTLength;

    ///////////////////////////// for rigidity control /////////////////////////////////////////////
    sofa::helper::vector< std::pair<Real, Real> > rigidCurveSegments, prevRigidCurvSegments;
    sofa::helper::vector< bool > rigidBeamList;
    sofa::helper::vector<Transform> vec_global_H_gravityCenter;
    std::map<Real, Transform> prevRigidTransforms;
};

} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGYCONTROLLER_H */
