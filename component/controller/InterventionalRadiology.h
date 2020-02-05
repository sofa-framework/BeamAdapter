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
// C++ Implementation : InterventionalRadiology
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

#ifndef SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGY_H
#define SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGY_H

#include <SofaUserInteraction/MechanicalStateController.h>
#include <SofaBaseTopology/EdgeSetTopologyModifier.h>
#include <sofa/defaulttype/SolidTypes.h>
#include <SofaBoundaryCondition/FixedConstraint.h>
#include <sofa/core/DataEngine.h>
#include <SofaMeshCollision/PointModel.h>
#include <SofaMeshCollision/LineModel.h>
#include <sofa/helper/AdvancedTimer.h>
#include "../WireBeamInterpolation.h"


///FORWARD DECLARATION
namespace sofa {
namespace component {
namespace topology {
template <class T> class EdgeSetGeometryAlgorithms;
class EdgeSetTopologyModifier;
}
}
}


namespace sofa
{
namespace component
{
namespace controller
{
namespace _interventionalradiology_
{

using sofa::defaulttype::Vec;
using sofa::defaulttype::Vec3d;
using namespace sofa::component::fem;
using namespace sofa::helper;
using sofa::core::topology::BaseMeshTopology;
using sofa::helper::vector;
using sofa::component::projectiveconstraintset::FixedConstraint;


/*!
 * \class InterventionalRadiology
 * @brief InterventionalRadiology Class
 *
 * Provides a Mouse & Keyboard user control on an EdgeSet Topology.
 */
template<class DataTypes>
class InterventionalRadiology : public MechanicalStateController<DataTypes>,
        public collision::PointActiver,
        public collision::LineActiver
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(InterventionalRadiology,DataTypes),SOFA_TEMPLATE(MechanicalStateController,DataTypes));
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
    InterventionalRadiology();
    virtual ~InterventionalRadiology(){}

    ////////////////////// Inherited from BaseObject ///////////////////////////////////////////////
    virtual void init() override ;
    virtual void bwdInit() override ;
    virtual void reinit() override;
    virtual void draw(const core::visual::VisualParams*) override {}
    virtual std::string getTemplateName() const override;
    static std::string templateName(const InterventionalRadiology<DataTypes>* = NULL);
    ////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////// Inherited from Controller ///////////////////////////////////////////////
    virtual void onBeginAnimationStep(const double dt) override ;
    ////////////////////////////////////////////////////////////////////////////////////////////////


    ////////////////////// For Point & Line Activer interface //////////////////////////////////////
    virtual bool activePoint(int index, core::CollisionModel * cm = nullptr) override;
    virtual bool activeLine(int index, core::CollisionModel * cm = nullptr) override;
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

    /// Conditional elements for construction of InterventionalRadiology
    Data< vector< std::string > >  d_instrumentsPath;
    vector< WBeamInterpolation * > m_instrumentsList;

    /// For point and line activer
    vector<bool> m_activatedPointsBuf;

    /// Interface for interventionalRadiology instruments:
    virtual void applyInterventionalRadiology(void);

    void processDrop(unsigned int &previousNumControlledNodes,  unsigned int &seg_remove);
    void interventionalRadiologyComputeSampling(vector<Real> &newCurvAbs, vector< vector<int> > &id_instrument_table, const vector<Real> &xBegin, const Real& totalLengthCombined);
    /// Sort the curv Abs in the ascending order and avoid doubloon
    void sortCurvAbs(vector<Real> &CurvAbs,  vector< vector<int> >& id_instrument_table);
    void totalLengthIsChanging(const vector<Real>& newNodeCurvAbs, vector<Real>& modifiedNodeCurvAbs, const vector< vector<int> >& newTable);
    void fixFirstNodesWithUntil(unsigned int first_simulated_Node);
    void activateBeamListForCollision( vector<Real> &curv_abs, vector< vector<int> > &id_instrument_table);
    void loadMotionData(std::string filename);
    void getTotalLengthCombined(vector<Real> &newCurvAbs, vector<Real>& xbegin, Real &totalLengthCombined){
        Real xend;
        for (unsigned int i=0; i<m_instrumentsList.size(); i++)
        {

            xend= d_xTip.getValue()[i];
            msg_info() << "======> xend :"<< xend << "  m_instrumentsList[i]->getRestTotalLength() :"<< m_instrumentsList[i]->getRestTotalLength();
            Real xb = xend - m_instrumentsList[i]->getRestTotalLength();
            xbegin.push_back(xb);

            if (xend> totalLengthCombined)
                totalLengthCombined=xend;

            // clear the present interpolation of the beams
            m_instrumentsList[i]->clear();

            // create the first node (on x=0)
            if( xend > 0.0)
                newCurvAbs.push_back(0.0);
        }

    }

    //    void findTotalLengthCombined (vector<Real>& newCurvAbs, vector<Real> & xbegin, Real totalLengthCombined);

    Data<int>            d_controlledInstrument;
    Data<vector<Real>>   d_xTip;
    vector<Real>         d_old_xTip;
    Data<vector<Real>>   d_rotationInstrument;
    vector<Real>         d_old_rotationInstrument;
    Data<Real>           d_step;
    Data<Real>           d_angularStep;
    Data<Real>           d_speed;
    Data<Coord>          d_startingPos;
    Data<Real>           d_threshold;
    //    Data<vector<Real>>   d_rigidCurvAbs; // Pairs (start - end)
    Data<unsigned int>   d_indexFirstNode; // First Node simulated
    Data<vector<Real>>   d_curvAbs;
    Data<vector<Real>>   d_totalLengths;

    FixedConstraint<DataTypes> *    m_fixedConstraint;
    vector<Real>                    m_nodeCurvAbs;
    vector< vector<int> >           m_idInstrumentCurvAbsTable;
    unsigned int                    m_numControlledNodes; // Excluding the nodes that are "dropped"
};

} // namespace _interventionalradiologycontroller_

using _interventionalradiology_::InterventionalRadiology;

} // namespace controller

} // namespace component

} // namespace sofa

#endif /* SOFA_COMPONENT_CONTROLLER_INTERVENTIONALRADIOLOGY_H */
