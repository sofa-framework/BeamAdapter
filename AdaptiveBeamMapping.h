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

#include "initBeamAdapter.h"
#include "BeamInterpolation.h"
#include "AdaptiveBeamController.h"
#include <sofa/core/Mapping.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <SofaBaseTopology/EdgeSetTopologyModifier.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/Event.h>
#include <sofa/core/visual/VisualParams.h>


using namespace sofa::component::fem;
using namespace sofa::core::objectmodel;

namespace sofa
{

namespace component
{
   /* namespace topology
    {
        template <class T>
        class EdgeSetGeometryAlgorithms;

        class EdgeSetTopologyModifier;
    }
    */

namespace mapping
{

/*!
 * \class AdaptiveBeamMapping
 * @brief AdaptiveBeamMapping Class
 */
template <class TIn, class TOut>
class AdaptiveBeamMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(AdaptiveBeamMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));

    typedef core::Mapping<TIn, TOut> Inherit;
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

    typedef helper::vector<unsigned int>                                         VecIndex    ;
    typedef sofa::core::topology::BaseMeshTopology::EdgeID                       ElementID   ;
    typedef sofa::helper::vector<sofa::core::topology::BaseMeshTopology::Edge>   VecEdges    ;
    typedef sofa::helper::vector<sofa::core::topology::BaseMeshTopology::EdgeID> VecElementID;

    typedef typename sofa::defaulttype::SolidTypes<InReal>::Transform      Transform       ;
    typedef std::pair<int, Transform>                                      IndexedTransform;
    typedef typename sofa::defaulttype::SolidTypes< InReal>::SpatialVector SpatialVector   ;

    typedef sofa::defaulttype::Vec<3, Real>   Vec3;
    typedef sofa::defaulttype::Vec<6, Real>   Vec6;
    typedef sofa::defaulttype::Mat<3,12,Real> Mat3x12;
    typedef sofa::defaulttype::Mat<12,3,Real> Mat12x3;
    typedef sofa::defaulttype::Mat<6,12,Real> Mat6x12;
    typedef sofa::defaulttype::Mat<12,6,Real> Mat12x6;
    typedef sofa::component::fem::BeamInterpolation< TIn > BInterpolation;

    typedef std::pair<unsigned int, Vec3> BeamIdAndBaryCoord;
    typedef struct
    {
       unsigned int beamId;
       //A bary point has 3 components
       // The first denote the curvilinear coordinate
       // The two followings denote the planar coordinate on the perpendicular cross section on the curve
       Vec3 baryPoint;
    } PosPointDefinition;


public:
    Data<bool> d_useCurvAbs;							/*!< true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same */
    Data< sofa::helper::vector< Vec3 > > d_points;	/*!< defines the mapped points along the beam axis (in beam frame local coordinates) */
    Data< double > d_proximity;						/*!< if positive, the mapping is modified for the constraints to take into account the lever created by the proximity */
    Data<bool> d_contactDuplicate;					/*!< if true, this mapping is a copy of an input mapping and is used to gather contact points (ContinuousFrictionContact Response) */
    Data<std::string> d_nameOfInputMap;				/*!< if contactDuplicate==true, it provides the name of the input mapping */
    Data<double> d_nbPointsPerBeam;					/*!< if non zero, we will adapt the points depending on the discretization, with this num of points per beam (compatible with useCurvAbs)*/
    Data< sofa::helper::vector< Real > > d_segmentsCurvAbs; /*!< (output) the abscissa of each created point on the collision model */
    SingleLink<AdaptiveBeamMapping<TIn, TOut>, BInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_adaptativebeamInterpolation;

    AdaptiveBeamMapping(core::State< In >* from=NULL, core::State< Out >* to=NULL,BeamInterpolation< TIn >* _interpolation=NULL, bool _isSubMapping=false)
        : Inherit(from, to)
    , d_useCurvAbs(initData(&d_useCurvAbs,true,"useCurvAbs","true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same"))
    , d_points(initData(&d_points, "points", "defines the mapped points along the beam axis (in beam frame local coordinates)"))
    , d_proximity(initData(&d_proximity, 0.0, "proximity", "if positive, the mapping is modified for the constraints to take into account the lever created by the proximity"))
    , d_contactDuplicate(initData(&d_contactDuplicate,false,"contactDuplicate","if true, this mapping is a copy of an input mapping and is used to gather contact points (ContinuousFrictionContact Response)"))
    , d_nameOfInputMap(initData(&d_nameOfInputMap,"nameOfInputMap", "if contactDuplicate==true, it provides the name of the input mapping"))
    , d_nbPointsPerBeam(initData(&d_nbPointsPerBeam, 0.0, "nbPointsPerBeam", "if non zero, we will adapt the points depending on the discretization, with this num of points per beam (compatible with useCurvAbs)"))
    , d_segmentsCurvAbs(initData(&d_segmentsCurvAbs, "segmentsCurvAbs", "the abscissa of each point on the collision model", true, true))
    , l_adaptativebeamInterpolation(initLink("interpolation", "Path to the Interpolation component on scene"), _interpolation)
    , m_inputMapping(NULL)
    , m_isSubMapping(_isSubMapping)
    , m_isBarycentricMapping(false)
    {
    }

    void printIstrumentInfo()const
    {
        if (m_isSubMapping)//ctn_DEV
        {
            std::cout<<"Instrument Named "<<l_adaptativebeamInterpolation->getName()<<std::endl
                    <<" MState1:"<<this->fromModel->getName()<< "  size:"<<this->fromModel->getSize()<<std::endl
                    <<" MState2:"<<this->toModel->getName()<< "  size:"<<this->toModel->getSize()<<std::endl
                    <<"idPointSubMap."<<m_idPointSubMap.size()<<std::endl
                    <<"pointBeamDistribution."<<m_pointBeamDistribution.size()<<std::endl<<std::endl;
        }
    }
    virtual ~AdaptiveBeamMapping(){}

    void apply(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<VecCoord>& out, const Data<InVecCoord>& in);

    void applyJ(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<VecDeriv>& out, const Data<InVecDeriv>& in);

    void applyJT(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<InVecDeriv>& out, const Data<VecDeriv>& in);

    void applyJT(const core::ConstraintParams *cparams /* PARAMS FIRST */, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in);


    int addPoint ( const Coord& c, int /*indexFrom*/ )
    {

        int i = d_points.getValue().size();
        Vec3 test(c[0],c[1],c[2]);

        d_points.beginEdit()->push_back(test);
        d_points.endEdit();
        return i;
    }

    int addContactPoint(const Vec3& bary);

    void setBarycentricMapping()
    {
        m_isBarycentricMapping=true;
        d_points.beginEdit()->clear();d_points.endEdit();
    }
    int addBaryPoint(const int& _beamId,const Vec3& _baryCoord,bool /*todo_straightline_spline_option*/)
    //TODO add parameter label for different cases : unactive, linear, spline
    {
        //attention, beamId here is not the edge Id, but the id of a vec_edge_list defined in BeamInterpolation
        int newpointId = m_pointBeamDistribution.size();
        m_pointBeamDistribution.resize(newpointId+1);
        m_pointBeamDistribution[newpointId].baryPoint=_baryCoord;
        m_pointBeamDistribution[newpointId].beamId=_beamId;
        return newpointId;
    }
    //void clear(){};////// CTN_DEV todo for ContactMapper
    //clear the mapping in functions of size given
    void clear(int size)
    {
        this->clearidPointSubMap();
        m_pointBeamDistribution.clear();
        if ( size>0 && !m_isSubMapping)
        {
            m_pointBeamDistribution.reserve ( size );
            d_points.beginEdit()->reserve ( size ); d_points.endEdit();
            this->getMechTo()[0]->resize(size);
        }
        else
        //case where this clear is call by a Multimapping, all component will be clear to null size
        {
            d_points.beginEdit()->resize(0); d_points.endEdit();
            this->getMechTo()[0]->resize(0);
        }
    }

    void computeIdxAndBaryCoordsForAbs(unsigned int &b, Real &x_bary, const Real &x_abs );

    void init();            // get the interpolation
    void bwdInit();        // get the points
    void reset(){init();  computeDistribution();}
    void reinit(){init(); computeDistribution();}

    void draw(const core::visual::VisualParams*);

    void beginAddContactPoint();

    void clearidPointSubMap(){m_idPointSubMap.clear();}
    void addidPointSubMap(unsigned int _id){m_idPointSubMap.push_back(_id);}
    void setuseCurvAbs(bool _value){d_useCurvAbs.setValue(_value);}

    const sofa::helper::vector< PosPointDefinition >& getPointBeamDistribution() const
    {
        return m_pointBeamDistribution;
    }

protected:

    void applyJonPoint(unsigned int i, SpatialVector& VNode0input, SpatialVector& VNode1input, Deriv& vOutput, const InVecCoord& x);

    void applyJTonPoint(unsigned int i, const Deriv& finput, SpatialVector& FNode0output, SpatialVector& FNode1output, const InVecCoord& x);

    void computeJacobianOnPoint(unsigned int i, const typename In::VecCoord& x);

    void computeDistribution();

    sofa::core::topology::TopologyContainer* m_topology;


    bool x_buf_used;
    typename In::VecCoord x_buf;

    sofa::helper::vector< PosPointDefinition > m_pointBeamDistribution;
    // for continuous_friction_contact:
    AdaptiveBeamMapping<TIn, TOut> *m_inputMapping;
    sofa::helper::vector<unsigned int> m_idPointSubMap;
    bool m_isSubMapping;
    bool m_isBarycentricMapping;

public :


};

} // namespace mapping

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_H */
