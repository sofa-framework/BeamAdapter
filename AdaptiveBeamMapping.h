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
#include <sofa/component/topology/EdgeSetTopologyModifier.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/objectmodel/Event.h>


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
    Data<bool> useCurvAbs;							/*!< true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same */
    Data< sofa::helper::vector< Vec3 > > points;	/*!< defines the mapped points along the beam axis (in beam frame local coordinates) */
    Data< double > proximity;						/*!< if positive, the mapping is modified for the constraints to take into account the lever created by the proximity */
    Data<bool> contactDuplicate;					/*!< if true, this mapping is a copy of an input mapping and is used to gather contact points (ContinuousFrictionContact Response) */
    Data<std::string> nameOfInputMap;				/*!< if contactDuplicate==true, it provides the name of the input mapping */
    SingleLink<AdaptiveBeamMapping<TIn, TOut>, BInterpolation, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> m_adaptativebeamInterpolation;
	
    AdaptiveBeamMapping()
        : Inherit()
    , useCurvAbs(initData(&useCurvAbs,true,"useCurvAbs","true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same"))
    , points(initData(&points, "points", "defines the mapped points along the beam axis (in beam frame local coordinates)"))
    , proximity(initData(&proximity, 0.0, "proximity", "if positive, the mapping is modified for the constraints to take into account the lever created by the proximity"))
    , contactDuplicate(initData(&contactDuplicate,false,"contactDuplicate","if true, this mapping is a copy of an input mapping and is used to gather contact points (ContinuousFrictionContact Response)"))
    , nameOfInputMap(initData(&nameOfInputMap,"nameOfInputMap", "if contactDuplicate==true, it provides the name of the input mapping"))
    , m_adaptativebeamInterpolation(initLink("interpolation", "Path to the Interpolation component on scene"))
    , m_inputMapping(NULL)
    , isSubMapping(false)
    , isBarycentricMapping(false)
    {
	}

    AdaptiveBeamMapping(core::State< In >* from, core::State< Out >* to,BeamInterpolation< TIn >* _interpolation=NULL,bool _isSubMapping=false)
        : Inherit(from, to)
    , useCurvAbs(initData(&useCurvAbs,true,"useCurvAbs","true if the curvilinear abscissa of the points remains the same during the simulation if not the curvilinear abscissa moves with adaptivity and the num of segment per beam is always the same"))
    , points(initData(&points, "points", "defines the mapped points along the beam axis (in beam frame local coordinates)"))
    , proximity(initData(&proximity, 0.0, "proximity", "if positive, the mapping is modified for the constraints to take into account the lever created by the proximity"))
    , contactDuplicate(initData(&contactDuplicate,false,"contactDuplicate","if true, this mapping is a copy of an input mapping and is used to gather contact points (ContinuousFrictionContact Response)"))
    , nameOfInputMap(initData(&nameOfInputMap,"nameOfInputMap", "if contactDuplicate==true, it provides the name of the input mapping"))
    , m_adaptativebeamInterpolation(initLink("interpolation", "Path to the Interpolation component on scene"))
    , m_inputMapping(NULL)
    , isSubMapping(_isSubMapping)
    , isBarycentricMapping(false)
    {
		if(_interpolation)
			m_adaptativebeamInterpolation.set(_interpolation);
    }

	void printIstrumentInfo()const
	{
		if (isSubMapping)//ctn_DEV
		{
			std::cout<<"Instrument Named "<<m_adaptativebeamInterpolation->getName()<<std::endl
					<<" MState1:"<<this->fromModel->getName()<< "  size:"<<this->fromModel->getSize()<<std::endl
					<<" MState2:"<<this->toModel->getName()<< "  size:"<<this->toModel->getSize()<<std::endl
					<<"idPointSubMap."<<idPointSubMap.size()<<std::endl
					<<"pointBeamDistribution."<<pointBeamDistribution.size()<<std::endl<<std::endl;
		}
	}
    virtual ~AdaptiveBeamMapping(){}
	
	void apply(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<VecCoord>& out, const Data<InVecCoord>& in);
	
	void applyJ(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<VecDeriv>& out, const Data<InVecDeriv>& in);
	
	void applyJT(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<InVecDeriv>& out, const Data<VecDeriv>& in);

    void applyJT(const core::ConstraintParams *cparams /* PARAMS FIRST */, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in);


    int addPoint ( const Coord& c, int /*indexFrom*/ )
    {

        int i = points.getValue().size();
        Vec3 test = c;

        points.beginEdit()->push_back(test);
        points.endEdit();
        return i;
    }

    void setBarycentricMapping()
    {
    	isBarycentricMapping=true;
    	points.beginEdit()->clear();points.endEdit();
    }
	int addBaryPoint(const int& _beamId,const Vec3& _baryCoord,bool /*todo_straightline_spline_option*/)
	//TODO add parameter label for different cases : unactive, linear, spline
	{
		//attention, beamId here is not the edge Id, but the id of a vec_edge_list defined in BeamInterpolation
		int newpointId = pointBeamDistribution.size();
		pointBeamDistribution.resize(newpointId+1);
		pointBeamDistribution[newpointId].baryPoint=_baryCoord;
		pointBeamDistribution[newpointId].beamId=_beamId;
		return newpointId;
	}
	//void clear(){};////// CTN_DEV todo for ContactMapper
	//clear the mapping in functions of size given
	void clear(int size)
	{
		this->clearidPointSubMap();
		pointBeamDistribution.clear();
		if ( size>0 && !isSubMapping)
		{
			pointBeamDistribution.reserve ( size );
	        points.beginEdit()->reserve ( size ); points.endEdit();
			this->getMechTo()[0]->resize(size);
		}
		else
		//case where this clear is call by a Multimapping, all component will be clear to null size
		{
	        points.beginEdit()->resize(0); points.endEdit();
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

    void clearidPointSubMap(){idPointSubMap.clear();}
    void addidPointSubMap(unsigned int _id){idPointSubMap.push_back(_id);}
    void setuseCurvAbs(bool _value){useCurvAbs.setValue(_value);}


protected:

    void applyJonPoint(unsigned int i, SpatialVector& VNode0input, SpatialVector& VNode1input, Deriv& vOutput, const InVecCoord& x);

    void applyJTonPoint(unsigned int i, const Deriv& finput, SpatialVector& FNode0output, SpatialVector& FNode1output, const InVecCoord& x);

    void computeJacobianOnPoint(unsigned int i, const typename In::VecCoord& x);

    void computeDistribution();

    sofa::core::topology::TopologyContainer* _topology;


    bool x_buf_used;
    typename In::VecCoord x_buf;

    sofa::helper::vector< PosPointDefinition > pointBeamDistribution;
    // for continuous_friction_contact:
    AdaptiveBeamMapping<TIn, TOut> *m_inputMapping;
    sofa::helper::vector<unsigned int> idPointSubMap;
    bool isSubMapping;
    bool isBarycentricMapping;

public :


};

} // namespace mapping

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_MAPPING_ADAPTIVEBEAMMAPPING_H */
