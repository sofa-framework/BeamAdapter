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
// C++ Implementation : UnifiedMultiMultiAdaptiveBeamMapping
//
// Description:
//
//
// Author: Christian Duriez, INRIA
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef SOFA_COMPONENT_MAPPING_MULTIADAPTIVEBEAMMAPPING_H
#define SOFA_COMPONENT_MAPPING_MULTIADAPTIVEBEAMMAPPING_H

#include "initBeamAdapter.h"
#include "WireBeamInterpolation.h"
#include "InterventionalRadiologyController.h"
#include "AdaptiveBeamMapping.h"
#include <sofa/core/Mapping.h>
#include <sofa/core/behavior/MechanicalState.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/helper/vector.h>
#include <sofa/defaulttype/Vec.h>
#include <sofa/defaulttype/Mat.h>
#include <SofaBaseTopology/EdgeSetTopologyModifier.h>
#include <sofa/core/objectmodel/Event.h>

using namespace sofa::component::controller;
using namespace sofa::component::fem;
using namespace sofa::core::objectmodel;

namespace sofa
{

namespace component
{

namespace mapping
{

/*!
 * \class MultiAdaptiveBeamMapping
 * @brief MultiAdaptiveBeamMapping Class
 */
template <class TIn, class TOut>
class MultiAdaptiveBeamMapping : public core::Mapping<TIn, TOut>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(MultiAdaptiveBeamMapping,TIn,TOut), SOFA_TEMPLATE2(core::Mapping,TIn,TOut));

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

    typedef helper::vector<unsigned int> VecIndex;
    typedef sofa::core::topology::BaseMeshTopology::EdgeID ElementID;
    typedef sofa::helper::vector<sofa::core::topology::BaseMeshTopology::Edge>   VecEdges    ;
    typedef sofa::helper::vector<sofa::core::topology::BaseMeshTopology::EdgeID> VecElementID;


    typedef typename sofa::defaulttype::SolidTypes<InReal>::Transform Transform;
    typedef std::pair<int, Transform> IndexedTransform;
    typedef typename sofa::defaulttype::SolidTypes< InReal>::SpatialVector SpatialVector;

    typedef sofa::defaulttype::Vec<3, Real> Vec3;
    typedef sofa::defaulttype::Vec<6, Real> Vec6;
    typedef sofa::defaulttype::Mat<3,12,Real> Mat3x12;
    typedef sofa::defaulttype::Mat<12,3,Real> Mat12x3;
    typedef sofa::defaulttype::Mat<6,12,Real> Mat6x12;
    typedef sofa::defaulttype::Mat<12,6,Real> Mat12x6;

    typedef std::pair<unsigned int, Vec3> BeamIdAndBaryCoord;
    typedef sofa::component::controller::InterventionalRadiologyController<TIn> TInterventionalRadiologyController;

public:

    Data<bool> useCurvAbs;
    Data< helper::vector< std::string > > m_controlerPath;

    MultiAdaptiveBeamMapping(core::State< In >* from, core::State< Out >* to,InterventionalRadiologyController<TIn>* _ircontroller);
    MultiAdaptiveBeamMapping();

    virtual ~MultiAdaptiveBeamMapping(){}

    void apply(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<VecCoord>& out, const Data<InVecCoord>& in);

    void applyJ(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<VecDeriv>& out, const Data<InVecDeriv>& in);

    void applyJT(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<InVecDeriv>& out, const Data<VecDeriv>& in);

    void applyJT(const core::ConstraintParams *cparams /* PARAMS FIRST */, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in);


    virtual void init();
    virtual void bwdInit();

    virtual void handleEvent(sofa::core::objectmodel::Event *);

    virtual void draw(const core::visual::VisualParams*);

    void setBarycentricMapping()
    {
        isBarycentricMapping=true;
        for(unsigned int i=0;i<m_subMappingList.size();i++)
        {
            m_subMappingList[i]->setBarycentricMapping();
        }
    }

    int addBaryPoint(const int& edgeId,const Vec3& _baryCoord,bool isStraight)
    //TODO add parameter label for different cases : unactive, linear, spline
    {
        int returnId=this->getMechTo()[0]->getSize();
        this->getMechTo()[0]->resize(returnId+1);

        /*
         * The idea is given an edgeId, how to find out what instrument is on this edge
         * in order to evaluate the addBaryPoint of the coresponding SubMapping
         *
         * --- --- --- --- --- --- --- B === id === === === A
         *                                |
         *                              edgeId
         *
         * The entire mesh denote the line above,
         * --- --- denote the edges which are not under control
         * --- --- denote the edges which are under control
         *
         * id_instrument_curvAbs_table contain informations about instruments id on nodes from A to B
         *
         *
         * */

        assert(m_ircontroller !=NULL && isBarycentricMapping);
        const vector<vector<int> >& id_instrument_curvAbs_table = m_ircontroller->get_id_instrument_curvAbs_table();
        int nbControlledEdge  = id_instrument_curvAbs_table.size() - 1;
        int totalNbEdges = m_ircontroller->getTotalNbEdges();
        int nbUnControlledEdges = totalNbEdges - nbControlledEdge;
        assert(nbUnControlledEdges>=0);

        if (edgeId < (totalNbEdges-nbControlledEdge) )
        {
            //if the edge in question is not under control, dont need to compute for collision
        }
        else
        {
            int controledEdgeId = edgeId-nbUnControlledEdges;
            const sofa::helper::vector<int>&  id_instrument_table_on_node = id_instrument_curvAbs_table[controledEdgeId+1];
            sofa::helper::vector< sofa::component::fem::WireBeamInterpolation<In>  *> m_instrumentsList;
            m_ircontroller->getInstrumentList(m_instrumentsList);
            Real radius = m_instrumentsList[id_instrument_table_on_node[0]]->getBeamSection(controledEdgeId)._r;
            int idInstrument  = id_instrument_table_on_node[0];
            for(unsigned int i=1;i<id_instrument_table_on_node.size();i++)
            {
                if(radius < m_instrumentsList[id_instrument_table_on_node[i]]->getBeamSection(controledEdgeId)._r)
                {
                    radius = m_instrumentsList[id_instrument_table_on_node[i]]->getBeamSection(controledEdgeId)._r;
                    idInstrument = id_instrument_table_on_node[i];
                }
            }
            //idPointMapp=m_subMappingList[idInstrument]->addBaryPoint(edgeId,_baryCoord,isStraight);
            m_subMappingList[idInstrument]->addBaryPoint(controledEdgeId,_baryCoord,isStraight);
            m_subMappingList[idInstrument]->addidPointSubMap(returnId);
        }
        return returnId;
    }

    void clear(int size)
    {
        for(unsigned int i=0;i<m_subMappingList.size();i++)
            m_subMappingList[i]->clear(size);
        this->getMechTo()[0]->resize(0);
        //this->getMechTo()[0]->resize(size);
    }


protected:

    sofa::core::topology::TopologyContainer* _topology;

    void assignSubMappingFromControllerInfo();

    // for fromSeveralInterpolations option
    sofa::helper::vector< sofa::component::fem::WireBeamInterpolation<TIn>  *> m_instrumentList;
    sofa::helper::vector<  AdaptiveBeamMapping<TIn, TOut>* > m_subMappingList;
    TInterventionalRadiologyController* m_ircontroller;
    sofa::component::topology::EdgeSetTopologyModifier* _edgeMod;
    sofa::helper::vector<InReal> _xPointList;     //=> for each mapped point provides the local position (curv. abs.)
    sofa::helper::vector<int> _idm_instrumentList; //=> for each mapped point provides the interpolation (in m_instrumentList)
    bool isBarycentricMapping;


};

} // namespace mapping

} // namespace component

} // namespace sofa

#endif  /* SOFA_COMPONENT_MAPPING_MULTIADAPTIVEBEAMMAPPING_H */
