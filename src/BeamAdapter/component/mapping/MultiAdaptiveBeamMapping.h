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
#pragma once

#include <BeamAdapter/config.h>
#include <BeamAdapter/component/controller/InterventionalRadiologyController.h>
#include <BeamAdapter/component/mapping/AdaptiveBeamMapping.h>

namespace sofa::component::mapping
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

    typedef type::vector<unsigned int> VecIndex;
    typedef sofa::core::topology::BaseMeshTopology::EdgeID ElementID;
    typedef sofa::type::vector<sofa::core::topology::BaseMeshTopology::Edge>   VecEdges    ;
    typedef sofa::type::vector<sofa::core::topology::BaseMeshTopology::EdgeID> VecElementID;


    typedef typename sofa::defaulttype::SolidTypes<InReal>::Transform Transform;
    typedef std::pair<int, Transform> IndexedTransform;
    typedef typename sofa::defaulttype::SolidTypes< InReal>::SpatialVector SpatialVector;

    typedef sofa::type::Vec<3, Real> Vec3;
    typedef sofa::type::Vec<6, Real> Vec6;
    typedef sofa::type::Mat<3,12,Real> Mat3x12;
    typedef sofa::type::Mat<12,3,Real> Mat12x3;
    typedef sofa::type::Mat<6,12,Real> Mat6x12;
    typedef sofa::type::Mat<12,6,Real> Mat12x6;

    typedef std::pair<unsigned int, Vec3> BeamIdAndBaryCoord;
    typedef sofa::component::controller::InterventionalRadiologyController<TIn> TInterventionalRadiologyController;

public:
    Data<bool> useCurvAbs;
    Data< type::vector< std::string > > m_controlerPath;
    Data<bool> d_parallelMapping;           /*!< flag to enable parallel internal computation of apply/applyJ for the submapping(s) AdaptiveBeamMapping */

    MultiAdaptiveBeamMapping(core::State< In >* from, core::State< Out >* to, TInterventionalRadiologyController* _ircontroller);
    MultiAdaptiveBeamMapping();

    virtual ~MultiAdaptiveBeamMapping() = default;

    void apply(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<VecCoord>& out, const Data<InVecCoord>& in) override;

    void applyJ(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<VecDeriv>& out, const Data<InVecDeriv>& in) override;

    void applyJT(const core::MechanicalParams *mparams /* PARAMS FIRST */, Data<InVecDeriv>& out, const Data<VecDeriv>& in) override;

    void applyJT(const core::ConstraintParams *cparams /* PARAMS FIRST */, Data<InMatrixDeriv>& out, const Data<OutMatrixDeriv>& in) override;


    virtual void init() override;
    virtual void bwdInit() override;

    virtual void handleEvent(sofa::core::objectmodel::Event *) override;

    void setBarycentricMapping();

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
     * TODO add parameter label for different cases : unactive, linear, spline
     * */
    int addBaryPoint(const int& edgeId,const Vec3& _baryCoord,bool isStraight);

    void clear(int size);


protected:

    sofa::core::topology::TopologyContainer* _topology{ nullptr };

    void assignSubMappingFromControllerInfo();

    // for fromSeveralInterpolations option
    sofa::type::vector< sofa::component::fem::WireBeamInterpolation<TIn>  *> m_instrumentList;
    sofa::type::vector<  typename AdaptiveBeamMapping<TIn, TOut>::SPtr > m_subMappingList;
    TInterventionalRadiologyController* m_ircontroller;
    sofa::component::topology::container::dynamic::EdgeSetTopologyModifier* _edgeMod{nullptr};
    sofa::type::vector<InReal> _xPointList;     //=> for each mapped point provides the local position (curv. abs.)
    sofa::type::vector<int> _idm_instrumentList; //=> for each mapped point provides the interpolation (in m_instrumentList)
    bool isBarycentricMapping;


};

#if !defined(SOFA_PLUGIN_BEAMADAPTER_MULTIADAPTIVEBEAMMAPPING_CPP)
extern template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamMapping<defaulttype::Rigid3Types, defaulttype::Vec3Types>;
#endif

} // namespace sofa::component::mapping
