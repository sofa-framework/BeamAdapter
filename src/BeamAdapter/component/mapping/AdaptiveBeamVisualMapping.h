#pragma once

#include <BeamAdapter/config.h>

#include <BeamAdapter/component/engine/WireRestShape.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/topology/TopologicalMapping.h>
#include <sofa/core/DataEngine.h>

#include <BeamAdapter/component/mapping/AdaptiveBeamMapping.h>
#include <BeamAdapter/component/engine/WireRestShape.h>
#include <BeamAdapter/component/WireBeamInterpolation.h>

#include <sofa/component/statecontainer/MappedObject.h>

namespace beamadapter
{

template<typename InputDataTypes, typename OutputDataTypes = sofa::defaulttype::Vec3Types>
class AdaptiveBeamVisualMapping : public Mapping<InputDataTypes, OutputDataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE2(AdaptiveBeamVisualMapping, InputDataTypes, OutputDataTypes), SOFA_TEMPLATE2(Mapping, InputDataTypes, OutputDataTypes));
    
    using Inherit = Mapping<InputDataTypes, OutputDataTypes>;
        
    using Index = sofa::Index;
    using Edge = sofa::core::topology::BaseMeshTopology::Edge;
    using SeqEdges = sofa::core::topology::BaseMeshTopology::SeqEdges;
    using Quad = sofa::core::topology::BaseMeshTopology::Quad;
    using SeqQuads = sofa::core::topology::BaseMeshTopology::SeqQuads;
    using VecIndex = sofa::type::vector<Index>;
    using InCoord = typename InputDataTypes::Coord;
    using InVecCoord = typename InputDataTypes::VecCoord;
    using InVecDeriv = typename InputDataTypes::VecDeriv;
    using InMatrixDeriv = typename OutputDataTypes::MatrixDeriv;
    using OutCoord = typename OutputDataTypes::Coord;
    using OutVecCoord = typename OutputDataTypes::VecCoord;
    using OutVecDeriv = typename OutputDataTypes::VecDeriv;
    using OutMatrixDeriv = typename OutputDataTypes::MatrixDeriv;
        
protected:
    AdaptiveBeamVisualMapping();
    ~AdaptiveBeamVisualMapping() override = default;
    
    SeqEdges m_internalEdges;
    typename sofa::component::statecontainer::MappedObject<InputDataTypes>::SPtr m_centerState; // state by itself should be sufficient...
    typename AdaptiveBeamMapping<InputDataTypes, InputDataTypes>::SPtr m_beamMapping;

public:
    void init() override;
    void bwdInit() override;
    void apply(const MechanicalParams* mparams, Data<OutVecCoord>& out, const Data<InVecCoord>& in) override;
    void applyJ(const MechanicalParams*, Data<OutVecDeriv>& , const Data<InVecDeriv>& ) override {};
    void applyJT(const MechanicalParams*, Data<InVecDeriv>& , const Data<OutVecDeriv>& ) override {};
    //void applyJT(const core::ConstraintParams*, Data<InMatrixDeriv>&, const Data<OutMatrixDeriv>&) override {};
        
    //inputs
    SingleLink<AdaptiveBeamVisualMapping<InputDataTypes, OutputDataTypes>, WireBeamInterpolation<InputDataTypes>, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_wireBeamInterpolation;
    SingleLink<AdaptiveBeamVisualMapping<InputDataTypes, OutputDataTypes>, sofa::core::topology::BaseMeshTopology, BaseLink::FLAG_STOREPATH|BaseLink::FLAG_STRONGLINK> l_quadTopology;
    
    // internal data
    Data<InVecCoord> d_localCenters;
    Data<OutVecCoord> d_circlePoints;
    
    // parameters
    Data<sofa::Size> d_nbPointsOnEachCircle; ///< Discretization of created circles
    Data<SReal> d_radius; ///< Radius of created circles in yz plan
    Data<SReal> d_thickness; ///< if not 0, add a second layer (where this value will be the distance between the two layers)
    ///<
    Data<bool> d_flipNormals; ///< Flip Normal ? (Inverse point order when creating quad)
};

#if !defined(BEAMADAPTER_COMPONENT_MAPPING_BEAM2QUADMAPPING_CPP)
extern template class SOFA_BEAMADAPTER_API Beam2QuadMapping<sofa::defaulttype::Rigid3Types>;
#endif

} // namespace beamadapter
