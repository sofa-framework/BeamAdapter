#pragma once

#include <BeamAdapter/component/mapping/AdaptiveBeamVisualMapping.h>

#include <BeamAdapter/component/engine/WireRestShape.h>
#include <sofa/core/topology/BaseMeshTopology.h>
#include <sofa/core/behavior/MechanicalState.h>

#include <sofa/component/topology/container/dynamic/QuadSetTopologyModifier.h>
#include <sofa/component/topology/container/dynamic/QuadSetTopologyContainer.h>

#include <sofa/helper/visual/DrawTool.h>
#include <sofa/helper/ScopedAdvancedTimer.h>

namespace beamadapter
{

using Vec3 = sofa::type::Vec3;

template<typename InputDataTypes, typename OutputDataTypes>
AdaptiveBeamVisualMapping<InputDataTypes, OutputDataTypes>::AdaptiveBeamVisualMapping()
    : Inherit()
    , m_beamMapping()
    , l_wireBeamInterpolation(initLink("interpolation", "Path to the Interpolation component on scene"))
    , l_quadTopology(initLink("quadTopology", "Path to the Quad topology"))
    , d_nbPointsOnEachCircle( initData(&d_nbPointsOnEachCircle, "nbPointsOnEachCircle", "Discretization of created circles"))
    , d_radius( initData(&d_radius, 1_sreal, "radius", "Radius of created circles in yz plan"))
    , d_thickness( initData(&d_thickness, 0_sreal, "thickness", "if not 0, add a second layer (where this value will be the distance between the two layers)"))
    , d_flipNormals(initData(&d_flipNormals, bool(false), "flipNormals", "Flip Normal ? (Inverse point order when creating quad)"))
{
}

template<typename InputDataTypes, typename OutputDataTypes>
void AdaptiveBeamVisualMapping<InputDataTypes, OutputDataTypes>::init()
{
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Invalid);
        
    if (!d_radius.isSet())
    {
        msg_error() << "No radius  defined";
        return;
    }

    if (d_radius.isSet() && d_radius.getValue() < std::numeric_limits<SReal>::min())
    {
        msg_error() << "Radius is zero or negative";
        return;
    }
    
    if (!d_nbPointsOnEachCircle.isSet())
    {
        msg_error() << "nbPointsOnEachCircle is not defined";
        return;
    }
        
    if (!l_wireBeamInterpolation)
    {
        msg_error() <<"No Beam Interpolation found, the component can not work.";
        return;
    }
    
    if (!l_quadTopology)
    {
        msg_error() <<"No Quad Topology found, the component can not work.";
        return;
    }
    
    if(d_radius.getValue() < d_thickness.getValue())
    {
        msg_error() << "Thickness ( " << d_thickness.getValue() << " ) cannot be bigger that the radius (" << d_radius.getValue() << " ).";
        return;
    }
  
    
    auto* wireShape = l_wireBeamInterpolation->m_restShape.get();
    assert(wireShape);
        
    // create internal topology at rest
    m_centerState = sofa::core::objectmodel::New<sofa::component::statecontainer::MappedObject<InputDataTypes>>();
    this->addSlave(m_centerState);
    
    auto localCenters = sofa::helper::getWriteOnlyAccessor(d_localCenters);
    localCenters.clear();
    m_internalEdges.clear();
    
    SReal prev_length = 0.0;
    sofa::Index prev_edges = 0;
    sofa::Index startPtId = 0;
    const auto& keyPts = wireShape->d_keyPoints.getValue();
    for (sofa::Size i = 0; i < wireShape->l_sectionMaterials.size(); ++i)
    {
        // Add topology of the material
        sofa::Size nbrVisuEdges = wireShape->l_sectionMaterials.get(i)->getNbVisualEdges();
        SReal length = fabs(keyPts[i + 1] - keyPts[i]);
        SReal dx = length / nbrVisuEdges;

        // add points from the material
        for (sofa::Index i = startPtId; i < nbrVisuEdges + 1; i++)
        {
            const InCoord p({prev_length + i * dx, 0.0, 0.0}, {0.0, 0.0, 0.0, 1.0});
            localCenters.push_back(p);
        }

        // add segments from the material
        for (sofa::Index i = prev_edges; i < prev_edges + nbrVisuEdges; i++)
        {
            m_internalEdges.emplace_back(i, i + 1);
        }

        prev_length = length;
        prev_edges = nbrVisuEdges;
        startPtId = 1; // Assume the last point of mat[n] == first point of mat[n+1]
    }
    
    auto mappedCenters = m_centerState->writeOnlyPositions();
    mappedCenters.wref() = localCenters; // initialize state
    
    // Cache the created vertices
    constexpr Vec3 Y0 {0.0_sreal, 1.0_sreal, 0.0_sreal};
    constexpr Vec3 Z0 {0.0_sreal, 0.0_sreal, 1.0_sreal};
    const auto N = d_nbPointsOnEachCircle.getValue();
    
    const auto localCentersSize = static_cast<sofa::Size>(localCenters.size());
    auto circlePoints = sofa::helper::getWriteOnlyAccessor(d_circlePoints);
    const auto thickness = d_thickness.getValue();
    const SReal rho = d_radius.getValue();
    if(thickness == 0.0)
    {
        circlePoints.resize(N * localCentersSize);
        for (sofa::Index i=0; i< localCentersSize; ++i)
        {
            const auto p0 = i;
            
            for(unsigned int j=0; j<N; ++j)
            {
                circlePoints[p0*N+j] = (Y0*cos((SReal) (2.0*j*M_PI/N)) + Z0*sin((SReal) (2.0*j*M_PI/N)))*((SReal) rho);
            }
        }
    }
    else
    {
        circlePoints.resize(N * 2 * localCentersSize);
        const SReal rho2 = d_radius.getValue() - thickness;
        for (sofa::Index i=0; i< localCentersSize; ++i)
        {
            const auto p0 = i;
            
            for(unsigned int j=0; j<N; ++j)
            {
                circlePoints[p0*2*N+j] = (Y0*cos((SReal) (2.0*j*M_PI/N)) + Z0*sin((SReal) (2.0*j*M_PI/N)))*((SReal) rho);
            }
            for(unsigned int j=0; j<N; ++j)
            {
                circlePoints[p0*2*N+j + N] = (Y0*cos((SReal) (2.0*j*M_PI/N)) + Z0*sin((SReal) (2.0*j*M_PI/N)))*((SReal) rho2);
            }
        }
    }
         
    m_beamMapping = sofa::core::objectmodel::New<AdaptiveBeamMapping<InputDataTypes, InputDataTypes>>(this->fromModel, m_centerState.get(), l_wireBeamInterpolation.get());
    m_beamMapping->init();
    this->addSlave(m_beamMapping);
    
    Inherit::init();
    
    this->d_componentState.setValue(sofa::core::objectmodel::ComponentState::Valid);
    
}


template<typename InputDataTypes, typename OutputDataTypes>
void AdaptiveBeamVisualMapping<InputDataTypes, OutputDataTypes>::bwdInit()
{
    if (!this->isComponentStateValid())
        return;
    
    m_beamMapping->bwdInit();
    
    const auto N = d_nbPointsOnEachCircle.getValue();
    const auto& edges = m_internalEdges;
    const bool flipNormals = d_flipNormals.getValue();
    
    if(l_quadTopology->getNbQuads() == 0)
    {
        if(d_thickness.getValue() == 0.0)
        {
            for(const auto& edge : edges)
            {
                const Index p0 = edge[0];
                const Index p1 = edge[1];
                
                for(unsigned int j=0; j<N; ++j)
                {
                    const Index q0 = p0*N+j;
                    const Index q1 = p1*N+j;
                    const Index q2 = p1*N+((j+1)%N);
                    const Index q3 = p0*N+((j+1)%N);
                    
                    if (flipNormals)
                    {
                        l_quadTopology->addQuad(q3, q2, q1, q0);
                    }
                    else
                    {
                        l_quadTopology->addQuad(q0, q1, q2, q3);
                    }
                }
            }
        }
        else
        {
            for(const auto& edge : edges)
            {
                const Index p0 = edge[0];
                const Index p1 = edge[1];
                
                for(unsigned int j=0; j<N; ++j)
                {
                    const Index q0 = p0*2*N+j;
                    const Index q1 = p1*2*N+j;
                    const Index q2 = p1*2*N+((j+1)%N);
                    const Index q3 = p0*2*N+((j+1)%N);
                    
                    const Index q4 = p0*2*N+j + N;
                    const Index q5 = p1*2*N+j + N;
                    const Index q6 = p1*2*N+((j+1)%N) + N;
                    const Index q7 = p0*2*N+((j+1)%N) + N;
                    
                    if (flipNormals)
                    {
                        l_quadTopology->addQuad(q3, q2, q1, q0);
                        l_quadTopology->addQuad(q4, q5, q6, q7); // normal of the inner layer must be inverted from the outer layer
                    }
                    else
                    {
                        l_quadTopology->addQuad(q0, q1, q2, q3);
                        l_quadTopology->addQuad(q7, q6, q5, q4);
                    }
                }
            }
            
            
            // closing the extremities
            const sofa::Size offset = static_cast<sofa::Size>(edges.size() * N * 2);
            if(flipNormals)
            {
                for(unsigned int j=0; j<N-1; ++j)
                {
                    l_quadTopology->addQuad(j+N, j+N+1, j+1, j);
                    l_quadTopology->addQuad(offset + j, offset+j+1, offset+j+N+1, offset+j+N);
                }
                l_quadTopology->addQuad(N+N-1, N, 0, N-1);
                l_quadTopology->addQuad(offset+N-1, offset+0, offset+N, offset+N+N-1);
            }
            else
            {
                for(unsigned int j=0; j<N-1; ++j)
                {
                    l_quadTopology->addQuad(j, j+1, j+N+1, j+N);
                    l_quadTopology->addQuad(offset+j+N, offset+j+N+1, offset+j+1, offset + j);
                }
                l_quadTopology->addQuad(N-1, 0, N, N+N-1);
                l_quadTopology->addQuad(offset+N+N-1, offset+N, offset+0, offset+N-1);
            }
        }
    }
}

template<typename InputDataTypes, typename OutputDataTypes>
void AdaptiveBeamVisualMapping<InputDataTypes, OutputDataTypes>::apply(const MechanicalParams* mparams, Data<OutVecCoord>& dOut, const Data<InVecCoord>& dIn)
{
    SCOPED_TIMER("AdaptiveBeamVisualMapping_Apply");
    
    if (!this->isComponentStateValid())
        return;
    
    const auto* beamState = this->l_wireBeamInterpolation->getMState();
    
    // workaround for l_wireBeamInterpolation->m_mstate is nullptr until bwdinit
    if (!beamState)
    {
        return;
    }
            
    // final quads positions
    auto out = sofa::helper::getWriteOnlyAccessor(dOut);
    
    // 1 - apply the mapping on the circle "centers" from the beams, using the classical AdaptiveBeamMapping
    auto toBeMappedCenters = m_centerState->writeOnlyPositions();
    auto localCenters = sofa::helper::getReadAccessor(d_localCenters);
    toBeMappedCenters.wref() = localCenters; // initialize state
    
    m_beamMapping->apply(mparams, *m_centerState->write(core::vec_id::write_access::position), dIn);
    
    // 2 - create the points around the centers
    
    auto N = d_nbPointsOnEachCircle.getValue();
    
    const auto& mappedCenters = m_centerState->readPositions();
        
    out.clear();
    
    auto circlePoints = sofa::helper::getReadAccessor(d_circlePoints);
    
    if(d_thickness.getValue() == 0.0)
    {
        out.resize(mappedCenters.size() * N);
        for (sofa::Index i=0; i< mappedCenters.size(); ++i)
        {
            const auto p0 = i;

            for(unsigned int j=0; j<N; ++j)
            {
                out[p0*N+j] = mappedCenters[p0].projectPoint(circlePoints[p0*N+j]);
            }
        }
    }
    else
    {
        out.resize(mappedCenters.size() * N * 2);
        for (sofa::Index i=0; i< mappedCenters.size(); ++i)
        {
            const auto p0 = i;

            for(unsigned int j=0; j<N; ++j)
            {
                out[p0*2*N+j] = mappedCenters[p0].projectPoint(circlePoints[p0*2*N+j]);
                out[p0*2*N+j + N] = mappedCenters[p0].projectPoint(circlePoints[p0*2*N+j + N]);
            }
        }
    }
    
    
    return;
}

} // namespace beamadapter
