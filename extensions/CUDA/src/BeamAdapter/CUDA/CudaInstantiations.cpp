/******************************************************************************
*                 SOFA, Simulation Open-Framework Architecture                *
*                    (c) 2006 INRIA, USTL, UJF, CNRS, MGH                     *
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
* Authors: The SOFA Team and external contributors (see Authors.txt)          *
*                                                                             *
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <BeamAdapter/CUDA/config.h>
#include <sofa/gpu/cuda/CudaTypes.h>

#include <BeamAdapter/component/WireBeamInterpolation.inl>
#include <BeamAdapter/component/BeamInterpolation.inl>
#include <BeamAdapter/component/engine/WireRestShape.inl>
#include <BeamAdapter/component/forcefield/AdaptiveBeamForceFieldAndMass.inl>
#include <BeamAdapter/component/controller/InterventionalRadiologyController.inl>
#include <BeamAdapter/component/mapping/AdaptiveBeamMapping.inl>
#include <BeamAdapter/component/mapping/MultiAdaptiveBeamMapping.inl>
#include <BeamAdapter/component/model/RodMeshSection.inl>
#include <BeamAdapter/component/model/RodSpireSection.inl>
#include <BeamAdapter/component/model/RodStraightSection.inl>

#include <sofa/core/behavior/Mass.inl>
#include <sofa/core/Mapping.inl>
#include <sofa/component/controller/MechanicalStateController.inl>
#include <sofa/component/constraint/projective/FixedProjectiveConstraint.inl> // for InterventionalRadiologyController

#include <sofa/core/ObjectFactory.h>

using namespace sofa::gpu::cuda;

namespace beamadapter::cuda
{
    // template class SOFA_BEAMADAPTER_CUDA_API BeamInterpolation<CudaRigid3fTypes>;
    // template class SOFA_BEAMADAPTER_CUDA_API WireBeamInterpolation<CudaRigid3fTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API WireRestShape<CudaRigid3fTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API InterventionalRadiologyController<CudaRigid3fTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API AdaptiveBeamMapping<CudaRigid3fTypes, defaulttype::Vec3Types>;
    template class SOFA_BEAMADAPTER_CUDA_API MultiAdaptiveBeamMapping<CudaRigid3fTypes, defaulttype::Vec3Types>;
    template class SOFA_BEAMADAPTER_CUDA_API RodMeshSection<CudaRigid3fTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API RodSpireSection<CudaRigid3fTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API RodStraightSection<CudaRigid3fTypes>;

#ifdef SOFA_GPU_CUDA_DOUBLE
    template class SOFA_BEAMADAPTER_CUDA_API BeamInterpolation<CudaRigid3dTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API WireBeamInterpolation<CudaRigid3dTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API WireRestShape<CudaRigid3dTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API InterventionalRadiologyController<CudaRigid3dTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API AdaptiveBeamMapping<CudaRigid3dTypes, defaulttype::Vec3Types>;
    template class SOFA_BEAMADAPTER_CUDA_API MultiAdaptiveBeamMapping<CudaRigid3dTypes, defaulttype::Vec3Types>;
    template class SOFA_BEAMADAPTER_CUDA_API RodMeshSection<CudaRigid3dTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API RodSpireSection<CudaRigid3dTypes>;
    template class SOFA_BEAMADAPTER_CUDA_API RodStraightSection<CudaRigid3dTypes>;
#endif

using namespace sofa::gpu::cuda;

void registerBeamAdapterCUDAComponents(sofa::core::ObjectFactory* factory)
{
#ifdef SOFA_GPU_CUDA_DOUBLE
    factory->registerObjects(sofa::core::ObjectRegistrationData("Adaptive Beam Interpolation - Supports GPU-side computations using CUDA")
        // .add< BeamInterpolation<CudaRigid3fTypes> >()
        .add< BeamInterpolation<CudaRigid3dTypes> >());
#endif

#ifdef SOFA_GPU_CUDA_DOUBLE
    factory->registerObjects(sofa::core::ObjectRegistrationData("Adaptive Wire Beam Interpolation - Supports GPU-side computations using CUDA")
        // .add< WireBeamInterpolation<CudaRigid3fTypes> >()
        .add< WireBeamInterpolation<CudaRigid3dTypes> >());
#endif

    factory->registerObjects(sofa::core::ObjectRegistrationData("Wire Shape - Supports GPU-side computations using CUDA")
        .add< WireRestShape<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< WireRestShape<CudaRigid3dTypes> >()
#endif
    );

    factory->registerObjects(sofa::core::ObjectRegistrationData("Adaptive Beam finite elements - Supports GPU-side computations using CUDA")
        .add< AdaptiveBeamForceFieldAndMass<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< AdaptiveBeamForceFieldAndMass<CudaRigid3dTypes> >()
#endif
    );

    factory->registerObjects(sofa::core::ObjectRegistrationData("Provides a Mouse & Keyboard user control on an EdgeSet Topology - Supports GPU-side computations using CUDA")
        .add< InterventionalRadiologyController<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< InterventionalRadiologyController<CudaRigid3dTypes> >()
#endif
    );

    factory->registerObjects(sofa::core::ObjectRegistrationData("Provides a Mouse & Keyboard user control on an EdgeSet Topology - Supports GPU-side computations using CUDA")
            .add< InterventionalRadiologyController<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
            .add< InterventionalRadiologyController<CudaRigid3dTypes> >()
#endif
    );

    factory->registerObjects(sofa::core::ObjectRegistrationData("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs - Supports GPU-side computations using CUDA")
        .add< AdaptiveBeamMapping<CudaRigid3fTypes, defaulttype::Vec3Types> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< AdaptiveBeamMapping<CudaRigid3dTypes, defaulttype::Vec3Types> >()
#endif
    );

    factory->registerObjects(sofa::core::ObjectRegistrationData("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs - Supports GPU-side computations using CUDA")
        .add< MultiAdaptiveBeamMapping<CudaRigid3fTypes, defaulttype::Vec3Types> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< MultiAdaptiveBeamMapping<CudaRigid3dTypes, defaulttype::Vec3Types> >()
#endif
    );

    factory->registerObjects(sofa::core::ObjectRegistrationData("Class defining a Rod Section using a MeshLoader and material parameters using CUDA.")
        .add< RodMeshSection<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< RodMeshSection<CudaRigid3dTypes> >()
#endif
    );

    factory->registerObjects(sofa::core::ObjectRegistrationData("Class defining a rod spire section, defining material and geometry parameters using CUDA.")
        .add< RodSpireSection<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< RodSpireSection<CudaRigid3dTypes> >()
#endif
    );

    factory->registerObjects(sofa::core::ObjectRegistrationData("Class defining a rod straight section Material, defining material and geometry parameters using CUDA.")
        .add< RodStraightSection<CudaRigid3fTypes> >()
#ifdef SOFA_GPU_CUDA_DOUBLE
        .add< RodStraightSection<CudaRigid3dTypes> >()
#endif
    );
}

} // namespace beamadapter
