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
#include <SofaCUDA/sofa/gpu/cuda/CudaTypes.h>

#include <BeamAdapter/component/WireBeamInterpolation.inl>
#include <BeamAdapter/component/BeamInterpolation.inl>
#include <BeamAdapter/component/engine/WireRestShape.inl>
#include <BeamAdapter/component/forcefield/AdaptiveBeamForceFieldAndMass.inl>
#include <BeamAdapter/component/controller/InterventionalRadiologyController.inl>
#include <BeamAdapter/component/mapping/AdaptiveBeamMapping.inl>
// #include <BeamAdapter/component/mapping/MultiAdaptiveBeamMapping.inl>

#include <sofa/core/behavior/Mass.inl>
#include <sofa/core/Mapping.inl>
#include <sofa/component/controller/MechanicalStateController.inl>
#include <sofa/component/constraint/projective/FixedConstraint.inl> // for InterventionalRadiologyController

#include <sofa/core/ObjectFactory.h>

using namespace sofa::gpu::cuda;

namespace sofa::component::fem::_beaminterpolation_
{
    template class SOFA_BEAMADAPTER_API BeamInterpolation<CudaRigid3dTypes>;
} // namespace sofa::component::fem::_beaminterpolation_

namespace sofa::component::fem::_wirebeaminterpolation_
{
    template class SOFA_BEAMADAPTER_API WireBeamInterpolation<CudaRigid3dTypes>;
} // namespace sofa::component::fem::_beaminterpolation_

namespace sofa::component::engine::_wirerestshape_
{
    template class SOFA_BEAMADAPTER_API WireRestShape<CudaRigid3dTypes>;
} // namespace sofa::component::engine::_wirerestshape_

namespace sofa::component::controller::_interventionalradiologycontroller_
{
    template class SOFA_BEAMADAPTER_API InterventionalRadiologyController<CudaRigid3dTypes>;
} // namespace sofa::component::controller::_interventionalradiologycontroller_

namespace sofa::component::mapping::_adaptivebeammapping_
{
//    template class SOFA_BEAMADAPTER_API AdaptiveBeamMapping<CudaRigid3dTypes, defaulttype::Vec3Types>;
} // namespace sofa::component::mapping::_adaptivebeammapping_

namespace sofa::component::mapping
{
    // template class SOFA_BEAMADAPTER_API MultiAdaptiveBeamMapping<CudaRigid3dTypes, defaulttype::Vec3Types>;
} // namespace sofa::component::mapping


namespace sofa::gpu::cuda
{

int CudaBeamInterpolationClass = core::RegisterObject("Adaptive Beam Interpolation")
    .add< sofa::component::fem::BeamInterpolation<CudaRigid3dTypes> >();

int CudaWireBeamInterpolationClass = core::RegisterObject("Adaptive Wire Beam Interpolation")
    .add< sofa::component::fem::WireBeamInterpolation<CudaRigid3dTypes> >();

int CudaWireRestShapenClass = core::RegisterObject("Wire Shape")
    .add< sofa::component::engine::WireRestShape<CudaRigid3dTypes> >();

int CudaAdaptiveBeamForceFieldAndMassClass = core::RegisterObject("Adaptive Beam finite elements")
    .add< sofa::component::forcefield::AdaptiveBeamForceFieldAndMass<CudaRigid3dTypes> >();

int CudaInterventionalRadiologyControllerClass = core::RegisterObject("Provides a Mouse & Keyboard user control on an EdgeSet Topology.")
    .add< sofa::component::controller::InterventionalRadiologyController<CudaRigid3dTypes> >();

//int CudaAdaptiveBeamMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs")
//    .add< sofa::component::mapping::AdaptiveBeamMapping<CudaRigid3dTypes, defaulttype::Vec3Types> >();

//int CudaMultiAdaptiveBeamMappingClass = core::RegisterObject("Set the positions and velocities of points attached to a beam using linear interpolation between DOFs")
//.add< sofa::component::mapping::MultiAdaptiveBeamMapping<CudaRigid3dTypes, defaulttype::Vec3Types> >();


} // namespace sofa::gpu::cuda
