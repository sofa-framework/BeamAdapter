#define BEAMADAPTER_COMPONENT_MAPPING_BEAM2QUADMAPPING_CPP

#include <BeamAdapter/component/mapping/AdaptiveBeamVisualMapping.inl>


#include <sofa/defaulttype/VecTypes.h>
#include <sofa/defaulttype/RigidTypes.h>
#include <sofa/core/ObjectFactory.h>

namespace beamadapter
{

void registerAdaptiveBeamVisualMapping(sofa::core::ObjectFactory* factory)
{
    factory->registerObjects(core::ObjectRegistrationData("Special mapping where quads are created from beams.")
             .add< AdaptiveBeamVisualMapping<sofa::defaulttype::Rigid3Types> >());
}

template class SOFA_BEAMADAPTER_API AdaptiveBeamVisualMapping<sofa::defaulttype::Rigid3Types>;

} // namespace beamadapter
