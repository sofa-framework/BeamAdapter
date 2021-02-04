#include "ImplicitSurfaceAdaptiveConstraint.inl"

#include <sofa/defaulttype/Vec3Types.h>
#include <SofaBaseMechanics/MechanicalObject.h>
#include <sofa/core/ObjectFactory.h>

namespace sofa
{

namespace component
{

namespace constraint
{

using namespace sofa::defaulttype;
using namespace sofa::helper;
using core::RegisterObject;

#ifdef SOFAEVE
SOFA_DECL_CLASS(ImplicitSurfaceAdaptiveConstraint)

//TODO(damien): Il faut remplacer les descriptions dans RegisterObject par un vrai description
int ImplicitSurfaceAdaptiveConstraintClass = RegisterObject("PROUT TODO-ImplicitSurfaceAdaptiveConstraint")
.add< ImplicitSurfaceAdaptiveConstraint<Rigid3Types> >(true)

;

template class ImplicitSurfaceAdaptiveConstraint<Rigid3Types>;


#endif  // SOFAEVE

} // namespace constraint

} // namespace component

} // namespace sofa
