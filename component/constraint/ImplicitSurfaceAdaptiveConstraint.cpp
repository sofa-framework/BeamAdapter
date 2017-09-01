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
#ifdef SOFA_WITH_FLOAT
.add< ImplicitSurfaceAdaptiveConstraint<Rigid3fTypes> >()
#endif
#ifdef SOFA_WITH_DOUBLE
.add< ImplicitSurfaceAdaptiveConstraint<Rigid3dTypes> >(true)
#endif
;

#ifdef SOFA_WITH_FLOAT
template class ImplicitSurfaceAdaptiveConstraint<Rigid3fTypes>;
#endif
#ifdef SOFA_WITH_DOUBLE
template class ImplicitSurfaceAdaptiveConstraint<Rigid3dTypes>;
#endif

#endif  // SOFAEVE

} // namespace constraint

} // namespace component

} // namespace sofa
