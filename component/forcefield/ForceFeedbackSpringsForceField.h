/***************************
* Initial software         *
* Authors: see Authors.txt *
* Copyright © Inria        *
* All rights reserved      *
* 2006-2018                *
* v1.0                     *
***************************/
#ifndef SOFA_COMPONENT_FORCEFIELD_FORCEFEEDBACKSPRINGFORCEFIELD_H
#define SOFA_COMPONENT_FORCEFIELD_FORCEFEEDBACKSPRINGFORCEFIELD_H


#include <sofa/core/behavior/ForceField.h>
#include <SofaDeformable/RestShapeSpringsForceField.h>
#include "../initBeamAdapter.h"

namespace sofa
{

namespace component
{

namespace forcefield
{

/**
* @brief This class describes a simple elastic springs ForceField between DOFs positions and rest positions.
*
* Springs are applied to given degrees of freedom between their current positions and their rest shape positions.
* An external MechanicalState reference can also be passed to the ForceField as rest shape position.
*/
template<class DataTypes>
class ForceFeedbackSpringsForceField : public sofa::component::forcefield::RestShapeSpringsForceField<DataTypes>
{
public:
    SOFA_CLASS(SOFA_TEMPLATE(ForceFeedbackSpringsForceField, DataTypes), SOFA_TEMPLATE(sofa::component::forcefield::RestShapeSpringsForceField, DataTypes));

    typedef core::behavior::ForceField<DataTypes> Inherit;
    typedef typename DataTypes::VecCoord VecCoord;
    typedef typename DataTypes::VecDeriv VecDeriv;
    typedef typename DataTypes::Coord Coord;
    typedef typename DataTypes::CPos CPos;
    typedef typename DataTypes::Deriv Deriv;
    typedef typename DataTypes::Real Real;
    typedef helper::vector< unsigned int > VecIndex;
    typedef helper::vector< Real >	 VecReal;

    typedef core::objectmodel::Data<VecCoord> DataVecCoord;
    typedef core::objectmodel::Data<VecDeriv> DataVecDeriv;

    Data< Deriv > d_forceFeedback;

protected:
    ForceFeedbackSpringsForceField();
public:

    /// init function
    void init() override;
    void bwdInit() override;

    /// Add the forces.
    void addForce(const core::MechanicalParams* mparams, DataVecDeriv& f, const DataVecCoord& x, const DataVecDeriv& v) override;

    bool local_useRestMState;
};


#if defined(SOFA_EXTERN_TEMPLATE) && !defined(SOFA_COMPONENT_FORCEFIELD_FORCEFEEDBACKSPRINGFORCEFIELD_CPP)
#ifndef SOFA_FLOAT
extern template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<sofa::defaulttype::Vec3dTypes>;
extern template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<sofa::defaulttype::Vec1dTypes>;
extern template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<sofa::defaulttype::Rigid3dTypes>;
#endif
#ifndef SOFA_DOUBLE
extern template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<sofa::defaulttype::Vec3fTypes>;
extern template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<sofa::defaulttype::Vec1fTypes>;
extern template class SOFA_BEAMADAPTER_API ForceFeedbackSpringsForceField<sofa::defaulttype::Rigid3fTypes>;
#endif

#endif

} // namespace forcefield

} // namespace component

} // namespace sofa

#endif // SOFA_COMPONENT_FORCEFIELD_FORCEFEEDBACKSPRINGFORCEFIELD_H
